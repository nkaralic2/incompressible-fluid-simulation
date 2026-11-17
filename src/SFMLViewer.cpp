#include "FLIPSimulator.hpp"

#include <SFML/Graphics.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#ifndef WATER_SHADER_PATH
#define WATER_SHADER_PATH "shaders/water.frag"
#endif

namespace {

struct Viewport {
    float scale;
    sf::Vector2f origin;
};

Viewport makeViewport(const SimConfig& config, sf::Vector2u windowSize) {
    const float padding = 48.0f;
    const float availableWidth = static_cast<float>(windowSize.x) - 2.0f * padding;
    const float availableHeight = static_cast<float>(windowSize.y) - 2.0f * padding;
    const float scale = std::min(availableWidth / config.domainWidth,
                                 availableHeight / config.domainHeight());
    const float contentWidth = scale * config.domainWidth;
    const float contentHeight = scale * config.domainHeight();
    return {scale, {0.5f * (static_cast<float>(windowSize.x) - contentWidth),
                    0.5f * (static_cast<float>(windowSize.y) - contentHeight)}};
}

sf::Vector2f toScreen(const Particle& particle, const SimConfig& config,
                      const Viewport& viewport) {
    return {viewport.origin.x + particle.x * viewport.scale,
            viewport.origin.y + (config.domainHeight() - particle.y) * viewport.scale};
}

sf::Color particleColor(const Particle& particle) {
    const float speed = std::sqrt(particle.vx * particle.vx + particle.vy * particle.vy);
    const float intensity = std::min(1.0f, speed / 8.0f);
    const auto red = static_cast<std::uint8_t>(80.0f + 175.0f * intensity);
    const auto blue = static_cast<std::uint8_t>(235.0f - 130.0f * intensity);
    return {red, 150, blue, 235};
}

sf::Vertex makeVertex(sf::Vector2f position, sf::Color color) {
    sf::Vertex vertex;
    vertex.position = position;
    vertex.color = color;
    return vertex;
}

class SmoothSurfaceRenderer {
public:
    void build(const FLIPSimulator& simulator, const Viewport& viewport,
         sf::VertexArray& triangles, sf::VertexArray& highlights,
          sf::Image& densityImage, sf::Image& velocityImage) const {
        const SimConfig& config = simulator.config;
        const int samplesX = config.gridNX * 2;
        const int samplesY = config.gridNY * 2;
        const int rowSize = samplesX + 1;
        std::vector<float> density(static_cast<std::size_t>(rowSize) * (samplesY + 1), 0.0f);
        const float supportRadius = 2.6f;
     const float supportWorldX = supportRadius / samplesX * config.domainWidth;
     const float supportWorldY = supportRadius / samplesY * config.domainHeight();

        auto densityAt = [&](int x, int y) -> float& {
            return density[static_cast<std::size_t>(y) * rowSize + x];
        };

        auto splat = [&](float particleX, float particleY) {
            const float sampleX = particleX / config.domainWidth * samplesX;
            const float sampleY = particleY / config.domainHeight() * samplesY;
            const int minX = std::max(0, static_cast<int>(std::floor(sampleX - supportRadius)));
            const int maxX = std::min(samplesX, static_cast<int>(std::ceil(sampleX + supportRadius)));
            const int minY = std::max(0, static_cast<int>(std::floor(sampleY - supportRadius)));
            const int maxY = std::min(samplesY, static_cast<int>(std::ceil(sampleY + supportRadius)));

            for (int y = minY; y <= maxY; ++y) {
                for (int x = minX; x <= maxX; ++x) {
                    const float dx = static_cast<float>(x) - sampleX;
                    const float dy = static_cast<float>(y) - sampleY;
                    const float distance = std::sqrt(dx * dx + dy * dy);
                    if (distance < supportRadius) {
                        densityAt(x, y) += 1.0f - distance / supportRadius;
                    }
                }
            }
        };

        for (const Particle& particle : simulator.particles) {
            std::array<float, 3> mirroredX = {particle.x, particle.x, particle.x};
            std::array<float, 3> mirroredY = {particle.y, particle.y, particle.y};
            int xCount = 1;
            int yCount = 1;
            if (particle.x < supportWorldX) {
                mirroredX[xCount++] = -particle.x;
            }
            if (config.domainWidth - particle.x < supportWorldX) {
                mirroredX[xCount++] = 2.0f * config.domainWidth - particle.x;
            }
            if (particle.y < supportWorldY) {
                mirroredY[yCount++] = -particle.y;
            }
            if (config.domainHeight() - particle.y < supportWorldY) {
                mirroredY[yCount++] = 2.0f * config.domainHeight() - particle.y;
            }
            for (int x = 0; x < xCount; ++x) {
                for (int y = 0; y < yCount; ++y) {
                    splat(mirroredX[x], mirroredY[y]);
                }
            }
        }

        const float normalization = std::max(1, config.particlesPerCell);
        for (float& value : density) {
            value /= normalization;
        }

        // Build textures representing actual particle density and local MAC velocities.
        // Images have their origin at the top-left; the shader uses the matching mapping.
        densityImage = sf::Image({static_cast<unsigned int>(samplesX + 1),
                                  static_cast<unsigned int>(samplesY + 1)}, sf::Color::Black);
        velocityImage = sf::Image({static_cast<unsigned int>(samplesX + 1),
                                   static_cast<unsigned int>(samplesY + 1)}, sf::Color(128, 128, 0));
        const auto encodeVelocity = [](float v) -> std::uint8_t {
            return static_cast<std::uint8_t>(std::round(
                std::clamp(0.5f + v / 16.0f, 0.0f, 1.0f) * 255.0f));
        };
        for (int y = 0; y <= samplesY; ++y) {
            for (int x = 0; x <= samplesX; ++x) {
                const float value = std::clamp(densityAt(x, y) / 0.8f, 0.0f, 1.0f);
                const auto byte = static_cast<std::uint8_t>(std::round(value * 255.0f));
                const sf::Vector2u imagePos(static_cast<unsigned int>(x),
                                             static_cast<unsigned int>(samplesY - y));
                densityImage.setPixel(imagePos, sf::Color(byte, byte, byte));
                const int i = std::clamp(x / 2, 0, simulator.grid.NX - 1);
                const int j = std::clamp(y / 2, 0, simulator.grid.NY - 1);
                const Grid& grid = simulator.grid;
                const float ux = 0.5f * (grid.u[grid.uIndex(i, j)] +
                                         grid.u[grid.uIndex(i + 1, j)]);
                const float vy = 0.5f * (grid.v[grid.vIndex(i, j)] +
                                         grid.v[grid.vIndex(i, j + 1)]);
                velocityImage.setPixel(imagePos, sf::Color(encodeVelocity(ux),
                                                           encodeVelocity(vy), 0));
            }
        }

        const float iso = 0.18f;
        std::vector<std::uint8_t> fluidMask(
            static_cast<std::size_t>(samplesX) * samplesY, 0);
        for (int y = 0; y < samplesY; ++y) {
            for (int x = 0; x < samplesX; ++x) {
                const float averageDensity = 0.25f * (
                    densityAt(x, y) + densityAt(x + 1, y) +
                    densityAt(x + 1, y + 1) + densityAt(x, y + 1));
                fluidMask[static_cast<std::size_t>(y) * samplesX + x] =
                    averageDensity >= iso ? 1 : 0;
            }
        }
        removeSmallComponents(fluidMask, samplesX, samplesY, 24);
        morphologicalCleanup(fluidMask, samplesX, samplesY);

        triangles = sf::VertexArray(sf::PrimitiveType::Triangles);
        highlights = sf::VertexArray(sf::PrimitiveType::Lines);

        // Clip each cell triangle against the SAME continuous iso-value.
        // Unlike the previous averaged-cell mask, this avoids disagreements
        // between the displayed mask and the reconstructed boundary.
        auto clipAndEmit = [&](const std::array<sf::Vector2f, 3>& positions,
                               const std::array<float, 3>& values,
                               sf::Color color) {
            std::array<sf::Vector2f, 6> poly{};
            int count = 0;
            for (int edge = 0; edge < 3; ++edge) {
                const int next = (edge + 1) % 3;
                const bool inside = values[edge] >= iso;
                const bool nextInside = values[next] >= iso;
                if (inside) poly[count++] = positions[edge];
                if (inside != nextInside) {
                    const float t = (iso - values[edge]) /
                                    (values[next] - values[edge]);
                    poly[count++] = positions[edge] + t *
                                    (positions[next] - positions[edge]);
                }
            }
            for (int k = 1; k + 1 < count; ++k) {
                triangles.append(makeVertex(poly[0], color));
                triangles.append(makeVertex(poly[k], color));
                triangles.append(makeVertex(poly[k + 1], color));
            }
        };

        for (int y = 0; y < samplesY; ++y) {
            for (int x = 0; x < samplesX; ++x) {
                if (fluidMask[static_cast<std::size_t>(y) * samplesX + x] == 0) {
                    continue;
                }
                const std::array<float, 4> values = {
                    densityAt(x, y), densityAt(x + 1, y),
                    densityAt(x + 1, y + 1), densityAt(x, y + 1)};
                const std::array<sf::Vector2f, 4> positions = {
                    samplePosition(x, y, samplesX, samplesY, config, viewport),
                    samplePosition(x + 1, y, samplesX, samplesY, config, viewport),
                    samplePosition(x + 1, y + 1, samplesX, samplesY, config, viewport),
                    samplePosition(x, y + 1, samplesX, samplesY, config, viewport)};
                const float densityValue = 0.25f * (values[0] + values[1] + values[2] + values[3]);
                const float heightFraction = (static_cast<float>(y) + 0.5f) / samplesY;
                const sf::Color fill = liquidColor(1.0f - heightFraction, densityValue,
                                                   0.0f, 0.0f);
                clipAndEmit({positions[0], positions[1], positions[2]},
                            {values[0], values[1], values[2]}, fill);
                clipAndEmit({positions[0], positions[2], positions[3]},
                            {values[0], values[2], values[3]}, fill);
            }
        }

    }

private:
    static void removeSmallComponents(std::vector<std::uint8_t>& mask,
                                      int width, int height, int minimumCells) {
        std::vector<std::uint8_t> visited(mask.size(), 0);
        const int offsets[8][2] = {
            {-1, -1}, {0, -1}, {1, -1}, {-1, 0},
            {1, 0}, {-1, 1}, {0, 1}, {1, 1}};

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const std::size_t start = static_cast<std::size_t>(y) * width + x;
                if (mask[start] == 0 || visited[start] != 0) {
                    continue;
                }

                std::queue<std::size_t> pending;
                std::vector<std::size_t> component;
                pending.push(start);
                visited[start] = 1;

                while (!pending.empty()) {
                    const std::size_t current = pending.front();
                    pending.pop();
                    component.push_back(current);
                    const int currentX = static_cast<int>(current % width);
                    const int currentY = static_cast<int>(current / width);

                    for (const auto& offset : offsets) {
                        const int neighborX = currentX + offset[0];
                        const int neighborY = currentY + offset[1];
                        if (neighborX < 0 || neighborX >= width ||
                            neighborY < 0 || neighborY >= height) {
                            continue;
                        }
                        const std::size_t neighbor =
                            static_cast<std::size_t>(neighborY) * width + neighborX;
                        if (mask[neighbor] != 0 && visited[neighbor] == 0) {
                            visited[neighbor] = 1;
                            pending.push(neighbor);
                        }
                    }
                }

                if (static_cast<int>(component.size()) < minimumCells) {
                    for (const std::size_t index : component) {
                        mask[index] = 0;
                    }
                }
            }
        }
    }

    static void morphologicalCleanup(std::vector<std::uint8_t>& mask,
                                     int width, int height) {
        std::vector<std::uint8_t> cleaned = mask;
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                const std::size_t index = static_cast<std::size_t>(y) * width + x;
                int neighbors = 0;
                for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                    for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                        if (offsetX == 0 && offsetY == 0) {
                            continue;
                        }
                        neighbors += mask[static_cast<std::size_t>(y + offsetY) * width +
                                          x + offsetX] != 0 ? 1 : 0;
                    }
                }

                if (mask[index] == 0 && neighbors >= 6) {
                    cleaned[index] = 1;
                } else if (mask[index] != 0 && neighbors <= 1) {
                    cleaned[index] = 0;
                }
            }
        }
        mask.swap(cleaned);
    }

    static float cellSpeed(const FLIPSimulator& simulator, int sampleX, int sampleY) {
        const Grid& grid = simulator.grid;
        const int i = std::clamp(sampleX / 2, 0, grid.NX - 1);
        const int j = std::clamp(sampleY / 2, 0, grid.NY - 1);
        const float u = 0.5f * (grid.u[grid.uIndex(i, j)] + grid.u[grid.uIndex(i + 1, j)]);
        const float v = 0.5f * (grid.v[grid.vIndex(i, j)] + grid.v[grid.vIndex(i, j + 1)]);
        return std::sqrt(u * u + v * v);
    }

    static sf::Color liquidColor(float depth, float density, float surface, float speed) {
        const float densityLight = std::clamp((density - 0.18f) / 1.2f, 0.0f, 1.0f);
        const float topLight = 1.0f - depth;
        const float highlight = surface * (0.25f + 0.5f * std::clamp(speed / 5.0f, 0.0f, 1.0f));
        const auto red = static_cast<std::uint8_t>(10.0f + 22.0f * topLight + 70.0f * highlight);
        const auto green = static_cast<std::uint8_t>(48.0f + 74.0f * topLight + 105.0f * highlight);
        const auto blue = static_cast<std::uint8_t>(120.0f + 76.0f * densityLight + 52.0f * topLight);
        const auto alpha = static_cast<std::uint8_t>(205.0f + 35.0f * densityLight);
        return {red, green, blue, alpha};
    }

    static sf::Vector2f samplePosition(int x, int y, int samplesX, int samplesY,
                                       const SimConfig& config, const Viewport& viewport) {
        const float worldX = static_cast<float>(x) / samplesX * config.domainWidth;
        const float worldY = static_cast<float>(y) / samplesY * config.domainHeight();
        return {viewport.origin.x + worldX * viewport.scale,
                viewport.origin.y + (config.domainHeight() - worldY) * viewport.scale};
    }
};

std::string windowTitle(const FLIPSimulator& simulator, bool paused, bool debugParticles) {
    std::ostringstream title;
    title << "PIC/FLIP | step " << simulator.currentStep
          << " | particles " << simulator.particles.size()
          << " | div " << std::scientific << std::setprecision(2)
          << simulator.divAfterProjection
          << (paused ? " | PAUSED" : "")
          << (debugParticles ? " | PARTICLES" : " | SMOOTH SURFACE");
    return title.str();
}

void resetSimulation(FLIPSimulator& simulator) {
    simulator.currentStep = 0;
    simulator.divBeforeProjection = 0.0f;
    simulator.divAfterProjection = 0.0f;
    simulator.initialize();
}

}

int main(int argc, char* argv[]) {
    using Clock = std::chrono::steady_clock;

    bool rebuildSurfaceEveryFrame = false;
    float profileSeconds = 0.0f;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string option = argv[argument];
        if (option == "--rebuild-surface-every-frame") {
            rebuildSurfaceEveryFrame = true;
        } else if (option == "--profile-seconds" && argument + 1 < argc) {
            profileSeconds = std::stof(argv[++argument]);
        }
    }
    SimConfig config;
    config.gridNX = 64;
    config.gridNY = 64;
    config.numSteps = 0;
    config.enableOutput = false;

    // Testiranje PIC/FLIP omjera i pritiska
    config.flipAlpha = 0.85f;
    config.pressureIterations = 250;

    FLIPSimulator simulator(config);
    simulator.initialize();

    sf::RenderWindow window(sf::VideoMode({1200, 850}), "PIC/FLIP");
    window.setFramerateLimit(60);

    if (!sf::Shader::isAvailable()) {
        std::cerr << "SFML shaders are not available on this system.\n";
        return 1;
    }
    sf::Shader waterShader;
    if (!waterShader.loadFromFile(WATER_SHADER_PATH, sf::Shader::Type::Fragment)) {
        std::cerr << "Could not load water shader: " << WATER_SHADER_PATH << "\n";
        return 1;
    }

    sf::RenderTexture backgroundLayer;
    sf::RenderTexture waterLayer;
    if (!backgroundLayer.resize(window.getSize()) || !waterLayer.resize(window.getSize())) {
        std::cerr << "Could not create SFML render textures.\n";
        return 1;
    }

    sf::Clock frameClock;
    float accumulator = 0.0f;
    bool paused = false;
    bool singleStep = false;
    bool debugParticles = false;
    sf::Texture densityTexture;
    sf::Texture velocityTexture;
    densityTexture.setSmooth(true);
    velocityTexture.setSmooth(true);
    SmoothSurfaceRenderer surfaceRenderer;
    sf::VertexArray surface;
    sf::VertexArray surfaceHighlights;
    sf::Image densityImage;
    sf::Image velocityImage;
    sf::Vector2u cachedWindowSize;
    int cachedSimulationStep = -1;

    double physicsMilliseconds = 0.0;
    double surfaceMilliseconds = 0.0;
    double textureUploadMilliseconds = 0.0;
    std::uint64_t measuredFrames = 0;
    std::uint64_t measuredSteps = 0;
    auto profileStart = Clock::now();
    const auto applicationStart = profileStart;

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window.close();
            }
            if (const auto* key = event->getIf<sf::Event::KeyPressed>()) {
                if (key->code == sf::Keyboard::Key::Space) {
                    paused = !paused;
                } else if (key->code == sf::Keyboard::Key::R) {
                    resetSimulation(simulator);
                    accumulator = 0.0f;
                    cachedSimulationStep = -1;
                } else if (key->code == sf::Keyboard::Key::N) {
                    singleStep = true;
                } else if (key->code == sf::Keyboard::Key::D) {
                    debugParticles = !debugParticles;
                } else if (key->code == sf::Keyboard::Key::Escape) {
                    window.close();
                }
            }
        }

        const float frameSeconds = std::min(frameClock.restart().asSeconds(), 0.1f);
        if (!paused) accumulator += frameSeconds;
        if (singleStep) {
            const auto physicsStart = Clock::now();
            simulator.step(config.dt);
            physicsMilliseconds += std::chrono::duration<double, std::milli>(
                Clock::now() - physicsStart).count();
            ++simulator.currentStep;
            ++measuredSteps;
            accumulator = 0.0f;
            singleStep = false;
        } else if (!paused) {
            constexpr int maxCatchUpSteps = 4;
            int catchUpSteps = 0;
            while (accumulator >= config.dt && catchUpSteps < maxCatchUpSteps) {
                const auto physicsStart = Clock::now();
                simulator.step(config.dt);
                physicsMilliseconds += std::chrono::duration<double, std::milli>(
                    Clock::now() - physicsStart).count();
                ++simulator.currentStep;
                ++measuredSteps;
                ++catchUpSteps;
                accumulator -= config.dt;
            }
            if (catchUpSteps == maxCatchUpSteps && accumulator > config.dt) {
                accumulator = config.dt;
            }
        }

        const Viewport viewport = makeViewport(config, window.getSize());

        sf::RectangleShape container({config.domainWidth * viewport.scale,
                                      config.domainHeight() * viewport.scale});
        container.setPosition(viewport.origin);
        container.setFillColor(sf::Color(10, 28, 43));
        container.setOutlineColor(sf::Color(188, 220, 230));
        container.setOutlineThickness(4.0f);

        backgroundLayer.clear(sf::Color(13, 20, 32));
        backgroundLayer.draw(container);
        backgroundLayer.display();

        if (rebuildSurfaceEveryFrame ||
            cachedSimulationStep != simulator.currentStep ||
            cachedWindowSize != window.getSize()) {
            const auto surfaceStart = Clock::now();
            surfaceRenderer.build(simulator, viewport, surface, surfaceHighlights,
                                  densityImage, velocityImage);
            surfaceMilliseconds += std::chrono::duration<double, std::milli>(
                Clock::now() - surfaceStart).count();

            const auto textureStart = Clock::now();
            if (!densityTexture.loadFromImage(densityImage) ||
                !velocityTexture.loadFromImage(velocityImage)) {
                std::cerr << "Could not update fluid field textures.\n";
                return 1;
            }
            textureUploadMilliseconds += std::chrono::duration<double, std::milli>(
                Clock::now() - textureStart).count();
            cachedSimulationStep = simulator.currentStep;
            cachedWindowSize = window.getSize();
        }
        densityTexture.setSmooth(true);
        velocityTexture.setSmooth(true);
        waterShader.setUniform("u_densityField", densityTexture);
        waterShader.setUniform("u_velocityField", velocityTexture);
        waterShader.setUniform("u_domainOrigin", viewport.origin);
        waterShader.setUniform("u_domainSize", sf::Vector2f(
            config.domainWidth * viewport.scale, config.domainHeight() * viewport.scale));
        waterShader.setUniform("u_fieldResolution", sf::Vector2f(
            static_cast<float>(densityImage.getSize().x),
            static_cast<float>(densityImage.getSize().y)));
        waterShader.setUniform("u_background", backgroundLayer.getTexture());
        waterShader.setUniform("u_resolution", sf::Vector2f(window.getSize()));
        waterShader.setUniform("u_lightPosition", sf::Vector2f(
            static_cast<float>(window.getSize().x) * 0.30f,
            static_cast<float>(window.getSize().y) * 0.18f));

        waterLayer.clear(sf::Color::Transparent);
        sf::RenderStates waterStates;
        waterStates.shader = &waterShader;
        waterLayer.draw(surface, waterStates);
        waterLayer.display();

        window.clear(sf::Color(13, 20, 32));
        window.draw(sf::Sprite(backgroundLayer.getTexture()));
        window.draw(sf::Sprite(waterLayer.getTexture()));

        if (debugParticles) {
            sf::CircleShape particle(2.0f);
            particle.setOrigin({2.0f, 2.0f});
            for (const Particle& value : simulator.particles) {
                particle.setPosition(toScreen(value, config, viewport));
                particle.setFillColor(particleColor(value));
                window.draw(particle);
            }
        }

        window.draw(surfaceHighlights);
        sf::RectangleShape border = container;
        border.setFillColor(sf::Color::Transparent);
        window.draw(border);
        window.setTitle(windowTitle(simulator, paused, debugParticles));
        window.display();

        ++measuredFrames;
        const double elapsedProfileSeconds = std::chrono::duration<double>(Clock::now() - profileStart).count();
        if (elapsedProfileSeconds >= 2.0) {
            const double averagePhysics = measuredSteps > 0
                ? physicsMilliseconds / static_cast<double>(measuredSteps) : 0.0;
            const double averageSurface = measuredSteps > 0
                ? surfaceMilliseconds / static_cast<double>(measuredSteps) : 0.0;
            const double averageUpload = measuredSteps > 0
                ? textureUploadMilliseconds / static_cast<double>(measuredSteps) : 0.0;
            std::cout << "[ViewerPerf] fps=" << measuredFrames / elapsedProfileSeconds
                      << " sim_steps_per_sec=" << measuredSteps / elapsedProfileSeconds
                      << " physics_ms/step=" << averagePhysics
                      << " surface_ms/step=" << averageSurface
                      << " texture_upload_ms/step=" << averageUpload
                      << (paused ? " paused" : "") << std::endl;
            measuredFrames = 0;
            measuredSteps = 0;
            physicsMilliseconds = 0.0;
            surfaceMilliseconds = 0.0;
            textureUploadMilliseconds = 0.0;
            profileStart = Clock::now();
        }
        if (profileSeconds > 0.0 &&
            std::chrono::duration<float>(Clock::now() - applicationStart).count() >= profileSeconds) {
            window.close();
        }
    }

    return 0;
}
