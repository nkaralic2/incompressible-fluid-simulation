#ifndef PARTICLE_HPP
#define PARTICLE_HPP

/**
 * @file Particle.hpp
 * @brief Particle data structure for PIC/FLIP fluid simulation.
 *
 * Each particle carries position and velocity in 2D.
 * The coordinate convention is:
 *   - Positive X points right.
 *   - Positive Y points upward.
 *   - Gravity acts in the negative Y direction.
 */

struct Particle {
    float x, y;   // Position in world coordinates
    float vx, vy; // Velocity components
};

#endif // PARTICLE_HPP
