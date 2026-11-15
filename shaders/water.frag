// GLSL compatibility profile, for SFML 3 sf::Shader.
// No time-dependent procedural noise: all changes come from the simulation.
uniform sampler2D u_background;
uniform sampler2D u_densityField;
uniform sampler2D u_velocityField;
uniform vec2 u_resolution;
uniform vec2 u_domainOrigin; // top-left in SFML pixels
uniform vec2 u_domainSize;   // width, height in pixels
uniform vec2 u_fieldResolution;
uniform vec2 u_lightPosition; // top-left pixel coordinates

void main() {
    vec2 pixelTop = vec2(gl_FragCoord.x, u_resolution.y - gl_FragCoord.y);
    vec2 fieldUV = (pixelTop - u_domainOrigin) / u_domainSize;
    // SFML texture coordinates are top-left based.
    vec2 uv = clamp(fieldUV, vec2(0.0), vec2(1.0));
    vec2 texel = 1.0 / u_fieldResolution;
    float density = texture2D(u_densityField, uv).r;
    float dx = texture2D(u_densityField, clamp(uv + vec2(texel.x, 0.0), 0.0, 1.0)).r
             - texture2D(u_densityField, clamp(uv - vec2(texel.x, 0.0), 0.0, 1.0)).r;
    float dy = texture2D(u_densityField, clamp(uv + vec2(0.0, texel.y), 0.0, 1.0)).r
             - texture2D(u_densityField, clamp(uv - vec2(0.0, texel.y), 0.0, 1.0)).r;
    vec2 outward = -vec2(dx, dy);
    float edgeStrength = clamp(length(outward) * 2.0, 0.0, 1.0);
    vec2 edgeNormal = outward / max(length(outward), 0.0001);

    // Encoded velocities from MAC u/v faces, range approx [-8,8] m/s.
    vec2 localVelocity = (texture2D(u_velocityField, uv).rg - 0.5) * 16.0;
    // Note: SFML image is top-down; physical y points up.
    vec2 velocityOnScreen = vec2(localVelocity.x, -localVelocity.y);

    // The surface contour controls the edge highlight. This is a 2D optical
    // approximation, not a physically correct 3D Fresnel/refraction model.
    vec2 towardsLight = normalize((u_lightPosition - pixelTop) / u_resolution + vec2(0.0001));
    float rim = pow(max(dot(edgeNormal, towardsLight), 0.0), 2.5) * edgeStrength;

    // Refraction is strongest at actual density gradients and changes with
    // local fluid velocity, not wall-clock time or scrolling procedural noise.
    vec2 distortionPixels = edgeNormal * edgeStrength * 4.0
                          + clamp(velocityOnScreen * 0.30, vec2(-2.5), vec2(2.5));
    vec2 refracted = clamp((pixelTop + distortionPixels) / u_resolution,
                           vec2(0.0), vec2(1.0));
    vec3 background = texture2D(u_background, refracted).rgb;

    // Density is a rendering estimate, not actual optical path length.
    float opticalThickness = smoothstep(0.08, 0.85, density);
    vec3 shallow = vec3(0.12, 0.69, 0.79);
    vec3 deep = vec3(0.018, 0.22, 0.44);
    vec3 water = mix(shallow, deep, opticalThickness * 0.78);
    float transmission = mix(0.45, 0.10, opticalThickness);
    water = mix(water, background, transmission);
    water += vec3(0.25, 0.41, 0.47) * rim;
    // No free-running texture motion: paused simulation -> identical output.
    gl_FragColor = vec4(clamp(water, 0.0, 1.0), 0.94);
}
