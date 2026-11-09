
#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
    double x;
    double y;

    double u;
    double v;

    Particle(double x, double y)
        : x(x), y(y), u(0.0), v(0.0) {
    }
};

#endif