#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

const double G = 6.67430e-11; // Gravitational constant
const double c = 299792458;   // Speed of light
const double c2 = c * c;

struct Vector3 {
    double x, y, z;

    Vector3() : x(0), y(0), z(0) {} 
    Vector3(double x, double y, double z) : x(x), y(y), z(z) {} 

    double length() const {
        return sqrt(x*x + y*y + z*z);
    }

    Vector3 normalize() const {
        double len = length();
        if (len > 0) {
            return Vector3(x/len, y/len, z/len);
        }
        return *this;
    }
};

Vector3 operator+(const Vector3& a, const Vector3& b) {
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector3 operator-(const Vector3& a, const Vector3& b) {
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3 operator*(const Vector3& a, double scalar) {
    return Vector3(a.x * scalar, a.y * scalar, a.z * scalar);
}

Vector3 operator*(double scalar, const Vector3& a) {
    return a * scalar;
}


class BlackHole {
    private:
        double mass; 
        Vector3 position; 

    public:
        BlackHole(double mass, const Vector3& pos): mass(mass), position(pos){}

            double getSchwarzschildRadius() const {
                return (2 * G * mass) / c2;
            }

            Vector3 getGravitationalAcceleration(const Vector3& point) const {
                Vector3 r_vec = position - point;
                double r = r_vec.length();

                if (r < 1e-10) {
                    return Vector3 (0,0,0);
                }

                double acceleration_magnitude = (G * mass) / (r * r);
                return r_vec.normalize() * acceleration_magnitude;
            }
       

        double getTimeDilationFactor(const Vector3& point) const {
            Vector3 r_vec = position - point;
            double r = r_vec.length();
            double rs = getSchwarzschildRadius();
            
            if (r <= rs) {
                return 0.0; // Inside event horizon, time stops
            }

            return sqrt(1 - rs / r);
        }

        double getLightDeflection(double impact_parameter) const {
            double rs = getSchwarzschildRadius();
            return (4 * G * mass) / (c2 * impact_parameter);
        }

        Vector3 getPosition() const {
            return position;
        }

        double getMass() const {
            return mass;
        }
};


class Particle {
    private: 
        Vector3 position;
        Vector3 velocity;
        double mass;
    public:
        Particle(const Vector3& pos, const Vector3& vel, double m)
            : position(pos), velocity(vel), mass(m) {}

    void update(double dt, const Vector3& acceleration) {
        velocity = velocity + acceleration * dt;
        position = position + velocity * dt;
    }

    Vector3 getPosition() const { return position; }
    Vector3 getVelocity() const { return velocity; }
    double getMass() const { return mass; }
};


class BlackHoleSimulation {
    private:
        BlackHole blackhole;
        std::vector<Particle> particles;
        double time_step;

    public:
        BlackHoleSimulation(const Blackhole& bh, double dt)
        : blackhole(bh), time_step(dt) {}

        void addParticle(const Particle& p) {
            particles.push_back(p);
        }

        void simulateStep() {
            for (auto& particle : particles)
        }
}