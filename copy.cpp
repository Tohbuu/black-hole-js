#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

// Constants
const double G = 6.67430e-11;  // Gravitational constant
const double c = 299792458.0;  // Speed of light
const double c2 = c * c;       // c squared

// Structure for 3D vector
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

// Operator overloads for Vector3
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

// Black hole class
class BlackHole {
private:
    double mass;        // Mass in kg
    Vector3 position;   // Position in space
    
public:
    BlackHole(double mass, const Vector3& pos) : mass(mass), position(pos) {}
    
    // Calculate Schwarzschild radius
    double getSchwarzschildRadius() const {
        return (2 * G * mass) / c2;
    }
    
    // Calculate gravitational acceleration at a point
    Vector3 getGravitationalAcceleration(const Vector3& point) const {
        Vector3 r_vec = position - point;
        double r = r_vec.length();
        
        // Avoid division by zero
        if (r < 1e-10) {
            return Vector3(0, 0, 0);
        }
        
        double acceleration_magnitude = (G * mass) / (r * r);
        return r_vec.normalize() * acceleration_magnitude;
    }
    
    // Calculate time dilation factor at a point (relative to infinity)
    double getTimeDilationFactor(const Vector3& point) const {
        Vector3 r_vec = position - point;
        double r = r_vec.length();
        double rs = getSchwarzschildRadius();
        
        if (r <= rs) {
            return 0;  // Inside event horizon
        }
        
        return sqrt(1 - rs / r);
    }
    
    // Calculate light deflection (angle of deflection in radians)
    double getLightDeflection(double impact_parameter) const {
        double rs = getSchwarzschildRadius();
        return (4 * G * mass) / (c2 * impact_parameter);
    }
    
    // Getter for position
    Vector3 getPosition() const {
        return position;
    }
    
    // Getter for mass
    double getMass() const {
        return mass;
    }
};

// Particle class for simulation
class Particle {
private:
    Vector3 position;
    Vector3 velocity;
    double mass;
    
public:
    Particle(const Vector3& pos, const Vector3& vel, double m) 
        : position(pos), velocity(vel), mass(m) {}
    
    // Update particle position based on velocity and acceleration
    void update(double dt, const Vector3& acceleration) {
        velocity = velocity + acceleration * dt;
        position = position + velocity * dt;
    }
    
    // Getters
    Vector3 getPosition() const { return position; }
    Vector3 getVelocity() const { return velocity; }
    double getMass() const { return mass; }
};

// Simulation class
class BlackHoleSimulation {
private:
    BlackHole blackhole;
    std::vector<Particle> particles;
    double time_step;
    
public:
    BlackHoleSimulation(const BlackHole& bh, double dt) 
        : blackhole(bh), time_step(dt) {}
    
    void addParticle(const Particle& p) {
        particles.push_back(p);
    }
    
    void simulateStep() {
        for (auto& particle : particles) {
            Vector3 accel = blackhole.getGravitationalAcceleration(particle.getPosition());
            particle.update(time_step, accel);
        }
    }
    
    void runSimulation(int steps, const std::string& output_filename) {
        std::ofstream outfile(output_filename);
        if (!outfile.is_open()) {
            std::cerr << "Error opening output file!" << std::endl;
            return;
        }
        
        // Write header
        outfile << "time,x,y,z,vx,vy,vz\n";
        
        for (int i = 0; i < steps; ++i) {
            double current_time = i * time_step;
            
            // Write particle data
            for (const auto& particle : particles) {
                Vector3 pos = particle.getPosition();
                Vector3 vel = particle.getVelocity();
                
                outfile << current_time << ","
                        << pos.x << "," << pos.y << "," << pos.z << ","
                        << vel.x << "," << vel.y << "," << vel.z << "\n";
            }
            
            simulateStep();
        }
        
        outfile.close();
    }
};

// Function to visualize spacetime curvature around black hole
void visualizeSpacetimeCurvature(const BlackHole& bh, double range, int resolution) {
    std::ofstream outfile("spacetime_curvature.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }
    
    outfile << "x,y,curvature,time_dilation\n";
    
    double step = range * 2 / resolution;
    Vector3 bh_pos = bh.getPosition();
    
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double x = -range + i * step;
            double y = -range + j * step;
            Vector3 point(x, y, 0);
            
            // Calculate curvature (simplified as gravitational potential)
            Vector3 r_vec = bh_pos - point;
            double r = r_vec.length();
            double curvature = (r > 0) ? -G * bh.getMass() / r : 0;
            
            // Calculate time dilation
            double time_dilation = bh.getTimeDilationFactor(point);
            
            outfile << x << "," << y << "," << curvature << "," << time_dilation << "\n";
        }
    }
    
    outfile.close();
}

int main() {
    // Create a black hole with 10 solar masses
    double solar_mass = 1.989e30;
    BlackHole bh(10 * solar_mass, Vector3(0, 0, 0));
    
    std::cout << "Black hole created with mass: " << bh.getMass() << " kg" << std::endl;
    std::cout << "Schwarzschild radius: " << bh.getSchwarzschildRadius() << " m" << std::endl;
    
    // Create simulation
    BlackHoleSimulation sim(bh, 0.1);
    
    // Add some particles
    sim.addParticle(Particle(Vector3(1e9, 1e9, 0), Vector3(-1e4, 0, 0), 1e10));
    sim.addParticle(Particle(Vector3(-1e9, 1e9, 0), Vector3(0, -1e4, 0), 1e10));
    sim.addParticle(Particle(Vector3(0, 2e9, 0), Vector3(1.5e4, 0, 0), 1e10));
    
    // Run simulation
    sim.runSimulation(1000, "blackhole_simulation.csv");
    
    // Visualize spacetime curvature
    visualizeSpacetimeCurvature(bh, 3e9, 100);
    
    std::cout << "Simulation completed. Results saved to files." << std::endl;
    
    return 0;
}