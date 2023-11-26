#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>

#include <random>
#include <vector>

namespace nbody {

using real = double;

struct particle {
  // Position
  real rx;
  real ry;
  real rz;

  // Velocity
  real vx;
  real vy;
  real vz;
};

struct ensemble {
  real mass;                       // Uniform mass for all particles
  std::vector<particle> particles; // Particle data

  // STL containers use unsigned integers for container sizes. That is
  // recognized as a design mistake by most C++ experts nowadays. Unsigned
  // integers convey the correct semantics but do not enforce it. Mixing signed
  // and unsigned arithmetic is dangerous, unsigned integers "wrap around" when
  // they overflow. Use signed integers when you can, use unsigned only for
  // low level bitwise operations.
  //
  int num_particles() const noexcept {
    return static_cast<int>(particles.size());
  }
};

// Setup an initial random state for the simulation.
//
// We set up the initial state as a uniform distribution
// over the unit box [0,1)^3.
//
std::vector<particle> generate_particles(int const n) {
  constexpr real vel_fact = 0.005;

  // Seed for the pseudo random number generator. To have a different seed every
  // time the application launches, use `std::random_device`.
  //
  constexpr int rseed = 3;

  // Generating random numbers with C++11's random requires an engine and a
  // distribution. This is an engine based on the Mersenne Twister 19937.
  //
  auto eng = std::mt19937{rseed};

  // Create a uniform distribution [0, 1).
  //
  auto dist = std::uniform_real_distribution<double>{0., 1.};

  std::vector<particle> particles;
  particles.reserve(n);

  for (int idx = 0; idx < n; ++idx) {
    // Create a particle with random position and velocity.
    //
    particles.push_back(particle{dist(eng), dist(eng), dist(eng),
                                 vel_fact * dist(eng), vel_fact * dist(eng),
                                 vel_fact * dist(eng)});
  }

  return particles;
}

// Write simulation to disk.
//
void write_state(ensemble const &sim, const int t) {
  char fn[128];
  sprintf(fn, "output/step-%06d.txt", t);

  auto fd = fopen(fn, "w");
  if (!fd) {
    fprintf(stderr, "[ERROR] Could not open file %s\n", fn);
    exit(42);
  }
  fprintf(fd, "# n rx ry rz vx vy vz\n");

  int n = sim.num_particles();
  for (int idx = 0; idx < n; ++idx) {
    auto const &p = sim.particles[idx];
    fprintf(fd, "%d %lf %lf %lf\t %lf %lf %lf\n", idx, p.rx, p.ry, p.rz, p.vx,
            p.vy, p.vz);
  }
  fclose(fd);
}

// Update velocities (Euler integrations)
//    v_i = v_i + dt * F_i / m_i
//
// Every particle experiences the force
//             ___
//            \        (r_i - r_j) m_i m_j
//     F_i =   >     ----------------------- = m_i a_i
//            /___        |r_i - r_j|^3
//           j =/= i
//
void update_velocities(ensemble &sim, const real dt) {

  // A minimal allowed distance to avoid particles getting infinitely close.
  // This has two consequences
  //
  // 1. Forces cannot become arbitrarily large
  // 2. An error is introduced
  //
  constexpr auto eps = 1e-9;

  // Loop over target particle
  //
  int n = sim.num_particles();
  for (auto idx = 0; idx < n; ++idx) {
    auto &pi = sim.particles[idx];
    auto Fx = 0.;
    auto Fy = 0.;
    auto Fz = 0.;

    // Loop over interaction partners
    //
    for (auto jdx = 0; jdx < n; ++jdx) {
      auto const &pj = sim.particles[jdx];
      // No self-force
      //
      if (idx != jdx) {
        // Distance to partner
        //
        auto dx = pi.rx - pj.rx;
        auto dy = pi.ry - pj.ry;
        auto dz = pi.rz - pj.rz;
        auto inv_dist2 = 1.0 / (dx * dx + dy * dy + dz * dz + eps);

        // Sum up force
        //
        Fx += std::copysign(inv_dist2, dx);
        Fy += std::copysign(inv_dist2, dy);
        Fz += std::copysign(inv_dist2, dz);
      }
    }

    // Update velocity
    //
    pi.vx += sim.mass * dt * Fx;
    pi.vy += sim.mass * dt * Fy;
    pi.vz += sim.mass * dt * Fz;
  }
}

// Update positions. (Euler integration)
//
//    r_i = r_i + dt * v_i
//
void update_positions(ensemble &sim, real const dt) {
  for (auto &p : sim.particles) {
    p.rx += p.vx * dt;
    p.ry += p.vy * dt;
    p.rz += p.vz * dt;
  }
}

// Run Simulation
//
// - Do n_run steps of force + move
// - Output every n_out
// - Finally show timings
//
void simulate(ensemble &sim, const real dt, const int n_step, const int n_out) {
  using clock = std::chrono::high_resolution_clock;
  using seconds = std::chrono::duration<double>;

  auto t_force = seconds::zero();
  auto t_move = seconds::zero();
  for (int t = 0; t <= n_step; ++t) {
    auto const t0 = clock::now();

    update_velocities(sim, dt);

    auto const t1 = clock::now();

    update_positions(sim, dt);

    auto const t2 = clock::now();

    t_force += t1 - t0;
    t_move += t2 - t1;

    if (t % n_out == 0) {
      auto avg_force = t_force / (t + 1);
      auto avg_move = t_move / (t + 1);
      printf("[INFO] Step %d:\n"
             " * avg force: %.4lf s\n"
             " * avg move:  %.4lf s\n",
             t, avg_force.count(), avg_move.count());
      write_state(sim, t);
    }
  }
  printf("[INFO] Done\n"
         " * tot force: %.4lf s\n"
         " * tot move:  %.4lf s\n",
         t_force.count(), t_move.count());
}

} // end namespace nbody
