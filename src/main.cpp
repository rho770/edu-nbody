#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <nbody.hpp>
#include <vector>

int main(int argc, char *argv[]) {
  using namespace nbody;

  if (argc != 4) {
    fprintf(stderr, " [ERROR] Wrong number of parameters.\n\n"
                    " Usage: ./exe <n_part> <n_steps> <n_output> \n");
    return 42;
  }

  if (!opendir("output")) {
    fprintf(stderr,
            " [ERROR] The `output` directory does not exist in the current "
            "working directory. \n\n"
            " Please create it and try again. \n");
    exit(42);
  }

  int n = atoi(argv[1]);
  int steps = atoi(argv[2]);
  int out_freq = atoi(argv[3]);

  constexpr real mass = 0.0025;
  constexpr real dt = 0.0025;
  ensemble sys{mass, generate_particles(n)};
  simulate(sys, dt, steps, out_freq);
}
