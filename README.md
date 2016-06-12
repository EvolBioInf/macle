# dnalc

Experimental local neighborhood complexity estimation for DNA,
based on match-length factorization (derived from Lempel-Ziv-factorization)
and periodicity counting (detection of repetitive DNA).

### Build

This program depends on [libdivsufsort](https://github.com/y-256/libdivsufsort)
and [GSL](https://www.gnu.org/software/gsl/). If you want to work with very long
sequences, you also need [SDSL](https://github.com/simongog/sdsl-lite).

The build itself is straightforward - obtain the source and run make:
```shell
git clone git@github.com:apirogov/dnalc.git
cd dnalc
make
```

If you don't plan to work with very long sequences (> 2^32), you can disable
SDSL in the Makefile or by invoking make like this: `make USE_SDSL=0`. Then a
faster 32-bit version will be built.

If you can't install a dependency using your package manager, you can
run `make [libdivsufsort|gsl|sdsl]` to automatically download and build it.

After the build is complete you can find the `dnalc` binary in the `build` directory.

(TODO: better description, usage examples)
