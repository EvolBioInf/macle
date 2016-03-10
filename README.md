# dnalc

Experimental local neighborhood complexity estimation for DNA,
based on match-length factorization (derived from Lempel-Ziv-factorization)
and periodicity counting (detection of repetitive DNA).

### Build

This program depends on [libdivsufsort](https://github.com/y-256/libdivsufsort)
and [GSL](https://www.gnu.org/software/gsl/).

The build itself is straightforward - obtain the source and run make:
```shell
git clone git@github.com:apirogov/dnalc.git
cd dnalc
make
```

This should leave you with a `dnalc` binary in the `build` directory.

(TODO: better description, usage examples)
