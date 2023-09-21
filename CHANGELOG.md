## version 0.1.2
- Switched to openmp for parallelization.
- Removed the number of threads argument from edt.

## version 0.1.1
- Fixed an issue causing the wrong results when anisotropic pixel
  sizes were used.
- Sets the number of threads automatically if not specified as `sysconf(_SC_NPROCESSORS_ONLN)/2`
