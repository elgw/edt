This repo contains an implementation of the Euclidean Distance
Transform (EDT) based on the work by Meijster et. al [^1] with
wrappers for MATLAB and Python.

The Euclidean distance transform, D, of a binary image, B, sets each
voxel in D to the shortest distance to any non-zero voxel in B. If the
voxel is considered to have the same length in each dimension it is
called *isotropic*, otherwise *anisotropic*. The computational
burden increases quadratically with the number of pixels in the obvious,
brute-force implementation. State of the art methods are linear in the
number of pixels and the number of dimensions, i.e., they are O(dN).

This repo contains an implementation of the Euclidean distance
transform for 2D and 3D images together with wrappers for MATLAB and
Python.


## Implementation
Differences to [:^1]

* Also 3D, however not nD :(
* Also handles anisotropic voxels.
* Only the Euclidean distance transform, not the other alternatives
  that are discussed in the paper.
* `lpthread` is used for parallelisation and hence the code should
  compile on both Linux and Mac.
* In the description of pass 2 and 3 on page ?, line ?, '<' is replaced by
'<='.


## Timings

For a `1042x1024x60` image with isotropic pixels, on a Intel
i5-4690K@3700 MHz, it took:

```
eudist:                  1.7 s
bwdist:                  2.2 s (MATLAB)
distance_transform_edt: 17.4 s (SCIPY)
```
timings were done in 2018... ish.

![2D timings](doc/timings_2D.png)
![3D timings](doc/timings_3D.png)

## TODO
* Write tests for corner cases and for invalid input.
* ND image implementation?
* Especially the Python wrapper needs some attention.
* The scheduling of the threads could for sure be optimized to handle
more exotic cases better.
* Switch to OMP to save some lines.
* Redo the documentation.
* Make a proper wheel or similar for Python

## References:
[^1]: Meijster, A., Roerdink, J.B.T.M., Hesselink, W.H. (2002). A General Algorithm for Computing Distance Transforms in Linear Time. In: Goutsias, J., Vincent, L., Bloomberg, D.S. (eds) Mathematical Morphology and its Applications to Image and Signal Processing. Computational Imaging and Vision, vol 18. Springer, Boston, MA. [doi:10.1007/0-306-47025-X_36](https://doi.org/10.1007/0-306-47025-X_36
