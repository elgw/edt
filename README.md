# edt

The Euclidean distance transform, D, of a binary image, B, sets each
voxel in D to the shortest distance to any non-zero voxel in B. If the
voxel is considered to have the same length in each dimension it is
called *isotropic*, otherwise *anisotropic*. The computational
burden increases quadratically with the number of pixels in the obvious, 
brute-force implementation. State of the art methods are linear in the
number of pixels and the number of dimensions, i.e., they are O(dN).

This repo contains an implementation of the Euclidean distance transform for 
 2D and 3D images together with wrappers for MATLAB and Python.

## Implementations

Here is a incomplete list of places where the Euclidean distance
transform can be found:

```
Environment, Implementation, Multi-Core, Anisotropic, Dimensions, Method
Matlab,      bwdist,         Yes,        No,          N           [3]
Matlab,      bwdistsc,       Yes,        Yes,         3           [1]
Python3 (1), ndimage,        No,         Yes,         N           [3]
C, (2)       eudist,         Yes,        Yes,         3           [2]
ImageJ                       ?           ?            2            ?
diplib                                                             ?
```

(1) `SciPy.ndimage.distance_transform_edt`
(2) with wrappers for Matlab and Python

## Implementation notes for eudist
This is an implementation of the algorithm presented by Meijster et al. [1], with the following major differences:
 * 3D, not 2D.
 * Not only isotropic but also anisotropic voxels are handled.
 * Only the Euclidean distance transform, not the other alternatives
   that are discussed in the paper.
 * `lpthread` is used for parallelisation and hence the code should
   compile on both Linux and Mac.
 * Each thread has a private copy of three line buffers of size
   `max(M,N,P)` but very little more overhead.
 * In the description of pass 2 and 3 on page ?, line ?, '<' is replaced by
   '<='.

## Usage:

### Standalone/debugging and timing

```
make eudist
./eudist -h
```

### Matlab
Run
```
df_eudist_ut
```
to compile and do some timings and testing.

### Python
```
make python
python3 scipy_timings.py
```
That will compile the python wrapper and do some timings.


## Timings

For a `1042x1024x60` image with isotropic pixels, on a Intel
i5-4690K@3700 MHz, it takes:

```
eudist:                  0.9 s
bwdist:                  2.2 s
bwdistsc                 8.8 s
distance_transform_edt: 17.4 s
```

![2D timings](timings_2D.png)
![3D timings](timings_3D.png)

## TODO
 * Write tests for corner cases and for invalid input, then handle
   those cases. For example, crashes when more threads than pixels for a dimension.
 * ND image implementation?
 * Finnish Python wrapper.

## References:
 1. Mishchenko, 2012. Code on [matlab file exchange](https://se.mathworks.com/matlabcentral/fileexchange/15455-3d-euclidean-distance-transform-for-variable-data-aspect-ratio).
 2. Meijster et al., 2000.
 3. Maurer, Calvin, Rensheng Qi, and Vijay Raghavan, "A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions," IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 25, No. 2, February 2003, pp. 265-270.
