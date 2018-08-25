# edt

Euclidean distance transform for 2D and 3D images together with a
MATLAB interface (mex). This is faster than `bwdist` for 3D images and
can also handle anisotropic voxels. Matlab's bwdist implements [3].

## Implementations:

Both MATLAB (`bwdist`) and Python3
(`SciPy.ndimage.distance_transform_edt`) uses [3].

```
Environment, Multi-Core, Anisotropic, N-dimensional, Method
Matlab,      Yes,        No,          Yes            [3]
Python3,     No,         Yes,         Yes            [3]
Eudist,      Yes,        Yes,         No             [2]
```

## Implementation
This is an implementation of the algorithm presented by Meijster et al. [1], with the following major differences:
 * 3D, not 2D.
 * Not only isotropic but also anisotropic voxels are handled.
 * Only the Euclidean distance transform, not the other alternatives
   that are discussed in the paper.
 * `lpthread` is used for parallelisation and hence the code should
   compile on both Linux and Mac.
 * Each thread has a private copy of three line buffers of size
   `max(M,N,P)` but very little more overhead.

## Notes
 * Mishchenko claims that their Matlab code for extending `eudist` in MATLAB from 2D to 3D is even faster than this. It is possible that it was when they wrote the paper but it clearly isn't today.

## Timings

For a `1042x1024x60` image with isotropic pixels:
```
bwdist: 2.178184
df_eudist: 0.897543
bwdistsc 8.824819
```

![2D timings](timings_2D.png)
![3D timings](timings_3D.png)

## TODO
 * Write tests for corner cases and for invalid input, then handle
   those cases. For example, crashes when more threads than pixels for a dimension.

## References:
 1. Mishchenko, 2012. Code on [matlab file exchange](https://se.mathworks.com/matlabcentral/fileexchange/15455-3d-euclidean-distance-transform-for-variable-data-aspect-ratio).
 2. Meijster et al., 2000.
 3. Maurer, Calvin, Rensheng Qi, and Vijay Raghavan, "A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions," IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 25, No. 2, February 2003, pp. 265-270.
