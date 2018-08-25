#!/bin/python3

import numpy as np
import time
from scipy import ndimage
import random

def test_size2D(M, N):
    B = np.zeros([M, N]);
    for kk in range(1,100):
        pos = random.randint(0, M*N-1);
        B.flat[pos] = 1;
    tic = time.time()
    D = ndimage.distance_transform_edt(B)
    toc = time.time()
    print(toc-tic);

def test_size3D(M, N, P):
    B = np.zeros([M, N, P]);
    for kk in range(1,100):
        pos = random.randint(0, M*N*P-1);
        B.flat[pos] = 1;
    tic = time.time()
    D = ndimage.distance_transform_edt(B)
    toc = time.time()
    print(toc-tic);



test_size2D(1024, 1024);
test_size3D(1024, 1024, 60);

