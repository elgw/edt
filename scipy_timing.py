#!/bin/python3

import numpy as np
import time
from scipy import ndimage
import random
import eudist
import matplotlib.pyplot as plt

def test_size2D(M, N):
    print("Testing {:d}x{:d}".format(M, N));
    B = np.zeros([M, N]);
    for kk in range(1,100):
        pos = random.randint(0, M*N-1);
        B.flat[pos] = 1;
    tic = time.time()
    D = ndimage.distance_transform_edt(B==0)
    toc = time.time()    

    tic2 = time.time()
    S = np.zeros([3,1]);
    S[0] = 1; S[1] = 1; S[2] = 1;
    D2 = eudist.eudist(np.asfortranarray(B), S); 
    print(D2.shape);
    toc2 = time.time()
    print("Scipy  took {:f} s".format(toc-tic));
    print("Eudist took {:f} s".format(toc2-tic2));
    max_diff = np.amax(np.abs(D-D2));
    print("Max difference: {:f}".format(max_diff));

    plt.subplot(221);
    plt.imshow(B);    
    plt.subplot(222);
    plt.imshow(D);
    plt.subplot(223);
    plt.imshow(D2);
    plt.show();

def test_size3D(M, N, P):
    print("Testing {:d}x{:d}x{:d}".format(M, N, P));
    B = np.zeros([M, N, P]);
    for kk in range(1,20):
        pos = random.randint(0, M*N*P-1);
        B.flat[pos] = 1;

    tic = time.time()
    D = ndimage.distance_transform_edt(B==0, sampling=[1,1,1],
                                       return_distances=True);
    toc = time.time()
    S = np.zeros([3,1]);
    S[0] = 1; S[1] = 1; S[2] = 1;
    tic2 = time.time()
    print(B.shape);
    D2 = eudist.eudist(np.asfortranarray(B), S); 
    print(D2.shape);
    toc2 = time.time()
    print("Scipy took  {:f} s".format(toc-tic));
    print("Eudist took {:f} s".format(toc2-tic2));
    max_diff = np.amax(np.abs(D-D2));
    print("Max difference: {:f}".format(max_diff));

    sl = 0;
    Bs = B[:,:,sl];
    Ds = D[:,:,sl];
    D2s = D2[:,:,sl];

    plt.subplot(131);
    plt.imshow(Bs);
    plt.subplot(132);
    plt.imshow(Ds);
    plt.subplot(133);
    plt.imshow(D2s);
    plt.show();


test_size2D(1024, 1024);
test_size3D(1024, 1024, 60);


