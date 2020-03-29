#cython: language_level=3

import numpy as np
cimport numpy as np
from skimage.filters import gaussian
from skimage.measure import label


cdef extern from "fastmm.cpp":
    cdef extern int cmm (const unsigned short* im, unsigned char* segmented, int m, int n, int nucNum)


def blur_frame(fr):
    return gaussian(fr, sigma=2, preserve_range=True).astype(np.uint16)


def blur_channel(channel):
    for frame in range(channel.shape[0]):
        channel[frame] = blur_frame(channel[frame])


cpdef mixture_model(const unsigned short[:, :] im):
    cdef int m, n, o, nuc_num
    cdef np.ndarray[np.uint8_t, ndim=2, mode="c"] segmented

    m = im.shape[0]
    n = im.shape[1]

    nuc_num = label(np.greater(im, 110), return_num=True)[1]
    segmented = np.zeros_like(im, dtype=np.uint8) # result

    cmm(&im[0, 0], &segmented[0, 0], m, n, nuc_num)

    return segmented, 0.
