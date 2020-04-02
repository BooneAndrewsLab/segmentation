#cython: language_level=3

import numpy as np
cimport numpy as npc
from skimage.filters import gaussian
from skimage.measure import label
from libcpp cimport bool  # for bool support

npc.import_array()

cdef extern from "fastmm.cpp":
    cdef struct CmmConfig:
        double learning_rate
        double var
        int m_nuc_size
        int m_cell_size
        char max_iterations
        bool debug
    ctypedef void (*cbfunc)(int step, double *x, int m, int n, int o)
    cdef extern int cmm (const unsigned short* im, unsigned char* segmented, int m, int n, int nucNum, CmmConfig config, cbfunc callback)


def blur_frame(fr):
    return gaussian(fr, sigma=2, preserve_range=True).astype(np.uint16)


def blur_channel(channel):
    for frame in range(channel.shape[0]):
        channel[frame] = blur_frame(channel[frame])

cdef void segment_step_callback(int step, double *x, int m, int n, int o):
    cdef int nd = 1
    cdef npc.npy_intp shape[1]

    shape[0] = m * n * o

    arr = npc.PyArray_SimpleNewFromData(nd, shape, np.NPY_DOUBLE, x)
    arr = arr.reshape(m, n, o)

    np.save('step%02d.npy' % step, arr)

    del arr

cpdef mixture_model(const unsigned short[:, :] im):
    cdef int m, n, nuc_num
    cdef npc.ndarray[npc.uint8_t, ndim=2, mode="c"] segmented
    cdef CmmConfig conf
    # conf.max_iterations = 10

    nuc_num = label(np.greater(im, 110), return_num=True)[1]
    segmented = np.zeros_like(im, dtype=np.uint8) # result

    cmm(&im[0, 0], &segmented[0, 0], im.shape[0], im.shape[1], nuc_num, conf, NULL)

    return segmented, 0.
