#distutils: language=c++
#cython: language_level=3

import numpy as np
cimport numpy as npc
from skimage.filters import gaussian
from skimage.measure import label
from libcpp cimport bool  # for bool support

npc.import_array()  # otherwise we get segfaults while wrapping double* in ndarray

cdef extern from "fastmm.h":
    ctypedef void (*cbfunc)(int step, double *x, int m, int n, int o, void *user_data)
    cdef struct CmmConfig:
        double learning_rate
        double var
        int m_nuc_size
        int m_cell_size
        char max_iterations
        bool debug
    cdef extern int cmm (const unsigned short* im, unsigned char* segmented, int dim1, int dim2, int nuc_num,
                         CmmConfig config, cbfunc callback, void *user_data)


def blur_frame(fr):
    return gaussian(fr, sigma=2, preserve_range=True).astype(np.uint16)


def blur_channel(channel):
    for frame in range(channel.shape[0]):
        channel[frame] = blur_frame(channel[frame])

cdef void segment_step_callback(int step, double *x, int m, int n, int o, void *user_data):
    cdef int nd = 1
    cdef npc.npy_intp shape[1]
    pycb = <object>user_data

    shape[0] = m * n * o
    arr = npc.PyArray_SimpleNewFromData(nd, shape, npc.NPY_DOUBLE, x)
    arr = arr.reshape(m, n, o)

    pycb(arr)

    del arr

def mixture_model(const unsigned short[:, :] im, learning_rate=None, var=None, nuc_size=None, cell_size=None,
                  max_iter=None, debug=None, intermediate_seg_callback=None):
    """ Segment im using mixture model implemented in C++

    :param im: image to segment
    :param learning_rate: cmm learning rate, default: 1e-6
    :param var: cmm variance, default: 1e8
    :param nuc_size: estimated average size of nucleus, default: 210
    :param cell_size: estimated average cell size, default: 1900
    :param max_iter: stop segmentation after this many iterations, default: 100
    :param debug: be verbose, default: False
    :param intermediate_seg_callback: receive segmentation progress at every iteration
    :type im: np.ndarray
    :type learning_rate: int
    :type var: int
    :type nuc_size: int
    :type cell_size: int
    :type max_iter: int
    :type debug: bool
    :type intermediate_seg_callback: function
    :return: segmented image
    """
    cdef int nuc_num
    # Preallocate output array
    cdef npc.ndarray[npc.uint8_t, ndim=2, mode="c"] segmented
    cdef CmmConfig conf

    # Set config variables
    if learning_rate:
        conf.learning_rate = learning_rate

    if var:
        conf.var = var

    if nuc_size:
        conf.m_nuc_size = nuc_size

    if cell_size:
        conf.m_cell_size = cell_size

    if max_iter:
        conf.max_iterations = max_iter

    if debug:
        conf.debug = debug

    # Estimate how many nuclei we have
    nuc_num = label(np.greater(im, 110), return_num=True)[1]
    segmented = np.empty_like(im, dtype=np.uint8) # don't initialize value

    if intermediate_seg_callback:
        cmm(&im[0, 0], &segmented[0, 0], im.shape[0], im.shape[1], nuc_num, conf,
            <cbfunc>segment_step_callback, <void*>intermediate_seg_callback)
    else:
        cmm(&im[0, 0], &segmented[0, 0], im.shape[0], im.shape[1], nuc_num, conf, NULL, NULL)

    return segmented, 0.
