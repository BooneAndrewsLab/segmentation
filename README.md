![Python package](https://github.com/BooneAndrewsLab/segmentation/workflows/Python%20package/badge.svg)

# Fast Mixture Model Segmentation

Fast mixture model segmentation used in Boone and Andrews labs

# Python requirements

  - Numpy
  - Cython
  - scikit-image

# C requirements
  - Boost

## Installation

Install C dependencies (ubuntu):

```sh
$ sudo apt install libboost-dev
```

Create a virtual environment (optional)
```sh
$ virtualenv -ppython3 segmentation-env
$ source segmentation-env/bin/activate
```

Install python requirements (needed to build the package)

```sh
$ pip install numpy cython
```

Install our library (pulls in all other dependencies)

```sh
$ pip install segmentation
```

## Usage

Read the image. 

read_channel takes an optional argument channel= that can be used
to specifiy which channel to read. "red" (default) is preset to begin 
with the second frame, "green" is first frame. The argument can
also be an integer, that can be used to manually specifiy at which
frame should we start reading.

```python
import segmentation as seg

tiff = seg.read_channel('./001001000.tiff', channel="red") # channel red == read every other frame, start from frame 2
im = seg.blur_frame(tiff[0]) # gaussian blur
labels, _ = seg.mixture_model(im) # second return argument is currently unused
```
