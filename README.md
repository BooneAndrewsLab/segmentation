![Python package](https://github.com/BooneAndrewsLab/segmentation/workflows/Python%20package/badge.svg)

# Fast Mixture Model Segmentation

Fast mixture model segmentation used in Boone and Andrews labs

# Python requirements

  - Numpy
  - Cython
  - scikit-image

## Installation

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

## Usage example

```python
import segmentation as seg
from skimage.io import imread

image = imread('./001001000.tiff', plugin='tifffile')[1]  # Read channel 1 of a tiff/flex 
im = seg.blur_frame(image) # gaussian blur
segmented, _ = seg.mixture_model(im, debug=True) # second return argument is currently unused
labels = seg.watershed(im, segmented)
```
