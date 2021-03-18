![Python package](https://github.com/BooneAndrewsLab/segmentation/workflows/Python%20package/badge.svg)
![Python package](https://github.com/BooneAndrewsLab/segmentation/workflows/Python%20package%20build%20and%20publish/badge.svg)
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
### Programmatically:
```python
from segmentation import segmentation
from segmentation import watershed
from skimage.io import imread

image = imread('./image.tiff', plugin='tifffile')[1]  # Read channel 1 of a tiff/flex 
im = segmentation.blur_frame(image) # gaussian blur
segmented, _ = segmentation.mixture_model(im, debug=True) # second return argument is currently unused
labels = watershed(im, segmented)
```

### Command line
```sh
$ segment -h # for usage information

$ segment -o segmented.data image.tiff
```

Output is a memmaped labels array. You can read it like this:
```python
from numpy import memmap
labels = memmap('segmented.data', dtype='int32', shape=(1010, 1346))  # shape is same as input image
```
