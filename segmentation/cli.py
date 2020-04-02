import warnings

from skimage.io import imsave, imread

from segmentation import segmentation


def main():
    import argparse
    parser = argparse.ArgumentParser('Segment yeast image')

    parser.add_argument('image')
    parser.add_argument('name')

    args = parser.parse_args()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # ignore classic "tags are not ordered by code"
        image = imread(args.image, plugin='tifffile')[1]

    image = segmentation.blur_frame(image)
    segmented, _ = segmentation.mixture_model(image, debug=True)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # ignore "tiff is a low contrast image"
        imsave('%s.tiff' % args.name, segmented)


if __name__ == '__main__':
    main()
