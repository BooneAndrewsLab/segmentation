import logging
import warnings

from numpy import memmap
from skimage.io import imread

from . import segmentation
from .watershed import watershed

log = logging.getLogger(__file__)
log.setLevel(logging.ERROR)


def main():
    import argparse
    parser = argparse.ArgumentParser('Segment yeast image')

    parser.add_argument('image')
    parser.add_argument('-o', '--output', help='Save raw segmentation output')
    parser.add_argument('-c', '--nuclei-channel', help='Channel with nuclei for segmentation, zero-based, default: 1',
                        default=1, type=int)
    parser.add_argument('-l', '--learning-rate', help='Segmentation learning rate, default: 0.000001 (1e-6)',
                        default=1e-6, type=float)
    parser.add_argument('-r', '--variance', help='Variance used in segmentation, default: 100000000 (1e8)',
                        default=1e8, type=int)
    parser.add_argument('-n', '--nucleus-size', help='Estimated average size of nucleus, default: 210',
                        default=210, type=int)
    parser.add_argument('-s', '--cell-size', help='Estimated average size of cell, default: 1900',
                        default=1900, type=int)
    parser.add_argument('-i', '--max-iterations', help='Max segmentation iterations, default: 100',
                        default=100, type=int)
    parser.add_argument("-v", "--verbose", dest="verbose_count",
                        action="count", default=0,
                        help="increases log verbosity for each occurence.")

    args = parser.parse_args()

    # Logging
    log.setLevel(max(3 - args.verbose_count, 1) * 10)
    logging.basicConfig(format="[%(asctime)-15s] %(name)s %(levelname)s %(message)s")

    # Read image
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # ignore classic "tags are not ordered by code"
        log.debug("Reading file %s", args.image)
        image = imread(args.image, plugin='tifffile')[args.nuclei_channel]
        log.info("Read channel %d of %s, image shape is %dx%d", args.nuclei_channel, args.image, *image.shape)

    log.debug("Blurring image %s", args.image)
    image = segmentation.blur_frame(image)

    log.info("Segmenting image %s", args.image)
    segmented, _ = segmentation.mixture_model(
        image,
        learning_rate=args.learning_rate,
        var=args.variance,
        nuc_size=args.nucleus_size,
        cell_size=args.cell_size,
        max_iter=args.max_iterations,
        debug=log.getEffectiveLevel() == logging.DEBUG)
    log.debug("Segmentation done")

    log.info("Watersheding segmented image %s", args.image)
    labels = watershed(image, segmented)
    log.debug("Watershed done")

    if args.output:
        log.info("Saving labels to %s", args.output)
        mm = memmap(args.output, dtype=labels.dtype, mode='w+', shape=labels.shape)
        mm[:] = labels[:]


if __name__ == '__main__':
    main()
