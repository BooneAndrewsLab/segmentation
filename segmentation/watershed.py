from skimage.measure import label
from skimage.morphology import closing, square
from skimage.segmentation import watershed as skwatershed


def watershed(image, segmentation):
    """ Run watershed on our cells to get labels

    :param image: The channel which shows cell outline
    :param segmentation: Segmented image, output of mixture_model
    :type image: ndarray
    :type segmentation: ndarray
    :return: labeled segmented cells, not filtered for anything
    :rtype: ndarray
    """
    mask = closing(segmentation > 0, square(5))  # try closing some gaps in segmentation
    markers = label(segmentation == 2)  # segmented nuclei are our markers/seeds
    labeled = skwatershed(-image, markers=markers, mask=mask)
    return labeled
