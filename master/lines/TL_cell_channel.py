import argparse
import sys
from skimage import feature
sys.path.insert(0, '.\\lib')
from render import *
from processing import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import cv2
from skimage.filters import median, gaussian, threshold_local,  roberts, sobel, scharr
from skimage.morphology import disk, watershed, dilation

''' The purpose of this script is to analyze Brightfield images z stack and segement them
for analysis.
'''

def get_args():
	parser = argparse.ArgumentParser(description = 'Script for processing transmitted cell outline images')
	parser.add_argument('-r',
						dest = 'img_path',
						help = 'location of the image',
						required = True)
	parser.add_argument('-w',
						dest = 'write_path',
						help = 'write directory for results',
						required = True)

	return parser.parse_args()


def analyze():
    args = get_args()
    cell_img = io.imread(args.img_path)
    cell_projection = max_projection(cell_img).astype(np.float64)
    threshold_filter = threshold_local(cell_projection, block_size = 31, offset = 0)
    mean, stdev = px_stats(threshold_filter)
    canny_mask = threshold_filter > mean + stdev

    struct_element = disk(15)
    wam = dilation(canny_mask, struct_element)



    edges1 = feature.canny(cell_projection, sigma = 5, mask = wam, low_threshold=0, high_threshold=10)
    edges2 = feature.canny(cell_projection, sigma = 3, mask = wam)
    montage_n_x((cell_projection, edges1, edges2, threshold_filter),
                (canny_mask * cell_projection, canny_mask, wam, wam * cell_projection))
    raise Exception
    # display results
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(8, 3),
                                        sharex=True, sharey=True)

    ax1.imshow(cell_projection, cmap=plt.cm.gray)
    ax1.axis('off')
    ax1.set_title('noisy image', fontsize=20)

    ax2.imshow(edges1, cmap=plt.cm.gray)
    ax2.axis('off')
    ax2.set_title('Canny filter, $\sigma=1$', fontsize=20)

    ax3.imshow(edges2, cmap=plt.cm.gray)
    ax3.axis('off')
    ax3.set_title('Canny filter, $\sigma=3$', fontsize=20)

    fig.tight_layout()

    plt.show()

if __name__ == "__main__":
    analyze()
