import sys
import argparse
from render import *
from processing import *
from math_funcs import *
from properties import properties
from read_write import *
from skimage import io
from skimage import measure
from scipy.ndimage.morphology import binary_fill_holes, binary_opening
from skimage.morphology import (disk, dilation, watershed,
								closing, opening, erosion, medial_axis)
import matplotlib.pyplot as plt
from skimage.filters import threshold_local, gaussian
from skimage.segmentation import active_contour
cell_prefix = "C_"
data_suffix = "_dat"
projection_suffix = "_avgP"
figure_suffix = ".png"

''' Purpose of script is to segment and binarize 3d cell images.
Cell images must be full of cytoplasmic fluorescent protein - detects cells based
on cell internal protein fluorescence.
Runs in line with main.py up one directory (under main project path)
'''

def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for 3d segmenting cells')
	parser.add_argument('-id',
						dest = 'UUID',
						help = 'Unique identifier (user generated)',
						required = True)
	parser.add_argument('-r',
						dest = 'read_path',
						help = 'read directory for data',
						required = True)
	parser.add_argument('-w',
						dest = 'write_path',
						help = 'write directory for results',
						required = True)

	options = vars(parser.parse_args())
	return options


def analyze(UID, read_path, write_path):
	cell = io.imread(read_path)
	noise_disk = disk(1)
	img_projection = avg_projection(cell)
	noise_processing = gamma_stabilize(img_projection, alpha_clean = 1.3)
	for i in xrange(1):
		# Remove noise
		noise_processing = smooth(noise_processing,
									smoothing_px = 0.5,
									threshold = 1)

	noise_processing = img_type_2uint8(noise_processing, func = 'floor')
	noise_processing = erosion(noise_processing, noise_disk)
	noise_processing = median(noise_processing, noise_disk)
	# Simple Thresholding
	smoothed_image = smooth(noise_processing, smoothing_px = 4, threshold = 1)
	threshold_filter = threshold_local(noise_processing,
									block_size = 31,
									offset = 0)
	# mean, stdev = px_stats(noise_processing)
	# noise_processing[noise_processing < mean + stdev] = 0
	# d5 = img_type_2uint8(noise_processing, func = 'floor')

	temp_binary = smoothed_image > threshold_filter * 0.9
	improved_binary = improved_watershed(temp_binary, smoothed_image, expected_separation = 1)
	# Attempt to smooth contours on the sides of cells
	smoothed_binary = np.zeros_like(improved_binary)
	for cell_label in range(1, np.max(improved_binary) + 1):
		mask = improved_binary == cell_label
		mask = gaussian(mask, 2) > 0.5
		smoothed_binary += mask * cell_label
	smoothed_binary = smoothed_binary > 0
	improved_binary2 = improved_watershed(smoothed_binary, smoothed_image, expected_separation = 1)
	removed_eccentric = rm_eccentric(improved_binary2,
						min_eccentricity = 0.62,
						max_ecc = 0.99,
						min_area = 500,
						max_area = 2500)

	write_stats(improved_binary2,
				removed_eccentric,
				UID,
				"single_cell_stats.txt",
				read_path,
				write_path)
	# montage_n_x((img_projection, noise_processing, temp_binary, improved_binary, removed_eccentric))
	save_data(img_projection, cell_prefix + UID + projection_suffix, write_path)
	save_data(removed_eccentric, cell_prefix + UID + data_suffix, write_path)
	save_figure(img_projection, cell_prefix + UID + "_maxP" + figure_suffix, write_path)
	save_figure(improved_binary2, cell_prefix + UID + "_binary" + figure_suffix, write_path)
	save_figure(removed_eccentric, cell_prefix + UID + "_rmECC" + figure_suffix, write_path)
	# return a16


if __name__ == "__main__":
	try:
		options = get_args(sys.argv)
		read_path = options['read_path']
		write_path = options['write_path']
		UID = options['UUID']
		analyze(sys.argv)
	except:
		raise Exception('Check Inputs into command line')
