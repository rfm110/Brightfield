import sys
from render import *
from processing import *
from math_funcs import *
from properties import properties
from read_write import *
from skimage import io
from skimage import measure
from scipy.ndimage.morphology import binary_fill_holes
from skimage.morphology import (disk, dilation, watershed,
								closing, opening, erosion, skeletonize, medial_axis)
from skimage.morphology import skeletonize_3d
from skimage.filters import threshold_local
from scipy.stats import iqr

mito_prefix = "M_"
skeleton_suffix = "_skel"
binary_suffix = "_bin"
figure_suffix = "_fig.png"

''' Purpose of script is to segment and binarize 3d mitochondrion images.
Runs in line with main.py up one directory (under main project path)
'''
def analyze(UID, read_path, write_path):
	# Script is designed to handle 512 by 512 sized images
	mito = io.imread(read_path)
	z, x, y = mito.shape
	binary = np.zeros_like(mito)
	max_P_d = avg_projection(mito)
	# Get average background px intensity from average projection
	px_dataset = max_P_d.flatten()
	# n_bins = int(2 * iqr(px_dataset) * (len(px_dataset) ** (1/3)))
	# n, bin_boundary = np.histogram(px_dataset, n_bins)
	# hist_peak_max_indx = np.argmax(n)
	# # threshold for the background
	# bin_midpt = (bin_boundary[hist_peak_max_indx] + bin_boundary[hist_peak_max_indx + 1]) / 2
	for layer in xrange(z):
		# clean up any background noise
		sel_elem = disk(1)
		layer_data = mito[layer,:,:]

		output1 = gamma_stabilize(layer_data,
									alpha_clean = 1,
									floor_method = 'min')
		output2 = smooth(output1,
							smoothing_px = 2,
							threshold = 1)
		median_filtered = median(output1, sel_elem)
		fft_filter_disk = disk_hole(median_filtered,
									radius = 5,
									pinhole = True)
		# Remove High frequency noise from image
		FFT_Filtered = fft_ifft(median_filtered, fft_filter_disk)
		# Convert image to 8 bit for median filter to work
		image_8bit = img_type_2uint8(FFT_Filtered, func = 'floor')
		test = median(image_8bit, sel_elem)

		test_px_dataset = test.flatten()
		n_bins = int(2 * iqr(px_dataset) * (len(px_dataset) ** (1/3)))
		n, bin_edge = np.histogram(test_px_dataset, n_bins)
		test_peak_max_indx = np.argmax(n)
		bin_midpt = (bin_edge[test_peak_max_indx] + bin_edge[test_peak_max_indx + 1]) / 2
		test_mask = test > bin_midpt
		test_masked = test * test_mask

		local_thresh = threshold_local(test_masked,
										block_size = 31,
										offset = -15)
		binary_local = test_masked > local_thresh
		# label individual elements and remove really small noise and background
		corrected_slice = label_and_correct(binary_local, test,
												min_px_radius = 1,
												max_px_radius = 100,
												min_intensity = 0,
												mean_diff = 11.9)
		corrected_slice[corrected_slice > 0] = 1
		binary[layer, :, :] = corrected_slice

	print 'OK\n'
	spooky = skeletonize_3d(binary)
	binary_projection = max_projection(binary)

	# montage_n_x((binary_projection, max_P_d))
	# stack_viewer(binary)
	# raise Exception
	save_figure(max_P_d, mito_prefix + UID + "_maxP" + figure_suffix, write_path)
	save_figure(binary_projection, mito_prefix + UID + "_maxPB" + figure_suffix, write_path)
	save_data(spooky, mito_prefix + UID + skeleton_suffix, write_path)
	save_data(binary, mito_prefix + UID + binary_suffix, write_path)

if __name__ == "__main__":
	# Run as standalone on single image
	arguments = sys.argv
	analyze(arguments[-1])
