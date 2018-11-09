import numpy as np
import os
import sys
import math
import copy
from point import *
from skimage.exposure import adjust_gamma
from skimage import io
from scipy.ndimage import gaussian_filter
from scipy import ndimage as ndi
from skimage.morphology import disk, closing, watershed
from skimage.filters import median, rank, threshold_otsu, gaussian
from skimage.segmentation import random_walker
from skimage.restoration import denoise_bilateral
# estimate_sigma
# from skimage.feature import peak_local_max
from properties import global_max, global_min
from skimage.transform import hough_circle
# hough_circle_peaks
from skimage.draw import circle_perimeter
from math_funcs import *
from render import *
from skimage import measure
from scipy.stats import iqr
from scipy.ndimage.morphology import binary_fill_holes
from skimage.segmentation import active_contour
from skimage import measure


dtype2bits = {'uint8': 8,
			  'uint16': 16,
			  'uint32': 32}

dtype2range = { 'uint8': 255,
				'uint16': 65535,
				'uint32': 4294967295,
				'uint64': 18446744073709551615}

def gamma_stabilize(image, alpha_clean = 5, floor_method = 'min'):
	"""Normalizes the luma curve. floor intensity becomes 0 and max allowed by the bit number - 1
	Borrowed from Andrei's Imagepipe

	:param image: [np.ndarray]
	:param alpha_clean: [int] size of features that would be removed if surrounded by a majority of
	:param floor_method: [str] ['min', '1q', '5p', 'median'] method of setting the floor intensity. 1q is first quartile, 1p is the first percentile
	:return: [np.ndarray]
	"""
	bits = dtype2bits[image.dtype.name]
	if floor_method == 'min':
		inner_min = np.min(image)
	elif floor_method == '1q':
		inner_min = np.percentile(image, 25)
	elif floor_method == '5p':
		inner_min = np.percentile(image, 5)
	elif floor_method == 'median':
		inner_min = np.median(image)
	else:
		raise PipeArgError('floor_method can only be one of the three types: min, 1q, 5p or median')
	stabilized = (image - inner_min) / (float(2 ** bits) - inner_min)
	stabilized[stabilized < alpha_clean*np.median(stabilized)] = 0
	return stabilized


def sum_projection(image, axis = 0):
	'''Axis is defined as the index of the image.shape output.
	By default it is the Z axis (z,x,y)
	Takes a 3d image, and sums it along a defined axis to form a 2d Image

	:param image: [np.ndarray] 3d stack image in [np.ndarray] format
	:param axis: [int] axis to sum along, z = 0, x = 1, y = 2
	:return: [np.ndarray] returns 2d image in the shape x,y summed along z axis (by default)
	'''
	try:
		return np.sum(image, axis)
	except ValueError:
		if not((axis >= 0 or axis <= 2) and isinstance(axis, int)):
			print "Axis value invalid"
		else:
			print "Image input faulty"
		sys.exit()


def max_projection(image, axis = 0):
	'''Axis is defined as the index of the image.shape output.
	By default it is the Z axis (z,x,y)
	Takes a 3d image, and projects it along a defined axis to form a 2d Image
	Takes the max value for each pixel along (default) x,y plane, and projects it
	to one plane

	:param image: [np.ndarray] 3d stack image in <np.ndarray> format
	:param axis: [int] axis to sum along, z = 0, x = 1, y = 2
	:return: [np.ndarray] returns 2d image in the shape x,y
	'''
	try:
		return np.amax(image, axis)
	except ValueError:
		if not((axis >= 0 or axis <= 2) and isinstance(axis, int)):
			print "Axis value invalid"
		else:
			print "Image input faulty"
		sys.exit()


def avg_projection(image, axis = 0):
	'''Axis is defined as the index of the image.shape output.
	By default it is the Z axis (z,x,y)
	Takes a 3d image, and projects it along a defined axis to form a 2d Image
	Takes the average value for each pixel in the (default) x,y plane along the
	z axis

	:param image: [np.ndarray] 3d stack image in <np.ndarray> format
	:param axis: [int] axis to sum along, z = 0, x = 1, y = 2
	:return: [np.ndarray] returns 2d image in the shape x,y
	'''
	try:
		# print axis
		z, x, y = image.shape
		return np.sum(image, axis)//z
	except ValueError:
		if not((axis >= 0 or axis <= 2) and isinstance(axis, int)):
			print "Axis value invalid"
		else:
			print "Image input faulty"
		sys.exit()


def disk_hole(image, radius, pinhole = False):
	'''Returns either an image of a pinhole or a circle in the
	middle for removing high frequency/ low frequency noise using FFT
	same dimensions as input image
	Pinhole =  True : pinhole filter, or high pass filter. Filters out low frequency
	content to yield edges
	Pinhole = False: single dot filter, preserves low frequency content

	:param image: [np.ndarray] 2d input image (filter will be applied eventually), used to get dims
	:param radius: [int] radius of pinhole/pinpoint
	:param pinhole: [bool] determines whether the filter will be a pinhole or pinpoint
	:return: [np.ndarray] 2d filter of same size of 2d image input
	'''
	x, y = image.shape
	structuring_element = np.zeros((x, y), dtype = long)
	center = x//2

	for rows in xrange(x):
		for cols in xrange(y):
			if (rows - center + 0.5) ** 2 + (cols - center + 0.5) ** 2 <= radius ** 2:
				structuring_element[rows, cols] = 1
	if pinhole:
		return 1 - structuring_element
	else:
		return structuring_element


def smooth(image, smoothing_px = 0.5, threshold = 1):
	"""Gaussian smoothing of the image
	Borrowed from Andrei's Imagepipe

	:param image: [np.ndarray] Input image
	:param smoothing_px: [float] size of smoothing pixel
	:param threshold: [int] threshold to filter image intensity
	:return: [np.ndarray]
	"""
	# print "> Filtering image with Gaussian filter"
	if len(image.shape) > 2:
		for i in range(0, image.shape[0]):
			image[i, :, :] = gaussian_filter(image[i, :, :],
											 smoothing_px, mode='constant')
			image[image < threshold * np.mean(image)] = 0
	else:
		image = gaussian_filter(image, smoothing_px, mode='constant')
		image[image < threshold * np.mean(image)] = 0
	return image


def smooth_tripleD(image, smoothing_px = 0.5, stdevs = 1):
	"""Gaussian smoothing of the image
	Borrowed from Andrei's Imagepipe

	:param image: [np.ndarray] 3d Input image
	:param smoothing_px: [float] size of smoothing pixel
	:param threshold: [float] threshold to filter image intensity
	:return: [np.ndarray]
	"""
	# print "> Filtering image with Gaussian filter"
	if len(image.shape) > 2:
		for i in range(0, image.shape[0]):
			image[i, :, :] = gaussian_filter(image[i, :, :],
											 smoothing_px, mode='constant')
			image[image < mean + stdev * stdevs] = 0
	else:
		image = gaussian_filter(image, smoothing_px, mode='constant')
		mean, stdev = px_hist_stats_n0(image)
		image[image < mean + stdev * stdevs] = 0
	return image


def fft_ifft(image, struct_element):
	'''Performs a fast fourier transform, removes certain frequencies highlighted by
	the structuring element, and returns the inverse fourier transform back.
	Helper function disk_hole

	:param image: [np.ndarray] Image to be filtered
	:param struct_element: [np.ndarray] filter to be applied to image in frequency space, should be same dimension as input image
	:return: [np.ndarray] filtered image
	'''
	# print "> Performing FFT>filter>IFFT transform"

	fft_transform = np.fft.fft2(image)
	f_shift = np.fft.fftshift(fft_transform)

	# magnitude_spectrum = 20*np.log(np.abs(f_shift))
	# view_2d_img(magnitude_spectrum)
	# view_2d_img(struct_element)

	f_shift_filtered = f_shift * struct_element

	# magnitude_spectrum_filtered = 20*np.log(np.abs(f_shift_filtered))
	# view_2d_img(magnitude_spectrum_filtered)

	f_inv_shift = np.fft.ifftshift(f_shift_filtered)
	recovered_img = np.fft.ifft2(f_inv_shift)
	recovered_img = np.abs(recovered_img)
	return recovered_img


def bandpass_disk(image, r_range = (10, 200), pinhole = False):
	'''
	Creates a bandpass disk for filtering FFT images, creates either a solid
	torus filter or negative image of that

	:param image: [np.ndarray] image that filter will be applied to
	:param r_range: [tuple] (inner, outer) radius of bandpass filter
	:param pinhole: [bool] torus (true) or inverse torus (false) filter
	:return: [np.ndarray] bandpass filter structuring element
	'''
	outer = disk_hole(image, r_range[1], pinhole)
	inner = disk_hole(image, r_range[0], pinhole)
	structuring_element = outer - inner
	return structuring_element


def median_layers(image, struct_disk_r = 5):
	'''
	Applies median filter over multiple layer of a 3d stack image

	:param image: [np.ndarray] input image
	:param struct_disk_r: [float] size of median filter kernel
	:return: [np.ndarray] de-noised median filtered image
	'''
	for i in range(0, image.shape[0]):
		image[i, :, :] = median(image[i, :, :], disk(struct_disk_r))
	return image


def img_type_2uint8(base_image, func = 'floor'):
	'''
	Converts a given image type to a uint8 image
	Rounding is done either via 'floor', 'ceiling', or 'fix' functions in numpy

	:param base_image: [np.ndarray] input image
	:param func: [str] function used for scaling image pixel intensity
	:return: [np.ndarray] uint8 image
	'''
	# print "> Converting Image to uin8"
	try:
		bi_max_val = global_max(base_image)
		bi_min_val = global_min(base_image)
		dt_max = dtype2range['uint8']
		dt_min = 0

		# scaled = dt_min * (1 - ((base_image - bi_min_val) / (bi_max_val - bi_min_val))) + dt_max * ((base_image - bi_min_val)/(bi_max_val - bi_min_val))
		scaled = (base_image - bi_min_val) * ((dt_max - dt_min) / (bi_max_val - bi_min_val)) + dt_min

		if func == 'floor':
			pre_int = np.floor(scaled)
		elif func == 'ceiling':
			pre_int = np.ceil(scaled)
		elif func == 'fix':
			pre_int = np.fix(scaled)
		else:
			raise IOError

		return np.uint8(pre_int)
	except IOError:
		print "Function '{}' not recognized ".format(func)
		sys.exit()


def binarize_image(base_image, _dilation = 0, feature_size = 2):
	'''
	Binarizes an image using local otsu and random walker
	Borrowed from Andrei's Imagepipe

	:param base_image: [np.ndarray] input image
	:param _dilation: [float] amount of dilation to implement in Binarization
	:param feature_size: [float] size of the structuring disk for random Walker
	:return: [np.ndarray] binarized image
	'''
	print "> Binarizing Image..."
	if np.percentile(base_image, 99) < 0.20:
		if np.percentile(base_image, 99) > 0:
			mult = 0.20 / np.percentile(base_image, 99)  # poissonean background assumptions
		else:
			mult = 1000. / np.sum(base_image)
		base_image = base_image * mult
		base_image[base_image > 1] = 1

	clustering_markers = np.zeros(base_image.shape, dtype=np.uint8)
	selem2 = disk(feature_size)
	print '> Performing Local Otsu'
	local_otsu = rank.otsu(base_image, selem2)
	# view_2d_img(local_otsu)
	clustering_markers[base_image < local_otsu * 0.9] = 1
	clustering_markers[base_image > local_otsu * 1.1] = 2
	print "> Performing Random Walker Binarization"
	binary_labels = random_walker(base_image, clustering_markers, beta = 10, mode = 'bf') - 1

	if _dilation:
		selem = disk(_dilation)
		binary_labels = dilation(binary_labels, selem)
	return binary_labels


def hough_num_circles(input_binary_img, min_r = 15, max_r = 35, step = 2):
	'''
	Helper function for cell_split
	Runs hough transform on a cell cluster in a binary image to determine where
	individual cells may lie. Does not take in a whole binary image, only takes in a contour converted
	to a binary image for a single cluster

	:param input_binary_img: [np.ndarray] Input binary image of the single cell group
	:param min_r: [float] minimum radius acceptable for a cell
	:param max_r: [float] maximum radius acceptable for a cell
	:param step: [float] rate at which minimum radius will be stepped up to maximum radius size
	:return: [np.ndarray] cropped and split version of input binary image
	'''
	print "> Performing Hough cell splitting"
	# Create a list of radii to test and perform hough transform to recover circle centers (x,y) and radii
	hough_radii = np.arange(min_r, max_r, 2)
	hough_res = hough_circle(input_binary_img, hough_radii)
	accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks = 3)
	circles = zip(cy, cx, radii);
	# remove any circles too close to each other

	no_duplicates = crop_close(circles, max_sep = 10)

	# HYPOTHETICAL # of cells
	N_cells = len(no_duplicates)
	# view_2d_img(input_binary_img)
	print "\t> Number cells in subsection: {}".format(N_cells)
	# print no_duplicates
	if N_cells > 1:
		# Set initial radius to 1
		for rows in no_duplicates:
			rows[-1] == 1
		# Create mask to divide both cells
		# Grow circle size until there is a collision followed by no more collisions
		actual_mask = np.zeros_like(input_binary_img)
		# Create Conditions for Whileloop
		collision = False
		end_collision = False
		stop_condition = False
		n = 0
		max_iter = 100

		# while end_collision == False or n < 10:
		while (n < max_iter and stop_condition == False):
			# Create empty mask to grow
			mask = np.zeros_like(input_binary_img)
			for center_y, center_x, radius in no_duplicates:
				# Create mask for each circle in no_duplicates
				sub_mask = np.zeros_like(input_binary_img)
				# List of coodinates for perimeter of a circle and remove any negative points to prevent edge looparound
				circy, circx = circle_perimeter(center_y, center_x, radius)
				no_negs = remove_neg_pts(zip(circy, circx))
				for y, x in no_negs:
					# for each circle point, if it falls within the dimensions of the submask, plot it.
					if y < sub_mask.shape[0] and x < sub_mask.shape[1]:
						sub_mask[y, x] = 1
				# Append submask to growing mask after dilation (aka making the circle boundaries wide as fuck)
				mask += dilation(sub_mask, disk(2))
			# Determine if there is a collision between circles (max element in grow mask is more than just one submask)
			# montage_n_x((actual_mask, mask))
			# print np.amax(mask.flatten())
			if np.amax(mask.flatten()) > 1:
				collision = True
				# collision_pt = np.where(mask >= 2)
				mask[mask < 2] = 0
				# view_2d_img(mask)
				actual_mask += mask
				actual_mask[actual_mask != 0 ] = 1
			# montage_n_x((actual_mask, mask))
			if collision == True and np.amax(mask.flatten()) <= 1:
				end_collision = True
				stop_condition = True

			# Grow circle radius by 8% per
			for rows in no_duplicates:
				rows[-1] *= 1.08
				rows[-1] = np.int(rows[-1])
			n += 1
			# print n, collision, end_collision, stop_condition
		if stop_condition == False:
			# montage_n_x((actual_mask, actual_mask + input_binary_img, filled_cells, filled_cells * (1 - actual_mask)))
			return np.ones_like(input_binary_img)
			# return binary_fill_holes(input_binary_img).astype(int)
		else:
			# Fill edges to create mask
			# filled_cells = binary_fill_holes(input_binary_img).astype(int)
			# montage_n_x((actual_mask, actual_mask + input_binary_img, filled_cells, filled_cells * (1 - actual_mask)))
			# view_2d_img(filled_cells * (1 - dm))
			# return filled_cells * (1 - actual_mask)
			return actual_mask
	else:
		# view_2d_img(input_binary_img)
		return np.ones_like(input_binary_img)
		# return binary_fill_holes(input_binary_img).astype(int)

	# Uncomment for visualization
	# montage_n_x((input_binary_img, filled_cells, dm,  filled_cells * (1 - dm)))
	# for center_y, center_x, radius in no_duplicates:
	# 	circy, circx = circle_perimeter(center_y, center_x, radius)
	# 	input_binary_img[circy, circx] = 10
	# view_2d_img(input_binary_img)


def just_label(binary_image):
	'''
	Just labels a binary image (segments everything)

	:params binary_image: [np.ndarray] input image
	:return: [np.ndarray] segmented image
	'''
	labeled_field, object_no = ndi.label(binary_image, structure = np.ones((3, 3)))
	return labeled_field


def label_and_correct(binary_channel, pre_binary, min_px_radius = 10, max_px_radius = 100, min_intensity = 0, mean_diff = 10):
	"""
	Labelling of a binary image, with constraints on minimal feature size, minimal intensity of area
	 covered by a binary label or minimal mean difference from background

	:param binary_channel: [np.ndarray] input image
	:param pre_binary: [np.ndarray] used to compute total intensity
	:param min_px_radius: [float] minimal feature size
	:param min_intensity: [float] minimal total intensity
	:param mean_diff: [float] minimal (multiplicative) difference from the background
	:return: [np.ndarray]
	"""
	labeled_field, object_no = ndi.label(binary_channel, structure = np.ones((3, 3)))

	# prebinary_px = pre_binary.flatten()
	# n_bins = int(2 * iqr(prebinary_px) * (len(prebinary_px) ** (1/3)))
	# n, bin_edge = np.histogram(prebinary_px, n_bins)
	# peak_max_indx = np.argmax(n)
	# background_val = (bin_edge[peak_max_indx] + bin_edge[peak_max_indx + 1]) / 2
	background_mean = np.mean(pre_binary[labeled_field == 0])
	for label in range(1, object_no+1):
		mask = labeled_field == label
		px_radius = np.sqrt(np.sum((mask).astype(np.int8)))
		total_intensity = np.sum(pre_binary[mask])
		label_mean = np.mean(pre_binary[labeled_field == label])
		if px_radius < min_px_radius or total_intensity < min_intensity or label_mean < mean_diff * background_mean or px_radius > max_px_radius:
			labeled_field[labeled_field == label] = 0
	# dbg.label_and_correct_debug(labeled_field)
	return labeled_field


def cell_split(input_img, contours, min_area = 100, max_area = 3500, min_peri = 100, max_peri = 1500):
	'''
	Function finds individual cluster of cells that have a contour fall within the bounds of area and perimeter
	and attempts to divide them by their constituents

	:param input_img: [np.ndarray] binary input image containing all segmented cells
	:param contours: [list] of [lists] a list of list of points for each of the contours for each cell
	:param min_area: [float] minimum acceptable area for a cell
	:param max_area: [float] max acceptable area for a cell
	:param min_peri: [float] minimum acceptable perimeter for a cell
	:param max_peri: [float] maximum acceptable perimeter for a cell
	:return: [np.ndarray] full image with cell clusters split
	'''
	print "> Starting Cell Split"
	output = np.zeros_like(input_img)
	output[input_img > 0] = 1
	for item_contour in contours:
		# remove cells that have a low circumference or too high circumference
		contour_len = item_contour.shape[0]
		contour_img = binary_fill_holes(points2img(item_contour)).astype(int)
		contour_area = sum(contour_img.flatten())
		if contour_len >= min_peri and contour_len <= max_peri and contour_area >= min_area and contour_area <= max_area:
		# if item_contour.shape[0] >= 100 and item_contour.shape[0] <= 350:
			# holding = points2img(item_contour)
			# holding_fill = binary_fill_holes(holding).astype(int)
			# if sum(holding_fill.flatten()) > 100:
			split_cells_mask = hough_num_circles(contour_img)

			tlx, tly, brx, bry = location(item_contour)
			if array_all_ones(split_cells_mask):
				for x in xrange(tlx, brx + 1):
					for y in xrange(tly, bry + 1):
						output[y, x] = output[y, x] * split_cells_mask[y - tly, x - tlx]
				# output[tly - 1:bry + 1, tlx - 1:brx + 1] = output[tly - 1:bry + 1, tlx - 1:brx + 1] * split_cells_mask
				# montage_n_x((output,output[tly-1:bry+1,tlx-1:brx+1], split_cells_mask, output[tly-1:bry+1,tlx-1:brx+1] * split_cells_mask))
			else:
				d_contour_img = dilation(contour_img, disk(1))
				specific_mask = split_cells_mask + d_contour_img
				specific_mask[specific_mask < 2] = 0
				specific_mask[specific_mask >= 2] = 1
				for x in xrange(tlx, brx + 1):
					for y in xrange(tly, bry + 1):
						output[y, x] = output[y, x] * (1- specific_mask[y - tly, x - tlx])


	return label_and_correct(output, input_img)


def rm_eccentric(input_img, min_eccentricity, max_ecc, min_area, max_area):
	'''Evaluates the eccentricity of single cells within an image with multiple cells, and throws away any cells that exhibit odd eccentricity
	Also chucks any cells that have an area larger than max_area

	:param input_img: [np.ndarray] segmented binary image
	:param min_eccentricity: [float] minimum acceptable eccentricity
	:param max_area: [float] maximum area acceptable for a cell
	'''

	max_cells = np.amax(input_img.flatten())
	output_img = copy.deepcopy(input_img)
	for x in xrange(max_cells):
		mask = np.zeros_like(input_img)
		mask[output_img == x + 1] = 1

		_, area, _, eccentricity = get_contour_details(mask)

		if eccentricity < min_eccentricity or eccentricity > max_ecc and area < min_area or area > max_area:
			output_img[output_img == x + 1] = 0
	return output_img


def improved_watershed(binary_base, intensity, expected_separation = 10):
	"""
	Improved watershed method that takes in account minimum intensity as well as minimal size of
	separation between the elements
	Borrowed from Andrei's Imagepipe

	:param binary_base: [np.ndarray] support for watershedding
	:param intensity: [np.ndarray] intensity value used to exclude  watershed points with too low of intensity
	:param expected_separation: [float] expected minimal separation (in pixels) between watershed centers
	:return: [np.ndarray]
	"""
	print "> Performing Improved Watershed"
	# sel_elem = disk(2)
	#
	# # changed variable name for "labels"
	# post_closing_labels = closing(binary_base, sel_elem)

	distance = ndi.distance_transform_edt(binary_base)
	local_maxi = peak_local_max(distance,
								indices = False,  # we want the image mask, not peak position
								min_distance = expected_separation,  # about half of a bud with our size
								threshold_abs = 10,  # allows to clear the noise
								labels = binary_base)
	# we fuse the labels that are close together that escaped the min distance in local_maxi
	local_maxi = ndi.convolve(local_maxi, np.ones((5, 5)), mode = 'constant', cval = 0.0)
	# finish the watershed
	struct = np.ones((3, 3))
	struct[0,0] = 0
	struct[0,2] = 0
	struct[2,0] = 0
	struct[2,2] = 0
	expanded_maxi_markers = ndi.label(local_maxi, structure = struct)[0]
	segmented_cells_labels = watershed(-distance, expanded_maxi_markers, mask = binary_base)

	unique_segmented_cells_labels = np.unique(segmented_cells_labels)
	unique_segmented_cells_labels = unique_segmented_cells_labels[1:]
	average_apply_mask_list = []
	# Gimick fix
	# intensity_array = intensity * np.ones_like(segmented_cells_labels)
	for cell_label in unique_segmented_cells_labels:
		my_mask = segmented_cells_labels == cell_label
		apply_mask = segmented_cells_labels[my_mask]
		average_apply_mask = np.mean(intensity[my_mask])
		if average_apply_mask < 0.005:
			average_apply_mask = 0
			segmented_cells_labels[segmented_cells_labels == cell_label] = 0
		average_apply_mask_list.append(average_apply_mask)
	# x_labels = ['cell13', 'cell1', 'cell7', 'cell2', 'cell14', 'cell6', 'cell3', 'cell5', 'cell4', 'cell11', 'cell12', 'cell8', 'cell10', 'cell9']
	# dbg.improved_watershed_debug(segmented_cells_labels, intensity)
	# dbg.improved_watershed_plot_intensities(x_labels, average_apply_mask_list.sort())
	return segmented_cells_labels


def get_contour_details(input_img):
	'''
	Input image must have only one cell in it with one contour, function extracts
	area, perimeter, radius, eccentricity from a given cell. radius is an
	approximation, assuming the cell is circular
	:param input_img: input image (binary) of just a single cell. Cannot contain
						multiple cells or multiple contours
	:return: radius, area, perimeter, and eccentricity in that order
	'''
	contours = measure.find_contours(input_img,
										level = 0.5,
										fully_connected = 'low',
										positive_orientation = 'low')
	Point_set = Point_set2D(contours[0])
	radius = (Point_set.shoelace() / math.pi) ** 0.5
	eccentricity = (4 * math.pi * Point_set.shoelace()) / (Point_set.perimeter() ** 2)
	return radius, Point_set.shoelace(), Point_set.perimeter(), eccentricity


def hist_peak(image):
	'''
	Returns the peak histogram count for a given image

	:param image: 2d-3d image for generating histogram map
	:return: peak histogram value
	'''
	px_dataset = image.flatten()
	n_bins = int(2 * iqr(px_dataset) * (len(px_dataset)) ** (1/3))
	n, bin_edges = np.histogram(px_dataset, n_bins)
	peak_max_indx = np.argmax(n)
	return (bin_edges[peak_max_indx] + bin_edges[peak_max_indx]) / 2

def write_stats(before_image, after_image, UID, filename, read_path, write_path, img_type = "cell"):
	'''
	Given two segmented binary images, determine the difference between the
	cells present on both images and save differences and cell stats to a file
	:param before_image: [np.ndarray] Image before cell deletion or addition
	:param after_image: [np.ndarray] Image after cell deletion or addition
	:param UID: [str] unique UID for the image, should be the same between the two Images
	:param filename: [str] name of the datafile to be written to
	:param write_path: [str] location of the datafile containing data
	'''
	write_file = open(os.path.join(write_path, filename),'a')
	# write_file.write("cell\tUID\tcell_num\tcell_updated_num\tdeleted\tRadius\tArea\tPerimeter\tEccentricity\tread_path\n")
	deletion = True
	before_cells = np.amax(before_image.flatten())
	after_cells = np.amax(after_image.flatten())
	if after_cells > before_cells:
		deletion = False
	max_cells = np.maximum(before_cells, after_cells)
	for cell_index in xrange(max_cells):
		current_label = cell_index + 1
		after_label = current_label
		if deletion:
			mask = np.zeros_like(before_image)
			mask[before_image == current_label] = current_label
		else:
			mask = np.zeros_like(after_image)
			mask[after_image == current_label] = current_label
		radius, area, perimeter, E = get_contour_details(mask)
		cell_isolation = after_image * mask

		cell_delete = False
		if np.amax(cell_isolation.flatten()) != current_label:
			if np.amax(cell_isolation.flatten()) == 0:
				after_label = 0
				cell_delete = True
			else:
				after_label = int(np.amax(cell_isolation.flatten()) / current_label)
		write_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(img_type,
																			UID,
																			current_label,
																			after_label,
																			cell_delete,
																			radius,
																			area,
																			perimeter,
																			E,
																			read_path))
	write_file.close()
