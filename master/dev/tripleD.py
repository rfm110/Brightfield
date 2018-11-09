import sys
from skimage import io
import argparse
from lib.render import *
from lib.processing import *
from skimage.morphology import (disk, dilation, erosion)
from skimage.filters import rank
import scipy.signal
import cv2
# from lread_write import *

def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for analyzing 3d Images without 2d compression')
	parser.add_argument('-r',
						dest = 'read_path',
						help = 'Raw data read directory',
						required = True)
	parser.add_argument('-w',
						dest = 'write_path',
						help = 'Results save directory',
						required = False)
	options = vars(parser.parse_args())
	return options


def main(args):
	options = get_args(args)
	read_path = options['read_path']
	save_path = options['write_path']

	image = io.imread(read_path)
	test = avg_projection(image)
	view_2d_img(test)
	z_dim, x_dim, y_dim = image.shape
	for slice_index in xrange(2):
		slice_data = image[slice_index, :, :]

		noise_disk = disk(1)
		output = gamma_stabilize(slice_data, alpha_clean = 0.05)
		output = smooth_tripleD(output, smoothing_px = 1, stdevs = 1)
		# output = median(output, disk(2))
		for num_dilate in xrange(1):
			output = erosion(output, disk(1))
		for num_erode in xrange(1):
			output = dilation(output, disk(1))
		# output = cv2.blur(output, (4,4))

		montage_n_x((slice_data, output))


		# output3 = dilation(test2, disk(1))
		# output3 = dilation(output3, disk(1))
		# output3 = dilation(output3, disk(1))
		# montage_n_x((output1, output2, output3))
		# output3 = erosion(output3, disk(3))
		# output4 = median(output3, noise_disk)
		# output5 = img_type_2uint8(output4, func = 'floor')
		# asdf = fft_ifft(output2, FFT_filter)
		# montage_n_x((slice_data, output1, output2),( output3, output4, output5))
		# a = rank.autolevel(slice_data, disk(5))
		# a2 = rank.autolevel(slice_data, disk(10))
		# b = rank.enhance_contrast_percentile(slice_data, disk(5))
		# b2 = rank.enhance_contrast_percentile(slice_data, disk(10))
		# montage_n_x((slice_data, output1),(a, a2), (b, b2))
	# stack_viewer(image)





if __name__ == "__main__":
	main(sys.argv)
