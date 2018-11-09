from skimage import io
from lib.read_write import *
from lib.render import *
from lib.processing import *
from lib.properties import *
from lib.math_funcs import *
import numpy as np
from skimage.morphology import disk
from skimage import measure
from dev.binary_cell import *
import copy
import math
def main():
	read_directory = "L:\\Common\\Gordon Jobs\\20180222 Josh Colony Size Analysis\\cropped_dataset"
	file_list = get_img_filenames(read_directory)

	results = open(os.path.join(read_directory, "results.txt"),'w')


	for img_file in file_list:
		print img_file[1]
		raw_image = io.imread(img_file[-1])[:, :, 1]

		mean_intensity = np.mean(raw_image[raw_image > 0].flatten())
		stdev_intensity = np.std(raw_image[raw_image > 0].flatten())
		print mean_intensity, stdev_intensity
		# Simple threshold by mean + stdev
		processed_image = copy.deepcopy(raw_image)
		processed_image[processed_image < (mean_intensity + stdev_intensity)] = 0
		processed_image = gamma_stabilize(processed_image, alpha_clean = 5, floor_method = "min")
		processed_image = smooth(processed_image, smoothing_px = 0.5, threshold = 1)

		structuring_element = disk(1)
		processed_image = median(processed_image, structuring_element)
		processed_image = erosion(processed_image, selem = disk(1))
		processed_image = median(processed_image, structuring_element)
		processed_image = dilation(processed_image, selem = disk(1))
		processed_image0 = img_type_2uint8(processed_image, func = 'floor')
		# view_2d_img(processed_image0)
		processed_image1 = binarize_image(processed_image0, _dilation = 0, feature_size = 2)
		processed_image2 = label_and_correct(processed_image1, raw_image, min_px_radius = 1, max_px_radius = 100, min_intensity = 0, mean_diff = 1)
		num_colonies = global_max(processed_image2)

		for colony_num in xrange(num_colonies):
			mask = np.zeros_like(processed_image2)
			mask[processed_image2 == colony_num + 1] = 1
			contour = measure.find_contours(mask, level = 0.5, fully_connected = 'low', positive_orientation = 'low')
			if len(contour) == 1:
				contour_data = Point_set2D(contour[0])
				E = (4 * math.pi * contour_data.shoelace()) / ((contour_data.perimeter()) ** 2)
				radius = np.sqrt(contour_data.shoelace() / math.pi)
				results.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(img_file[1],
															colony_num,
															contour_data.perimeter(),
															contour_data.shoelace(),
															E,
															radius))
		save_figure(processed_image2, img_file[1] + "bin.png", read_directory)
	results.close()




if __name__ == "__main__":
	main()
