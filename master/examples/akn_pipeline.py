from skimage import io
import sys
from lib.render import *
from lib.processing import *
import lib.read_write as rw
from PIL import Image
from skimage import feature
from skimage.filters import threshold_otsu, threshold_adaptive
import numpy as np

suffix = ".tif"

def single(image_path):
	image = io.imread(image_path)
	DAPI_raw = image[:, :, 0]
	YFP_raw = image[:, :, 1]

	x, y, z = image.shape

	DAPI = gamma_stabilize(DAPI_raw, alpha_clean = 1.3)
	YFP = gamma_stabilize(YFP_raw, alpha_clean = 1.3)
	#
	DAPI = smooth(DAPI)
	YFP = smooth(YFP)


	sel_elem = disk(2)
	DAPI = median(DAPI, sel_elem)
	YFP = median(YFP, sel_elem)

	a5 = erosion(DAPI, selem = disk(1))
	a6 = median(a5, sel_elem)
	a7 = dilation(a6, selem = disk(1))
	a8 = img_type_2uint8(a7, func = 'floor')

	b5 = erosion(YFP, selem = disk(1))
	b6 = median(b5, sel_elem)
	b7 = dilation(b6, selem = disk(1))
	YFP_foreground = img_type_2uint8(b7, func = 'floor')

	a9 = binarize_image(a8, _dilation = 0, feature_size = 100)
	# a10 = binary_fill_holes(a9).astype(int)
	# a11 = label_and_correct(a10, a8, min_px_radius = 10, max_px_radius = 500, min_intensity = 0, mean_diff = 10)
	a11 = label_and_correct(a9, a8, min_px_radius = 10, max_px_radius = 500, min_intensity = 0, mean_diff = 10)
	# a111 = label_and_correct(a10, a8, min_px_radius = 10, max_px_radius = 500, min_intensity = 0, mean_diff = 15)
	segmented_binary = improved_watershed(a11, a8, expected_separation = 1)

	# Get background and segmented region values
	positive_binary = np.zeros_like(segmented_binary)
	positive_binary[segmented_binary >= 1] = 1
	negative_binary = 1 - positive_binary

	# Apply positive and negative (feature and background) masks
	DAPI_pos_seg = positive_binary * DAPI
	DAPI_neg_seg = negative_binary * DAPI
	YFP_pos_seg = positive_binary * YFP
	YFP_neg_seg = negative_binary * YFP

	return DAPI, DAPI_raw, YFP, YFP_raw, YFP_foreground, segmented_binary, a11


def get_img_stats(binary, raw, channel_name, background = False):
	internal_binary = np.zeros_like(binary)
	if background:
		internal_binary[binary == 0] = 1
		img_type = "Background"
	else:
		internal_binary[binary != 0] = 1
		img_type = "Foreground"
	num_px = np.sum(internal_binary.flatten()).astype(float)
	px_mask = np.multiply(internal_binary, raw)
	avg = np.sum(px_mask.flatten().astype(float)) / num_px
	stdev = np.std(px_mask.flatten().astype(float))
	print "\t> {} Channel {} average: {}".format(channel_name, img_type, avg)
	print "\t> {} Channel {} stdev: {}".format(channel_name, img_type, stdev)
	print "\t> {} Channel {} Pixel Count: {}".format(channel_name, img_type, num_px)
	return num_px, avg, stdev, internal_binary, px_mask


def analyze(image_path, data_save_dir, txt_save_filename):
	filename = os.path.basename(image_path)
	img_filename = filename.replace(suffix, '')
	print img_filename, filename
	# sys.exit()
	print "> Analysing {}".format(filename)

	DAPI, DAPI_raw, YFP, YFP_raw, YFP_foreground, segmented_binary, labeled_binary = single(image_path)

	write_data = open(os.path.join(data_save_dir, txt_save_filename), 'a')
	write_data.write("Results for file \"{}\"====================\n\n".format(filename))
	write_data.write("File location: \n\t\"{}\"\n".format(image_path))
	write_data.write("Population Metrics:\n")
	DAPI_num_bkgd_px, DAPI_avg_bkgd, DAPI_bkgd_std, DAPI_neg_binary, DAPI_neg_mask = get_img_stats(segmented_binary,
																									DAPI_raw,
																									channel_name = "DAPI",
																									background = True)

	DAPI_num_fgd_px, DAPI_avg_fgd, DAPI_fgd_std, DAPI_pos_binary, DAPI_pos_mask  = get_img_stats(segmented_binary,
																									DAPI_raw,
																									channel_name = "DAPI",
																									background = False)

	YFP_num_bkgd_px, YFP_avg_bkgd, YFP_bkgd_std, YFP_neg_binary, YFP_neg_mask = get_img_stats(YFP_foreground,
																								YFP_raw,
																								channel_name = "YFP",
																								background = True)

	YFP_num_fgd_px, YFP_avg_fgd, YFP_fgd_std, YFP_pos_binary, YFP_pos_mask = get_img_stats(YFP_foreground,
																							YFP_raw,
																							channel_name = "YFP",
																							background = False)

	# montage_n_x((DAPI_neg_binary, DAPI_raw, DAPI_neg_mask),
	# 			(DAPI_pos_binary, DAPI_raw, DAPI_pos_mask),
	# 			(YFP_neg_binary, YFP_raw, YFP_neg_mask),
	# 			(YFP_pos_binary, YFP_raw, YFP_pos_mask))

	write_data.write("\tForeground:\n")
	write_data.write("\t\tDAPI Pixel Count: \t{}\n".format(DAPI_num_fgd_px))
	write_data.write("\t\tDAPI Intensity Average: \t{}\n".format(DAPI_avg_fgd))
	write_data.write("\t\tDAPI Intensity StDeviation: \t{}\n".format(DAPI_fgd_std))
	write_data.write("\t\tYFP Pixel Count: \t{}\n".format(YFP_num_fgd_px))
	write_data.write("\t\tYFP Intensity Average: \t{}\n".format(YFP_avg_fgd))
	write_data.write("\t\tYFP Intensity StDeviation: \t{}\n".format(YFP_fgd_std))
	write_data.write("\tBackground:\n")
	write_data.write("\t\tDAPI Pixel Count: \t{}\n".format(DAPI_num_bkgd_px))
	write_data.write("\t\tDAPI Intensity Average: \t{}\n".format(DAPI_avg_bkgd))
	write_data.write("\t\tDAPI Intensity StDeviation: \t{}\n".format(DAPI_bkgd_std))
	write_data.write("\t\tYFP Pixel Count: \t{}\n".format(YFP_num_bkgd_px))
	write_data.write("\t\tYFP Intensity Average: \t{}\n".format(YFP_avg_bkgd))
	write_data.write("\t\tYFP Intensity StDeviation: \t{}\n".format(YFP_bkgd_std))

	rw.save_figure(DAPI_neg_binary, img_filename + "_DAPI_neg_bin.tif", data_save_dir)
	rw.save_figure(DAPI_pos_binary, img_filename + "_DAPI_pos_bin.tif", data_save_dir)
	rw.save_figure(DAPI_raw, img_filename + "_DAPI_raw.tif", data_save_dir)
	rw.save_figure(DAPI_neg_mask, img_filename + "_DAPI_neg_mask.tif", data_save_dir)
	rw.save_figure(DAPI_pos_mask, img_filename + "_DAPI_pos_mask.tif", data_save_dir)
	rw.save_figure(YFP_neg_binary, img_filename + "_YFP_neg_bin.tif", data_save_dir)
	rw.save_figure(YFP_pos_binary, img_filename + "_YFP_pos_bin.tif", data_save_dir)
	rw.save_figure(YFP_raw, img_filename + "_YFP_raw.tif", data_save_dir)
	rw.save_figure(YFP_neg_mask, img_filename + "_YFP_neg_mask.tif", data_save_dir)
	rw.save_figure(YFP_pos_mask, img_filename + "_YFP_pos_mask.tif", data_save_dir)
	rw.save_figure(segmented_binary, img_filename + "_seg_bin.tif", data_save_dir)
	rw.save_figure(labeled_binary, img_filename + "_lab_bin.tif", data_save_dir)
	rw.save_data(segmented_binary, img_filename + "_seg_bin", data_save_dir)
	rw.save_data(labeled_binary, img_filename + "_lab_bin", data_save_dir)

	write_data.write("\nCell Resolution Metrics:\n")
	num_cell_frags = np.amax(segmented_binary.flatten())
	num_connected_component = np.amax(labeled_binary.flatten())
	print num_cell_frags
	print num_connected_component
	# view_2d_img(segmented_binary)
	write_data.write("\tNumber of Cell Fragments: {}\n".format(num_cell_frags))
	write_data.write("\tNumber of Cell Connected Components: {}\n".format(num_connected_component))

	write_data.write("\tSegmented Component Attributes: \n")
	write_data.write("\t\tCell ID \t DAPI_avg \t DAPI_std \t YFP_avg \t YFP_std \n")

	for fragment in xrange(1, num_cell_frags + 1):
		mask = np.zeros_like(segmented_binary)
		mask[segmented_binary == fragment] = 1
		mask_area = np.sum(mask.flatten().astype(float))
		if mask_area == 0:
			pass
		else:
			DAPI_masked = np.multiply(mask, DAPI_raw)
			YFP_masked = np.multiply(mask, YFP_raw)
			DAPI_avg = np.sum(DAPI_masked.flatten().astype(float)) / mask_area
			DAPI_std = np.std(DAPI_masked.flatten().astype(float))
			YFP_avg = np.sum(YFP_masked.flatten().astype(float)) / mask_area
			YFP_std = np.std(YFP_masked.flatten().astype(float))
			write_data.write("\t\t {} \t {} \t {} \t {} \t {} \n".format(fragment, DAPI_avg, DAPI_std, YFP_avg, YFP_std))

	write_data.write("\tConnected Component Attributes: \n")
	write_data.write("\t\tCell ID \t DAPI_avg \t DAPI_std \t YFP_avg \t YFP_std \n")

	for cc_fragment in xrange(1, num_connected_component + 1):
		mask = np.zeros_like(labeled_binary)
		mask[labeled_binary == cc_fragment] = 1
		mask_area = np.sum(mask.flatten().astype(float))
		if mask_area == 0:
			pass
		else:
			DAPI_masked = np.multiply(mask, DAPI_raw)
			YFP_masked = np.multiply(mask, YFP_raw)
			DAPI_avg = np.sum(DAPI_masked.flatten().astype(float)) / mask_area
			DAPI_std = np.std(DAPI_masked.flatten().astype(float))
			YFP_avg = np.sum(YFP_masked.flatten().astype(float)) / mask_area
			YFP_std = np.std(YFP_masked.flatten().astype(float))
			write_data.write("\t\t {} \t {} \t {} \t {} \t {} \n".format(cc_fragment, DAPI_avg, DAPI_std, YFP_avg, YFP_std))


	write_data.close()



if __name__	 == "__main__":
	read_dir = "L:\\Common\\AKN\\2-2-18 HCTP53YFP image analysis -andrei-gordon\\2-2-18 HCTp53YFP untreated vs MPS1i24w1uM 60x with z0.5uM step and 10uM range\\single images\\TIFF\\"
	save_dir = "L:\\Common\\AKN\\2-2-18 HCTP53YFP image analysis -andrei-gordon\\2-2-18 HCTp53YFP untreated vs MPS1i24w1uM 60x with z0.5uM step and 10uM range\\single images\\TIFF\\"
	filelist = rw.get_img_filenames(read_dir)
	for image in filelist:
		analyze(image[-1], save_dir, "holla.txt")
	# print "asdf"
	# analyze("L:\\Common\\AKN\\2-2-18 HCTP53YFP image analysis -andrei-gordon\\2-2-18 HCTp53YFP untreated vs MPS1i24w1uM 60x with z0.5uM step and 10uM range\\single images\\TIFF\\2-2-18_MPS1i24w_60x.tif")
