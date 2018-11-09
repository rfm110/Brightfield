import sys
import argparse
from mito_counter import *
from lib.processing import *
from lib.read_write import *
from lib.render import *
from lib.UUID  import *
'''
The objective of this script is to analyze mitochondria channel images only
without consideration for cell outlines. Intended to be run alone on a single folder with mitochondria imagesself.

Program will also match genotype with each mitochondria characterized
'''

def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for 3d segmenting mitochondria')
	parser.add_argument('-r',
						dest = 'results_dir',
						help = 'Main results directory',
						required = True)
	parser.add_argument('-z',
						dest = 'z_height',
						help = 'z stack height',
						type = float,
						required = True)

	options = vars(parser.parse_args())
	return options

def main(args):
	try:
		options = get_args(sys.argv)
		save_dir = options['results_dir']
		z_height = options['z_height']
	except:
		raise Exception('Provide results directory')
	save_dir_analysis =  os.path.join(save_dir, 'analysis')
	MO_save_dir =  os.path.join(save_dir, 'analysis_MOnly')
	mkdir_check(MO_save_dir)
	mito_img_database = os.path.join(save_dir, 'mito')
	# make directory for saving 3D_segmentation
	save_dir_3D = os.path.join(save_dir, '3D_seg')
	mkdir_check(save_dir_3D)


	UUID_database = read_txt_file(os.path.join(save_dir, "UUID_LUT.txt"))

	# Generate pairs of skeletonization and binary names
	mito_filenames = get_just_filenames(mito_img_database, suffix = '.mat')
	skel_names = []
	fill_names = []
	for img_filename in mito_filenames:
		if '_skel' in img_filename.lower():
			skel_names.append(img_filename)
		elif '_bin' in img_filename.lower():
			fill_names.append(img_filename)
	skel_names.sort()
	fill_names.sort()
	paired_namelist = list(zip(fill_names, skel_names))
	mito_stats = open(os.path.join(save_dir, "mitochondria_statistics_WT.txt"), "w")
	# for each pair of images (skeleton and binary), acquire attributes
	for pair in paired_namelist:
		# Get Heat shock status for the image
		UUID_ID = get_UUID(pair[0])
		image_info = UUID2Info(UUID_ID, UUID_database)
		complete_path = image_info[-1]
		if "hs" in complete_path:
			if "_hs" in complete_path:
				temp_status = "BeforeHS"
			else:
				temp_status = "AfterHS"
		elif "rec" in complete_path:
			temp_status = "Rec"
		# Start image processing
		binary = scipy.io.loadmat(os.path.join(mito_img_database, pair[0]))['data']
		skeleton = scipy.io.loadmat(os.path.join(mito_img_database, pair[1]))['data']
		# segment 3d image
		labeled_binary = layer_comparator(binary)
		labeled_skeleton = stack_stack_multply(labeled_binary, skeleton)
		# acquire skeleton length
		unique, counts = np.unique(labeled_skeleton, return_counts = True)
		skel_length = dict(zip(unique, counts))

		# Remove background pixel count (0)
		skel_length.pop(0, None)

		for element_num in set(labeled_binary.flatten()):
			if element_num == 0:
				pass
			else:
				mask = np.zeros_like(labeled_binary)
				mask[labeled_binary == element_num] = element_num
				volume, nTriangle, SArea = get_attributes(mask, x = 1.0, y = 1.0, stack_height = z_height)
				if element_num not in skel_length:
					length = 0
				else:
					length = skel_length[element_num]
				mito_stats.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(UUID_ID,
																					element_num,
																					temp_status,
																					volume,
																					nTriangle,
																					SArea,
																					length))
		save_data(labeled_binary, UUID_ID + "_M_bin", save_dir_3D)
		save_data(labeled_skeleton, UUID_ID + "_M_skel", save_dir_3D)

	print "> Mitochondria Statistics written to {}".format(os.path.join(save_dir, "mitochondria_statistics_WT.txt"))
	mito_stats.close()




if __name__  == "__main__":
	main(sys.argv)
