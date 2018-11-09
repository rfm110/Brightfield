
import sys
import os, errno
import argparse
import shutil
sys.path.insert(0, '.\\lib')
sys.path.insert(0, '.\\lines')

import scipy.io
from skimage import io
import time, string

from lib.UUID import *
from lib.render import *
from lib.read_write import *
from lib.processing import *
import lines.process_cell_channel as cell_line
import lines.process_mito_channel as mito_line
from lines import mito_counter

def blockPrint():
	sys.stdout = open(os.devnull, 'w')


def enablePrint():
	sys.stdout = sys.__stdout__


def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for analyzing mitochondria skeletonization')
	parser.add_argument('-r',
						dest = 'read_dir',
						help = 'Raw data read directory',
						required = False,
						default = ".\\testing_environment\\")
	parser.add_argument('-w',
						dest = 'save_dir',
						help = 'Save directory for segmentation and skeletonization data',
						required = False,
						default = ".\\testing_environment\\test_run\\")

	options = vars(parser.parse_args())
	if os.path.exists(options['save_dir']):
		if options['save_dir'] == ".\\test_run":
			shutil.rmtree(options['save_dir'])
	return options


def main(args):
	os.system('cls' if os.name == 'nt' else 'clear')
	options = get_args(args)
	root_read_dir = options['read_dir']
	save_dir = options['save_dir']

	print "> Parent Read Directory : {}\r".format(root_read_dir)
	print "> Save Directory : {}\r".format(save_dir)

	save_dir_cell = os.path.join(save_dir, 'cell')
	save_dir_mito = os.path.join(save_dir, 'mito')
	save_dir_anal = os.path.join(save_dir, 'analysis')

	mkdir_check(save_dir_cell)
	mkdir_check(save_dir_mito)
	mkdir_check(save_dir_anal)

	start = time.time()
	filenames = get_img_filenames(root_read_dir)
	num_images = len(filenames)
	end = time.time()
	print "> {} images detected, time taken: {}".format(num_images, end - start)

	mito_stats = []
	cell_stats = []

	print "> Processing IDs saved here: {}\r".format(save_dir)

	file_list_ID = open(os.path.join(save_dir, "UUID_LUT.txt"),'w')
	for UID, img_name, img_fname, path_diff, img_loc, img_path in filenames:
		file_list_ID.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(UID, img_name, img_fname, path_diff, img_loc, img_path))
	file_list_ID.close()


	img_num = 1
	for UID, img_name, _, _, _, img_path in filenames:
		print "> ==========================================================================================\r"
		print "\r> Currently Processing : {}\r".format(img_name)
		print "> \tImage Unique ID: {}\r".format(UID)
		print "> \tImage/Total Number of Images: {}/{}\r".format(img_num, num_images)
		if '1488' in img_name:
			# continue
			print "> Image ID: 1488 - Cell TD\r"
			blockPrint()
			start = time.time()
			cell_line.analyze(UID, img_path, save_dir_cell)
			end = time.time()
			enablePrint()
			print "> Time to Compete: {}".format(end - start)
			mito_stats.append(end - start)
			img_num += 1

		elif '2561' in img_name:
			# continue
			print "> Image ID: 2561 - Mitochondria\r"
			blockPrint()
			start = time.time()
			mito_line.analyze(UID, img_path, save_dir_mito)
			end = time.time()
			enablePrint()
			print "> Time to Compete: {}".format(end - start)
			cell_stats.append(end - start)
			img_num += 1

	print "> ==========================================================================================\r"
	print "> Prelim Analysis completed"
	save_data(mito_stats, "mito_processing_RT", save_dir)
	save_data(cell_stats, "cell_processing_RT",  save_dir)

	# Start merge of MC_analyzer
	cell_filelist = get_just_filenames(save_dir_cell, suffix = '_dat.mat')
	mito_filelist = get_just_filenames(save_dir_mito, suffix = '.mat')

	UUID_datatable = read_txt_file(os.path.join(save_dir, "UUID_LUT.txt"))

	C_M_UUID_pairs, UUID_pairs = create_pairTable(cell_filelist, UUID_datatable, save_dir)

	filename_pairs = []
	for cell_UUID, mito_UUID in UUID_pairs:
		filename_pairs.append(["C_" + cell_UUID + "_dat.mat",
								"M_" + mito_UUID + "_bin.mat",
								"M_" + mito_UUID + "_skel.mat"])
	print "> Creating UUID Filename Pairs"
	write_list_txt(save_dir_anal, "Cell_mito_UUID_Pairs.txt", C_M_UUID_pairs)
	write_list_txt(save_dir_anal, "UUID_paired_filenames.txt", filename_pairs)

	for filenames in filename_pairs:
		save_fileID = get_UUID(filenames[0])
		cell_img = scipy.io.loadmat(os.path.join(save_dir_cell, filenames[0]))['data']
		mito_stack = scipy.io.loadmat(os.path.join(save_dir_mito, filenames[1]))['data']
		mito_skel = scipy.io.loadmat(os.path.join(save_dir_mito, filenames[2]))['data']

		labeled_mito_bin = stack_multiplier(cell_img, mito_stack)
		labeled_mito_skel = stack_multiplier(cell_img, mito_skel)

		save_data(labeled_mito_bin, "CM_" + save_fileID + "_bin", save_dir_anal)
		save_data(labeled_mito_skel, "CM_" + save_fileID + "_skel", save_dir_anal)

	# Start merge of mitocounter
	mito_counter.main(save_dir)
if __name__ == "__main__":
	main(sys.argv)
