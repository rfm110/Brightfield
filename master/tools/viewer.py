from lib.render import *
import scipy.io
import sys
import argparse
import os
from lib.processing import avg_projection
from lib.read_write import get_img_filenames

'''
Given a path directory to a single image, script displays image on screen (2d and 3d compatible)
'''
def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for image visualization')
	parser.add_argument('-i', dest = 'image', help = 'Image path', required = True)
	options = vars(parser.parse_args())
	return options


def show_img(image):
	if len(image.shape) == 2:
		view_2d_img(image)
	elif len(image.shape) == 3:
		stack_viewer(image)
		avg_p = avg_projection(image)
		view_2d_img(avg_p)
	else:
		print "Too many dimensions to resolve"


def main(args):
	options = get_args(args)
	path = options['image']
	dir_read = os.path.isdir(path)
	# print dir_read
	if not dir_read:
		image = scipy.io.loadmat(path)['data']
		show_img(image)
	else:
		filenames = get_img_filenames(path, suffix = '.mat')
		num_images = len(filenames)
		print "> {} Images detected".format(num_images)
		for _, _, _, _, _, img_path in filenames:
			print "> Image Name: {}".format(img_path)
			image = scipy.io.loadmat(img_path)['data']
			show_img(image)





if __name__ == "__main__":
	main(sys.argv)
