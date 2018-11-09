import sys
from skimage import io
import argparse
from lib.render import *
from lib.processing import max_projection, avg_projection, sum_projection
from lib.read_write import *

'''
Given a directory path, traverses directory and finds all images, then displays their max , average, and sum projections for comparison
'''

def get_args(args):
	parser = argparse.ArgumentParser(description = 'Script for analyzing 3d Images without 2d compression')
	parser.add_argument('-r',
						dest = 'read_path',
						help = 'Raw data read directory',
						required = True)
	options = vars(parser.parse_args())
	return options


def main(args):
	options = get_args(args)
	read_path = options['read_path']
	filenames = get_img_filenames(read_path, suffix = '.mat')
	num_images = len(filenames)
	print "> {} Images detected".format(num_images)
	for _, _, _, _, _, img_path in filenames:
		# print img_path
		image = io.imread(img_path)
		max_p = max_projection(image)
		avg_p = avg_projection(image)
		sum_p = sum_projection(image)
		montage_n_x((max_p, avg_p, sum_p))

if __name__ == "__main__":
	main(sys.argv)
