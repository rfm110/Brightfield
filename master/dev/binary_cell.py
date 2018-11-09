
from lib.render import *
import numpy as np

from skimage import measure
import numpy as np
import math

def euclid_dist_nD(p0, p1):
	return np.sum((p1 - p0) ** 2)


class Point_set2D(object):
	def __init__(self, point_list):
		self.point_list = np.array([[float(coordinate) for coordinate in point] for point in point_list])


	def num_pts(self):
		return len(self.point_list)


	def center_mass(self):
		return np.sum(self.point_list, 0) / self.num_pts()


	def perimeter(self):
		peri_distance = 0
		for pt_indx in xrange(self.num_pts()):
			peri_distance += euclid_dist_nD(self.point_list[pt_indx],
											self.point_list[pt_indx - 1])
		return peri_distance


	def shoelace(self):
		area = 0
		for pt in xrange(len(self.point_list)):
			area += self.point_list[pt - 1][0] * self.point_list[pt][1]
			area -= self.point_list[pt - 1][1] * self.point_list[pt][0]
		return abs(area) / 2.0


class Binary_cell_img(object):
	def __init__(self, single_cell_img):
		self.single_cell_img = single_cell_img
		if self.num_contours() != 1:
			print "> Too many contours detected in object"

	def img_shape(self):
		return self.single_cell_img.shape


	def contours(self):
		return measure.find_contours(self.single_cell_img,
										level = 0,
										fully_connected = 'low',
										positive_orientation = 'low')


	def num_contours(self):
		return len(self.contours())


	def contour_perimeter(self):
		peri_table = []
		for contour_num in xrange(self.num_contours()):
			contour = Point_set2D(self.contours()[contour_num])
			peri_table.append(contour.perimeter())
		return np.array(peri_table)


	def cell_area(self):
		return np.sum(self.single_cell_img.flatten())



points = [(0,0),
		(0,10),
		(10,10),
		(10,0)]
# print q.cell_area()
# view_2d_img(q.return_peri())
rec = Point_set2D(points)
print rec.shoelace()
