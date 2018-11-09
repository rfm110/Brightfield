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
