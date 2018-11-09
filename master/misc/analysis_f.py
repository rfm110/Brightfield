from render import *
from processing import *

def verify_shape(img_2d, stack_3d):
	'''
	Function verifies that the shape of a 2d image matches with a single slice of a 3d stack image.
	Helper function for stack_multiplier

	:param img_2d: 2d image input
	:param stack_3d: 3d stack image input
	:return: boolean indicating whether a single slice of the 3d stack matches in dimension w/ the 2d image
	'''
	z3, x3, y3 = stack_3d.shape
	x2, y2 = img_2d.shape
	if x2 == x3 and y2 == y3:
		return True
	else:
		return False


def stack_multiplier(image, stack):
	'''
	Multiplies each layer of a 3d stack image (3d image) with a 2d image after verifying shape fit

	:param image: 2d Image to be multiplied
	:param stack: 3d stack image to have 2d image convoluted w/ along all slices
	:return: returns a convoluted 3d image
	'''
	z, x, y = stack.shape
	composite = np.zeros_like(stack)
	if verify_shape(image, stack):
		for layer in xrange(z):
			composite[layer, :, :] = stack[layer, :, :] * image
	return composite
