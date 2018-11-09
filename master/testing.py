from PIL import Image
from PIL.TiffTags import TAGS
from lib.render import *
path = "L:\\Common\\AKN\\2-2-18 HCTP53YFP image analysis -andrei-gordon\\2-2-18 HCTp53YFP untreated vs MPS1i24w1uM 60x with z0.5uM step and 10uM range\\with z stacks\\TIFF\\2-2-18_MPS1i24w_60x.ome.tif"

import lib.tifffile as tiff

tf = tiff.TiffFile("test.tif")

print len(tf.pages)
# try:
# 	n = 0
# 	while 1:
# 		# Convert PIL image to OpenCV
# 		print(tf.pages[n])
# 		n+=1
# except EOFError:
# 	pass
