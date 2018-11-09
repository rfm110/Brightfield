from lib.read_write import *
import re

file_object = open("L:\\users\\gordon\\0007 - Running Projects\\20180126 Mito quantification for Gordon\\20181014 WT Mito Only\\UUID_LUT.txt", 'r')
content = file_object.readlines()
file_object.close()
test = re.sub('\\r', '\n', content[0])
womp = open("L:\\users\\gordon\\0007 - Running Projects\\20180126 Mito quantification for Gordon\\20181014 WT Mito Only\\UUID_LUT_correction.txt", 'w')
womp.write(test)
womp.close()
