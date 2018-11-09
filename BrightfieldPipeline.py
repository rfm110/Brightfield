import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# import torch
# import torch.nn as nn
# import torchvision
import random
from PIL import Image
import skimage.io
from scipy import ndimage as ndi
from skimage.morphology import watershed as wtrsh
from skimage.feature import peak_local_max
# from skimage.filters import threshold_local
# from skimage.morphology import disk, closing, watershed
# from skimage.filters import median, rank, threshold_otsu, gaussian, threshold_local
from skimage.segmentation import random_walker
from scipy.ndimage.filters import gaussian_filter
import sklearn.preprocessing as skp
from skimage.morphology import watershed
from skimage.feature import peak_local_max

import sys
# CREDIT: Lab's Git
sys.path.insert(0, '/Users/rubabmalik786/ronglilab/Fall 2018/')
sys.path.insert(0, '/master/')
sys.path.insert(0, '/master/lib/')
sys.path.insert(0, '/master/lib/processing/')
sys.path.insert(0, '/master/lines')
# from master.lib.processing import binarize_image
from master.lib.pathfinder import *
import master.lib.read_write
from master.lines.mito_counter import imglattice2graph
import master.lib.processing
import master.lines.mito_counter

def gammaStabilize(image, alpha = 5, method = 'min'):
    bits = {'uint8':8,
            'uint16':16,
            'uint32':32
            }
    # bits = bits[image.dtype.name]
    bits = 64
    if method == 'min':
        inner = np.min(image)
        stabilized = (image - inner) / (float(2 ** bits) - inner)
        stabilized[stabilized < alpha * np.median(stabilized)] = 0
        return stabilized

def plotImage(image):
    plt.figure(figsize=(10, 7))
    plt.suptitle('Current')
    main_ax = plt.subplot(121)
    plt.title('Current')
    plt.imshow(image, interpolation='nearest', cmap=plt.cm.viridis)
    plt.show()

def plotCellAndMito(cellArray, mitoArray, binCellArray, binMitoArray):
    plt.figure(figsize=(10, 7))
    plt.suptitle('Linked Cell and Mitochondria')
    main_ax = plt.subplot(121)
    plt.subplot(2, 2, 1)
    plt.title('Cell Image')
    plt.imshow(cellArray, interpolation='nearest', cmap=plt.cm.viridis)
    plt.subplot(2, 2, 2)
    plt.title('Mitochondria Image')
    plt.imshow(mitoArray, interpolation='nearest', cmap=plt.cm.viridis)
    plt.subplot(2, 2, 3)
    plt.title('Binarized Cell Image')
    plt.imshow(binCellArray, interpolation='nearest', cmap=plt.cm.viridis)
    plt.subplot(2, 2, 4)
    plt.title('Binarized Mitochondria Image')
    plt.imshow(binMitoArray, interpolation='nearest', cmap=plt.cm.viridis)
    plt.show()
def cellBinarization(cell):
    pass
def watershedCell(image):
    distance = ndi.distance_transform_edt(image)
    localMaxi = peak_local_max(distance, indices=False)
    markers = ndi.label(localMaxi)[0]
    labels = wtrsh(-distance, markers, mask = image)
    fig, axes = plt.subplots(ncols=3, figsize=(9, 3))
    ax = axes.ravel()

    ax[0].imshow(image, cmap= plt.cm.binary, interpolation='nearest')
    ax[0].set_title('Overlapping objects')
    ax[1].imshow(-distance, cmap=plt.cm.binary, interpolation='nearest')
    ax[1].set_title('Distances')
    ax[2].imshow(labels, cmap=plt.cm.nipy_spectral, interpolation='nearest')
    ax[2].set_title('Separated objects')
    print 'hello'

    for a in ax:
        a.set_axis_off()

    fig.tight_layout()
    plt.show()

def findImages(pathname):
    print '>> in find images'
    testImagesPath = []
    for path, directory, files in os.walk(pathname):
        print path
        print directory
        print files
        for file in files:
            print file
            # find cell files (contain the word 'Brightfield')
            if file.lower().endswith('.tiff') or file.lower().endswith('.tif') and 'Brightfield' in file:
                imageNumber = file[2]
                print('image: ', imageNumber)
                print(file)
                print(str(os.path.join(path, file)))

                for altFile in files:
                    if altFile[2] == imageNumber and altFile != file:
                        testImagesPath.append([str(os.path.join(path, file)), str(os.path.join(path, altFile))])
    return testImagesPath
def binarize3D(image, layers, dimensions):
    print 'image layers: ', layers
    splitImage = np.split(image, layers, axis=0)
    print splitImage
    base = None
    counter = 0
    intensity = np.full((dimensions),10000)
    for stack in splitImage:
        stack = stack[0]
        print stack.shape

        # TODO: binarize stack here
        plotImage(stack)

        # binLayer = master.lib.processing.improved_watershed(stack, intensity)
        binLayer = watershedCell(stack)
        plotImage(binLayer)
        if counter == 1:
            base = binLayer
        else:
            base = np.concatenate((base, binLayer), axis=0)
    return base

def otsu_binar
def mainProcessing(pathname):

    # Find and import images as numpy arrays
    testImagesPath = findImages(pathname)

    # Initiate text file that will store data
    # dataAnalysis = open(os.path.join(pathname, "Cell and Mitochondria Statistics.txt"), "w")

    # Begin processing
    testImages = []
    print testImagesPath
    # TODO: remove this once finished testing
    testImagesPath.append(['C:\\Users\\rmalik9\Downloads\\WT1_w1Brightfield to Confocal.TIF','C:\\Users\\rmalik9\\Downloads\\WT1_w2561 Laser.TIF'])
    print testImagesPath

    for imagePath in testImagesPath:

        # Work simultaneously with corresponding cells and mitochondria
        cellImage = imagePath[0]
        print cellImage
        mitoImage = imagePath[1]
        cellImage2 = skimage.io.imread(cellImage)
        mitoImage2 = skimage.io.imread(mitoImage)

        imageCell = np.array(cellImage2)
        imageMito = np.array(mitoImage2)

        imageDimensionsCell = np.shape(imageCell[0])
        print 'image dimensions for cell are ',imageDimensionsCell
        print imageCell.shape
        imageMito = np.array(imageMito)
        imageDimensionsMito = np.shape(imageMito[0])
        z1 = imageCell.shape[0]
        imageCell2d = np.sum(imageCell, axis=0) // z1
        z2 = imageMito.shape[0]
        imageMito2d = np.sum(imageMito, axis=0) // z2
        # w=watershedCell(imageCell2d)


        binaryCell3D = None
        # Binarize each stack, then concatenate
        binaryCell3D = binarize3D(imageCell, imageCell.shape[0], imageDimensionsCell)
        # Convert 3D image into graph for segmentation
        item, graph = imglattice2graph(binaryCell3D)

        # Cell Smoothing - noise processing
        plotImage(imageCell2d)
        # imageCell2d = gammaStabilize(imageCell2d)
        # plotImage(imageCell2d)

        # tempBinary = imageCell2d > thresoldFilter * 0.9

        # For BFS, need to binarize each stack and then concatenate:
        for layer in imageCell:
            pass
        # Cell Binarization
        # imageCellBin = master.lib.processing.improved_watershed(imageCell2d, 1000)
        # plotCellAndMito(imageCell2d,imageMito2d,imageCellBin,imageMito)
        #  Mito Binarization

        # segmentation using layer comparator
        # master.lines.mito_counter.layer_comparator(imageCell)
        # master.lines.mito_counter.layer_comparator(imageMito)

        cellsBefore = None
        cellsAfter = None
        mitoBefore = None
        mitoAfter = None
        #
    return None

# Import images and find filenames
# pathname = '/Users/rubabmalik786/Downloads/Lab Images - BF'
pathnamegpu = 'C:\Users\rmalik9\Downloads'
print pathnamegpu
images = mainProcessing(pathnamegpu)


print images
# for path, directory, files in os.walk(pathname):
#     for file in files:
#         # find cell files (contain the word 'Brightfield')
#         if file.lower().endswith('.tiff') or file.lower().endswith('.tif') and 'Brightfield' in file:
#             imageNumber = file[2]
#             print('image: ', imageNumber)
#             print(file)
#             print(str(os.path.join(path,file)))
#
#
#             for altFile in files:
#                 if altFile[2] == imageNumber and altFile != file:
#                     testImagesPath.append([str(os.path.join(path,file)),str(os.path.join(path,altFile))])
# # print(testImagesPath)
#
#
# testImages = []
# # find matching mito images based on WTX, where X is the image number
# for imagePath in testImagesPath:
#     print(imagePath)
#     cellImage = imagePath[0]
#     mitoImage = imagePath[1]
#     cellImage2 = skimage.io.imread(cellImage)
#     mitoImage2 = skimage.io.imread(mitoImage)
#
#     imageArrayCell = [np.array(cellImage2)]
#     imageArrayMito = [np.array(mitoImage2)]
#     #     print('cell image shape: ',imageArrayCell.shape)
#     testImages.append((imageArrayCell, imageArrayMito))
#
#
# # input is 512x512xZStack
# def binarize(image, threshold=2200):
#     return skp.binarize(image, threshold=threshold)
#
#
# def plotImage(image):
#     plt.figure(figsize=(20.0, 15.0))
#     plt.suptitle('Current')
#     main_ax = plt.subplot(121)
#     plt.title('Current')
#     plt.imshow(image, interpolation='nearest', cmap=plt.cm.viridis)
#     plt.show()
#
#
# def plotCellAndMito(cellArray, mitoArray, binCellArray, binMitoArray):
#     plt.figure(figsize=(10, 7))
#     plt.suptitle('Linked Cell and Mitochondria')
#     main_ax = plt.subplot(121)
#     plt.subplot(2, 2, 1)
#     plt.title('Cell Image')
#     plt.imshow(cellArray, interpolation='nearest', cmap=plt.cm.viridis)
#     plt.subplot(2, 2, 2)
#     plt.title('Mitochondria Image')
#     plt.imshow(mitoArray, interpolation='nearest', cmap=plt.cm.viridis)
#     plt.subplot(2, 2, 3)
#     plt.title('Binarized Cell Image')
#     plt.imshow(binCellArray, interpolation='nearest', cmap=plt.cm.viridis)
#     plt.subplot(2, 2, 4)
#     plt.title('Binarized Mitochondria Image')
#     plt.imshow(binMitoArray, interpolation='nearest', cmap=plt.cm.viridis)
#     plt.show()
#
# # following function borrowed from lab's cell segmentation pipeline
# # tailored for mitochondrial binarization
# def binarize_image(base_image, _dilation = 0, feature_size = 2):
# 	'''
# 	Binarizes an image using local otsu and random walker
# 	Borrowed from Andrei's Imagepipe
#
# 	:param base_image: [np.ndarray] input image
# 	:param _dilation: [float] amount of dilation to implement in Binarization
# 	:param feature_size: [float] size of the structuring disk for random Walker
# 	:return: [np.ndarray] binarized image
# 	'''
# 	print "> Binarizing Image..."
# 	if np.percentile(base_image, 99) < 0.8:
# 		if np.percentile(base_image, 99) > 0:
# 			mult = 0.2 / np.percentile(base_image, 99)  # poissonean background assumptions
# 		else:
# 			mult = 1000. / np.sum(base_image)
# 		base_image = base_image * mult
# 		base_image[base_image > 1] = 1
#
# 	clustering_markers = np.zeros(base_image.shape, dtype=np.uint8)
# 	selem2 = disk(feature_size)
# 	print '> Performing Local Otsu'
# 	local_otsu = rank.otsu(base_image, selem2)
#
# 	# view_2d_img(local_otsu)
# 	clustering_markers[base_image < local_otsu * .9] = 1
# 	clustering_markers[base_image > local_otsu * 1.1] = 2
# 	print "> Performing Random Walker Binarization"
# 	binary_labels = random_walker(base_image, clustering_markers, beta = 10, mode = 'bf') - 1
#
# 	# if _dilation:
# 	# 	selem = disk(_dilation)
# 	# 	binary_labels = dilation(binary_labels, selem)
# 	return binary_labels
# def binarize_3d(floats_volume, cutoff):
#     """
#     Performs a 3d binarization
#
#     :param floats_volume:
#     :param cutoff:
#     :return:
#     """
#     binary_volume = np.zeros_like(floats_volume)
#     binary_volume[floats_volume > cutoff] = 1
#     return binary_volume.astype(np.bool)
#
# def yetAnotherBinarizeFunction(base):
#     value = threshold_otsu(base)
#     print value
#
#
#
# for image in testImages:
#     imageCell = image[0][0]
#
#     imageCell = np.array(imageCell)
#     imageDimensionsCell = np.shape(imageCell[0])
#     imageMito = image[1][0]
#     imageMito = np.array(imageMito)
#     imageDimensionsMito = np.shape(imageMito[0])
#
#
#
#     # Step 1: Flatten data images into a single z stack (2D)
#     #      the following may NOT HOLD FOR ALL INPUTS
#     # AVERAGE PROJECTION : CREDIT - Lab Git
#     z1 = imageCell.shape[0]
#     imageCell2d = np.sum(imageCell, axis=0) // z1
#     z2 = imageMito.shape[0]
#     imageMito2d = np.sum(imageMito, axis=0) // z2
#     # plotImage(imageMito)
#     print imageCell2d
#     print imageMito2d
#
#     # binarizedCell = binarize_image(skp.normalize(imageCell2d))
#     # binarizedMito = binarize_image(skp.normalize(imageMito2d))
#     # binarizedCell = skp.binarize(imageCell2d, threshold = 8000)
#     # binarizedMito = skp.binarize(imageMito2d, threshold = 1000)
#     binarizedCell = binarize_image(yetAnotherBinarizeFunction(imageCell2d))
#     binarizedMito = binarize_image(yetAnotherBinarizeFunction((imageMito2d)))
#     print binarizedCell
#     print binarizedMito
#
#     binarizedCell3d = None
#     binarizedMito3d = None
#
#     plotCellAndMito(imageCell2d, imageMito2d, binarizedCell, binarizedMito)
#
# # Step 2: Binarize image
#
#
# # Step 3: apply BFS (Breadth First Search)
#     cellBFS = Graph()
#     print cellBFS.get_self()
#     mitoBFS = Graph()
#     print cellBFS
#
#     break

# testImagesPath = []
#





