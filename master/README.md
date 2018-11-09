# Bootlegged Image Analysis Pipeline

the most bootlegged image processing pipeline for mitochondrial 3d skeletonization ever

more smashed together than smashmouth
## Running the Script

To run this script, must be entered into console in the form:
```python
python main.py "<READ DIRECTORY>" "<CELL DATA WRITE DIRECTORY>" "<MITO DATA WRITE DIRECTORY>"
```
with their respective fields filled appropriately.

## Dataset Prerequisites

The read directory must contain paired images with either w1488 or w2561 included in the filename.
* w1488 corresponds to transmitted light images of yeast cells
* w2561 corresponds to 3d stack images of mitochondrial structures

Each image set should be a pair, and each pair should be located within the same folder.
There is no way for the algorithm to determine otherwise which images are paired with each other.
## Main data analysis (cell segmentation and mitochondrial segmentation)
### Data Organization
Before running the script, it would be preferable for you to make three folders in the same directory titled "cell", "mito", and "analysis"
The "cell" folder would contain your cell verification images, "mito" folder mitochondria verification images, and "analysis" paired data results
This will allow for the second module, mc_analyzer.py to be run without modification

When done running, in the write directory , there will be a file called “**filename_list.txt**” < **DO NOT DELETE THIS FILE**
This file contains the associations of each image file to their corresponding location and experiment, since each image is given a new unique hexadecimal identifier (such as 000a8b8873a7497aa839c98a2ea55d6b)

Under the folder ‘**cell**’ or the directory you specified for saving your cell data, you’ll find .mat files and .png files for each image that was binarized
Under the folder ‘**mito**’ or the directory you specified for saving your mito data, you’ll find individual skeletonizations, and stack binary images, and max projection images
### Data Output Naming
The prefix ‘**C_**’ designates cell image
* The suffix ‘**_dat.mat**’ means the file is a data file with the labeled binarized cells
* The suffix ‘**_fig.png**’ means the file is an image of the binarized cells
* The suffix ‘**_maxP_fig.png**’ means the file is the max projection of the original image

The prefix ‘**M_**’ designates mitochondrial data
* The suffix ‘**_bin.mat**’ means the file contains data on the binary mitochondria in 3d
* The suffix ‘**_skel.mat**’ means the file contains data on the skeleton of the mitochondria in 3d
* The suffix ‘**_maxP_fig.png**’ means the file is the max projection of the original image

## Analysis Function
The Analysis module is included in mc_analyzer, where individual cell and mitochondrial images are paired with their respective counterparts from data stored in the lookuptable “**filename_list.txt**.”
### Data Organization
The lookuptable containing the UUIDs must be included in the same directory as the '**cell**' and '**mito**' results folders from before.
Segmented mitochondrial skeletons will be loaded to the desired save directory, or the directory '**analysis**'

### Running the Script
```python
python mc_analyzer "<LOOKUPTABLE LOCATION>" "<WRITE DIRECTORY>"
```
### Data output naming convention
Under the folder ‘analysis’, there are three types of files you’ll see there:
* **Cell_mito_UUID_Pairs.txt** < This file links each cell image with its corresponding mitochondria image via its UUID
* **UUID_paired_filenames.txt** < This file links each binary cell image filename with its corresponding mitochondria skeleton and mito binary image filename

You’ll see ‘**.mat**’  files with names like “CM_ff751e114f834eaabda6c7640715c473_skel.mat”
* The prefix ‘**CM_**’ means that the file contains data concerning cell juxtaposed mitochondrial data
* The suffix ‘**_bin.mat**’ means the file contains labeled mitochondria with each cell
* The suffix ‘**_skel.mat**’ means that the file contains labeled mito skeletons with each cell

The UUID of each cell file is used to ID each of the '**CM__**' files, so each '**CM__**' can be traced directly to the corresponding cell image '**C_**', and the corresponding '**M_**' file can be derived from the lookuptable




Data will be written to the folder '**analysis**', or whatever save directory you specify
