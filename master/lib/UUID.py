import os


def get_UUID(filename):
	'''
	Given a string filename in the form X_UUID_XXX.XXX
	Return the UUID element in the filename
	:param filename: [str] of the filename in question
	:return: UUID [str] fraction of the input filename
	'''
	split_filename = filename.split('_')
	return split_filename[1]


def get_file_datatype(filename):
	'''
	Given a string filename in the form X_UUID_XXX.TYPE
	Return the TYPE element in the filename
	:param filename: [str] of the filename in question
	:return: type [str] fraction of the input filename
	'''
	split_filename = filename.split('_')
	return split_filename[-1]


def get_M_C(filename):
	'''
	Given a string filename in the form MC_UUID_XXX.TYPE
	Return the MC element in the filename
	:param filename: [str] of the filename in question
	:return: M or C [str] fraction of the input filename, single letter string
	'''
	split_filename = filename.split('_')
	return split_filename[0]


def filename2info(filename, lookupTable):
	'''
	Given a filename (MC_UUID_TYPE), retrieve the UUID and look for its presence in the LUT
	:param filename: [str] complete filename in format MC_UUID_TYPE
	:param lookupTable: [list] LUT to search in
	:return: [list] if UUID is found, returns the row that the element was found in, does not return anything if nothing is found
	'''
	UUID = get_UUID(filename)
	found = False
	for row in lookupTable:
		if UUID == row[0]:
			found = True
			return row
	if found == False:
		print "UUID ({}, {}) not found in LUT".format(UUID, filename)


def UUID2Info(UUID, LUT):
	'''
	Given a UUID,  look for its presence in the LUT
	:param filename: [str] complete UUID
	:param lookupTable: [list] LUT to search in
	:return: [list] if UUID is found, returns the row that the element was found in, does not return anything if nothing is found
	'''
	found = False
	for row in LUT:
		if UUID == row[0]:
			found = True
			return row
	if found == False:
		print "UUID ({}, {}) not found in LUT".format(UUID, filename)


def get_partner(filename, lookupTable):
	'''
	Given a UUID for a cell or mitochondria image, determine the UUID for the mitochondria or cell data partner
	:param filename: [str] original filename of the file being queried (full filename)
	:param lookupTable: [list] table to look up query for the filename
	'''
	attribute_type = get_M_C(filename)
	input_info = filename2info(filename, lookupTable)
	if input_info == None:
		return
	else:
		input_fname = input_info[2]
		shared_loc = input_info[-2]
		if attribute_type == "M":
			partner_fname = input_fname.replace("2561", "1488")
		elif attribute_type == "C":
			partner_fname = input_fname.replace("1488", "2561")
		for row_num, row in enumerate(lookupTable):
			if partner_fname == row[2] and shared_loc == row[-2]:
				return row


def create_pairTable(filelist, lookupTable, save_dir):
	'''
	given a list of files, create a pairing between M and C file counterparts based on data from lookupTable
	:param filelist: list of files with their locations and filenames
	:param lookupTable: [list] LUT to be queried to find matches
	:param save_dir: [str] location to save the pair table
	:return: [list] list of UUID pairs with the partner information, also returns condensed version (UUID_pairs_no_info)
	'''
	UUID_pairs = []
	UUID_pairs_no_info = []
	for row in filelist:
		try:
			input_info = filename2info(row, lookupTable)
			partner_info = get_partner(row, lookupTable)
			UUID_pairs.append([input_info[0], partner_info[0], input_info[2], partner_info[2], input_info[-2]])
			UUID_pairs_no_info.append([input_info[0], partner_info[0]])
		except:
			print "> Data not found in LUT"
			print row
			unfounds = open(os.path.join(save_dir,"Unfound_Imgs.txt"), 'a')
			unfounds.write(row+"\n")
			unfounds.close()
	return UUID_pairs, UUID_pairs_no_info
