
import os.path

data_folder = os.path.join("source_data", "text_files")

file_to_open = os.path.join(data_folder, "raw_data.txt")

f = open(file_to_open)

print(f.read())
