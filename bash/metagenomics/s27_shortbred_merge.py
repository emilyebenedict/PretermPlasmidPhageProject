#!/usr/bin/env python3
"""
Merges shortBRED output files together
NOTE: Run this script in the same directory as the shortbred outputs that you want to merge! Make sure that there aren't any other files in that directory 
other than this script or else the script will error out!
Usage: python3 shortbred_merge.py
"""

import os

samples=[]

for sample in os.listdir():
	if ".py" not in sample:
		sample = sample.split("_ShortBRED")[0]
		samples.append(sample)
print(samples)
joint_output = []
joint_output.append("Sample\tFamily\tCount\tHits\tTotMarkerLength")

paths = []
for sample in samples:
	path = os.getcwd() + "/" + sample + "_ShortBRED.txt"
	paths.append(path)

string_to_match = "shortbred_merge.py_ShortBRED.txt"
# delete path in the list paths that matches string_to_match because this is wrongly added to the list in line 20
paths = [item for item in paths if string_to_match not in item]

for path in paths:	
	for line in open(path, "r"):
		if not line.startswith("Family"):
			entry = sample + "\t" + line.strip()
			joint_output.append(entry)

joint_output_file = open(os.getcwd() + "/joint_output2.txt", "w")

for line in joint_output:
	print(line, file = joint_output_file)

joint_output_file.close()
