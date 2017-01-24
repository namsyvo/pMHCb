'''
Download numerical data (Last position-specific scoring matrix) from MHC Motif Viewer (http://www.cbs.dtu.dk/biotools/MHCMotifViewer/Home.html)
@Nam Sy Vo 2016
'''

import urllib
import sys
import os
import time

mhc_name = []
f = open(sys.argv[1])
for line in f:
	name = line.strip()
	mhc_name.append(name)

for dir_name in ["HLA-A", "HLA-B07_B39", "HLA-B40_B95", "HLA-C", "HLA-G", "DRB1", "DRB3", "DRB4", "DRB5"]:
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

for name in mhc_name:
	if "HLA-A" in filename:
		filename = name.replace("*", "")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/HLA-A_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "HLA-A/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "HLA-B" in filename:
		filename = name.replace("*", "")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/B07_-_B39_files/Media/" \
				+ filename + "/" + name + ".mat?disposition=download", "HLA-B07_B39/" + filename + ".txt")
			time.sleep(0.1)
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/B40_-_B95_files/Media/" \
				+ filename + "/" + name + ".mat?disposition=download", "HLA-B40_B95/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "HLA-C" in name:
		filename = name.replace("*", "")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/HLA-c_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "HLA-C/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "HLA-G" in name:
		filename = name.replace("*", "")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/HLA-g_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "HLA-G/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "DRB1" in name:
		filename = name.replace("-", "_")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/DRB1_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "DRB1/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "DRB3" in name:
		filename = name.replace("-", "_")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/DRB3_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "DRB3/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "DRB4" in name:
		filename = name.replace("-", "_")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/DRB4_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "DRB4/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)

	if "DRB5" in name:
		filename = name.replace("-", "_")
		try:
			urllib.urlretrieve ("http://www.cbs.dtu.dk/biotools/MHCMotifViewer/DRB5_files/Media/" \
				+ filename + "/" + filename + ".mat?disposition=download", "DRB5/" + filename + ".txt")
			time.sleep(0.1)
		except:
			print("Error ", sys.exc_info()[0], " for id ", filename)
