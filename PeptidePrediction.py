'''
Predicting peptides from mass specs data and PSSM data of HLA alleles
Input: mass specs data, PSSM data of HLA alleles, amino acids and their mass
Output: peptide sequences which most likely bind to HLA alleles
@Nam Sy Vo 2016

'''
import PeptideGeneration

import os
import sys
import itertools
import collections

#Analyzing PSSM (position-specific scoring matrix) data
#Data obtained from MHC Motif Viewer(http://www.cbs.dtu.dk/biotools/MHCMotifViewer/Home.html)
def getPSSM(path):
        aa_names = []
	PSSM_HLA_1 = {}
	for dir_name in ["HLA-A", "HLA-B", "HLA-C", "HLA-G"]:
		for (dirpath, dirnames, filenames) in os.walk(os.path.join(path, dir_name)):
			for filename in filenames:
				f = open(os.path.join(dirpath, filename))
				hla_name = f.readline().strip()
				f.readline()
				f.readline()
				aa_names = f.readline().strip().split()
				if len(aa_names) != 20:
					continue
				PSSM_HLA_1[hla_name] = []
				for line in f:
					tmp = line.strip().split()
					PSSM_HLA_1[hla_name].append(map(float, tmp[2:]))
        PSSM_HLA_2 = {}
	for dir_name in ["DRB1", "DRB3", "DRB4", "DRB5"]:
		for (dirpath, dirnames, filenames) in os.walk(os.path.join(path, dir_name)):
			for filename in filenames:
				f = open(os.path.join(dirpath, filename))
				hla_name = f.readline().strip()
				f.readline()
				f.readline()
				aa_names = f.readline().strip().split()
				if len(aa_names) != 20:
					continue
				PSSM_HLA_2[hla_name] = []
				for line in f:
					tmp = line.strip().split()
					PSSM_HLA_2[hla_name].append(map(float, tmp[2:]))
	return aa_names, PSSM_HLA_1, PSSM_HLA_2

#Parsing mass spec data in MGF format
def parseMGF(mgf_file):
	mgf_file = open(mgf_file)
	pep_mass = []
	for line in mgf_file:
		if "PEPMASS" in line:
			tmp = line.strip().split("=")
			pep_mass.append(float(tmp[1].split(" ")[0]))

	sorted_pep_mass = sorted(pep_mass, key=float)
	uni_pep_mass = [sorted_pep_mass[0]]
	for m in sorted_pep_mass:
		if (m - uni_pep_mass[-1])/m > pow(10, -5):
			uni_pep_mass.append(m)
	return uni_pep_mass

#Predicting peptides based on pssm for each HLA allele
def predict_hla_bound_peptides(pssm_path, top_num, pep_seq, pep_mass):

	pep_seq_set = {}
	for k, v in pep_seq.iteritems():
		pep_set = []
		for p in v:
			pep_set.append(set(p))
		pep_seq_set[k] = pep_set

	aa_names, PSSM_HLA_1, PSSM_HLA_2 = getPSSM(pssm_path)
	PSSM_HLA = PSSM_HLA_1.copy()
	PSSM_HLA.update(PSSM_HLA_2)
	print "# of HLA alleles:", len(PSSM_HLA)

	pep_db = {}
	print "indexing pssm..."
	for hla_name, hla_pssm in PSSM_HLA.iteritems():
		top_aa = []
		top_aa_set = set()
		path_num = 1
		for pos_val in hla_pssm:
			pssm_dict = {}
			for i in range(len(pos_val)):
				pssm_dict[aa_names[i]] = pos_val[i]
			top_val, top_name = [], []
			i = 0
			for k in sorted(pssm_dict, key=pssm_dict.get, reverse=True):
				i += 1
				if i > top_num or pssm_dict[k] <= 0: #take top_num AAs but with positive scores
					break
				top_val.append(pssm_dict[k])
				top_name.append(k)
				top_aa_set.add(k)
			top_aa.append(top_name)
			path_num *= (i-1)

		print "#of paths of ", hla_name, ":", path_num

		count = 0
		for t in itertools.product(*top_aa):
			aa_col = collections.Counter(t)
			aa_score = ""
			for aa in aa_names:
				aa_score += str(aa_col.get(aa, 0))
				if aa_score not in pep_db:
					pep_db[aa_score] = []
				pep_db[aa_score].append([hla_name, t])
			count += 1
			if count % 100000 == 0:
				print "processed", count, "paths"

	print "# of items in pssm index:", len(pep_db)
	'''
	cand_pep_seq = {}
	for k, v in pep_seq_set.iteritems():
		cand_pep = []
		for i in range(len(v)):
			if v[i].issubset(top_aa_set):
				cand_pep.append(pep_seq[k][i])
		if len(cand_pep) != 0:
			cand_pep_seq[k] = cand_pep
	if len(cand_pep_seq) == 0:
		continue

	print hla_name, top_aa
	for k, v in cand_pep_seq.iteritems():
		print k, len(v), len(pep_seq[k]), v

	bound_pep_seq = {}
	for k, v in cand_pep_seq.iteritems():
		bound_pep = []
		for pep in v:
			aa_col = collections.Counter(pep)
			aa_score = ""
			for aa in aa_names:
				aa_score += str(aa_col.get(aa, 0))
			if aa_score in pep_db:
				bound_pep.append([pep, hla_name, pep_db[aa_score]])
        if len(bound_pep) != 0:
			if k not in bound_pep_seq:
				bound_pep_seq[k] = []
			bound_pep_seq[k].append(bound_pep)

	if len(bound_pep_seq) == 0:
		continue

	for k, v in bound_pep_seq.iteritems():
		f = open(os.path.join("predicted_peptides_bound_top" + str(top_num), str(k_mer) + "-mer_" + str(k) + "-mass_bound.txt"), "a")
		for bound_pep in v:
			for pep in bound_pep:
				f.write(pep[0] + "\t" + pep[1] +"\t" + str(pep[2]) + "\n")
		f.close()
	'''

#Main program
if __name__ == "__main__":

	print "reading input data..."
	uni_pep_mass = parseMGF(sys.argv[1])
	print "# of unique items in mass data", len(uni_pep_mass)

	aaName, aaMass = [], []
	f = open(sys.argv[2])
	for line in f:
		tokens = line.strip().split("\t")
		aaName.append(tokens[0])
		aaMass.append(float(tokens[1]))
	f.close()

	pssm_path = sys.argv[3]
	top_num = int(sys.argv[4])
	if not os.path.exists("predicted_peptides_bound_top" + str(top_num)):
		os.makedirs("predicted_peptides_bound_top" + str(top_num))

	k_mass = int(sys.argv[5])
	ppm = 10.0
	k_mer = 9
	pep_seq, pep_mass = {}, {}

	'''
	# Test
	for pm in [385.23495, 402.22269, 443.23782]:
		if pm == 385.23495:
			pep_seq[pm] = ["KIVGAGPGA"]
			pep_mass[pm] = ["769.46263"]
		if pm == 402.22269:
			pep_seq[pm] = ["LGGSGSGLR"]
			pep_mass[pm] = ["803.4381"]
		if pm == 443.23782:
			pep_seq[pm] = ["AVATEAPNL"]
			pep_mass[pm] = ["885.46837"]

		f = open(os.path.join("predicted_peptides_top" + str(top_num), str(k_mer) + "-mer_" + str(pm) + "-mass_predicted.txt"), "w")
		f.write("predicted_peptides\taa_mass\n")
		for i in range(len(pep_seq[pm])):
			f.write(pep_seq[pm][i] + "\t" + pep_mass[pm][i] + "\n")
		f.close()

	bound_pep_seq = predict_hla_bound_peptides(pssm_path, top_num, pep_seq, pep_mass)
	'''

	print "generating peptide sequences..."
	count = 0
	for pm in uni_pep_mass:
	#for pm in uni_pep_mass[70*k_mass:70*k_mass+70]:
	#for pm in uni_pep_mass[70*k_mass:]:
		mass = 2*pm - 1.007825 #charge=2, and exclusive H+? need to clarify...
		pep_seq_list, pep_mass_list = [], []
		solPep, solMass = PeptideGeneration.generate_peptides(mass, ppm, aaName, aaMass)
		if k_mer not in solPep:
			continue
		count += len(solPep[k_mer])
		for i in range(len(solPep[k_mer])):
			pep_seq_list.append(solPep[k_mer][i])
			pep_mass_list.append(solMass[k_mer][i])
		pep_seq[pm] = pep_seq_list
		pep_mass[pm] = pep_mass_list
		print pm, len(pep_seq[pm])
		'''
		f = open(os.path.join("predicted_peptides_bound_top" + str(top_num), str(k_mer) + "-mer_" + str(pm) + "-mass_predicted.txt"), "w")
		f.write("peptides\taa_mass\n")
		for i in range(len(pep_seq[pm])):
			f.write(pep_seq[pm][i] + "\t" + pep_mass[pm][i] + "\n")
		f.close()
		'''
	print "# of generated peptide sequences:", count

	print "searching for bound peptides..."
	#bound_pep_seq = predict_hla_bound_peptides(pssm_path, top_num, pep_seq, pep_mass)
