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
cmp_list = lambda x, y: collections.Counter(x) == collections.Counter(y)

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

        hla_num = 0
	for hla_name, hla_pssm in PSSM_HLA.iteritems():
		if hla_name not in ["#HLA-B4001", "#HLA-B4004", "#HLA-B4004", "#HLA-B4011", "#HLA-B4014", \
		"#HLA-B4033", "#HLA-B4038", "#HLA-B4042", "#HLA-B4054", "#HLA-B4055", "#HLA-B4062", "#HLA-B4064", \
		"#HLA-B4064", "#HLA-B4065", "#HLA-B4066", "#HLA-B4067", "#HLA-B4069", "#HLA-B4074", "#HLA-B4076", "#HLA-B4077"]:
			continue
		top_aa = []
		top_aa_set = set()
		for pos_val in hla_pssm:
			pssm_dict = {}
			for i in range(len(pos_val)):
				pssm_dict[aa_names[i]] = pos_val[i]
			top_val, top_name = [], []
			i = 0
			for k in sorted(pssm_dict, key=pssm_dict.get, reverse=True):
				i += 1
				#if i > top_num or pssm_dict[k] <= 0: #take top_num AAs but with positive scores
				if i > top_num: #take top_num AAs
					break
				top_val.append(pssm_dict[k])
				top_name.append(k)
				top_aa_set.add(k)
			top_aa.append(top_name)

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

		f = open(os.path.join("test_mass_top" + str(top_num), str(k_mer) + "-mer_log.txt"), "a")
		hla_num += 1
		f.write(str(hla_num) + "\t" + str(hla_name) + "\t" + str(top_aa) + "\n")
		print "hla", hla_num, hla_name
		for k, v in cand_pep_seq.iteritems():
			f.write(str(k) + "\t" + str(len(v)) + "\t" + str(len(pep_seq[k])) + "\t" + str(v) + "\n")
			print "cand", k, len(v), len(pep_seq[k])
		f.close()

		bound_pep_seq = {}
		for t in itertools.product(*top_aa):
			for k, v in cand_pep_seq.iteritems():
				bound_pep = []
				for pep in v:
					if cmp_list(pep, t):
						bound_pep.append([pep, hla_name, t])
				if len(bound_pep) != 0:
					if k not in bound_pep_seq:
						bound_pep_seq[k] = []
						bound_pep_seq[k].append(bound_pep)
		if len(bound_pep_seq) == 0:
			continue
		for k, v in bound_pep_seq.iteritems():
			f = open(os.path.join("test_mass_top" + str(top_num), str(k_mer) + "-mer_" + str(k) + "-mass_bound.txt"), "a")
			for bound_pep in v:
				for pep in bound_pep:
					f.write(pep[0] + "\t" + pep[1] +"\t" + str(pep[2]) + "\n")
					print "bound", pep[0], pep[1], pep[2]
			f.close()

#Main program
if __name__ == "__main__":

	uni_pep_mass = parseMGF(sys.argv[1])

	aaName, aaMass = [], []
	f = open(sys.argv[2])
	for line in f:
		tokens = line.strip().split("\t")
		aaName.append(tokens[0])
		aaMass.append(float(tokens[1]))
	f.close()

	pssm_path = sys.argv[3]
	top_num = int(sys.argv[4])
	if not os.path.exists("test_mass_top" + str(top_num)):
		os.makedirs("test_mass_top" + str(top_num))

	k_mass = int(sys.argv[5])
	ppm = 10.0
	k_mer = 9
	pep_seq, pep_mass = {}, {}

	# Test
	seq_list = ["KIVGAGPGA","LGGSGSGLR","AVATEAPNL","IGIAPLAQL","GLLGTLVQL","DGGNESDPM","ALASHLIEA","AIVDKVPSV","SLLDKIIGA","SLLGGNIRL","ALLDSAHLL","SLLEKSLGL","AVLTELRAV","TLIEDILGV","LSAEKIQAL","KLIANNTTV","ISRALVTTL","SLGLPQDVPG","SLFPGKLEV"]
	mz_list = ["385.23495","402.22269","443.23782","448.28342","457.28906","461.16895","462.76172","464.27966","465.28659","471.7916","476.77765","480.29266","486.29782","486.78546","486.79044","487.28766","487.30121","491.76447","495.28641"]
	mass_list = ["769.46263","803.4381","885.46837","895.55956","913.57085","921.33061","924.51616","927.55205","929.5659","942.57591","952.54802","959.57805","971.58837","972.56365","972.57359","973.56804","973.59514","982.52165","989.56554"]
	#for i in range(len(seq_list)):
	for i in [0]:
		pep_seq[mz_list[i]] = [seq_list[i]]
		pep_mass[mz_list[i]] = [mass_list[i]]

		f = open(os.path.join("test_mass_top" + str(top_num), str(k_mer) + "-mer_" + str(mz_list[i]) + "-mass_predicted.txt"), "w")
		f.write("predicted_peptides\taa_mass\n")
		for j in range(len(pep_seq[mz_list[i]])):
			f.write(pep_seq[mz_list[i]][j] + "\t" + pep_mass[mz_list[i]][j] + "\n")
		f.close()

	bound_pep_seq = predict_hla_bound_peptides(pssm_path, top_num, pep_seq, pep_mass)

	'''
	for pm in uni_pep_mass[:1]:
	#for pm in uni_pep_mass[70*k_mass:70*k_mass+70]:
	#for pm in uni_pep_mass[70*k_mass:]:
		mass = 2*pm - 1.007825 #charge=2, and exclusive H+? need to clarify...
		pep_seq_list, pep_mass_list = [], []
		solPep, solMass = PeptideGeneration.generate_peptides(mass, ppm, aaName, aaMass)
		if k_mer not in solPep:
			continue
		for i in range(len(solPep[k_mer])):
			pep_seq_list.append(solPep[k_mer][i])
			pep_mass_list.append(solMass[k_mer][i])
		pep_seq[pm] = pep_seq_list
		pep_mass[pm] = pep_mass_list

		f = open(os.path.join("predicted_peptides_bound_top" + str(top_num), str(k_mer) + "-mer_" + str(pm) + "-mass_predicted.txt"), "w")
		f.write("peptides\taa_mass\n")
		for i in range(len(pep_seq[pm])):
			f.write(pep_seq[pm][i] + "\t" + pep_mass[pm][i] + "\n")
		f.close()

	bound_pep_seq = predict_hla_bound_peptides(pssm_path, top_num, pep_seq, pep_mass)
	'''
