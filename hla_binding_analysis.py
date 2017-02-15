import sys
import os
import itertools

import PeptidePrediction

hla_binding = {}
seq_mass = {}
for (dirpath, dirnames, filenames) in os.walk(os.path.join(sys.argv[1])):
	for filename in filenames:
		if ".txt" not in filename:
			continue
		f = open(os.path.join(dirpath, filename))
		for line in f:
			tokens = line.strip().split()
			if tokens[2] not in hla_binding:
				hla_binding[tokens[2]] = set()
			if tokens[2] not in seq_mass:
				seq_mass[tokens[2]] = ""
			hla_binding[tokens[2]].add(tokens[1])
			seq_mass[tokens[2]] = filename
		f.close()

'''
hla_alleles = ["A0201", "A2902", "B0702", "B4403", "C0702"]
hla_num = {}
print "#peptides", len(hla_binding)
#for k, v in hla_binding.iteritems():
#	print k, seq_mass[k], len(v), "\t".join(v)
for k, v in hla_binding.iteritems():
	for hla in hla_alleles:
		for t in v:
			if hla in t:
				print k, t, seq_mass[k]
'''
aa_names, pssm_hla_1, pssm_hla_2 = PeptidePrediction.getPSSM(sys.argv[2])
pssm_hla = pssm_hla_1.copy()
pssm_hla.update(pssm_hla_2)
#print "# of HLA alleles:", len(pssm_hla)

check_hla = ["#HLA-A0201", "#HLA-A2902", "#HLA-B0702", "#HLA-B4403", "#HLA-C0702"]

top_num = 5
hla_path_score = {}
count = 0
for hla, pssm in pssm_hla.iteritems():
		top_score = []
		path_num = 1
		for pos_val in pssm:
			pssm_dict = {}
			for i in range(len(pos_val)):
				pssm_dict[aa_names[i]] = pos_val[i]
			top_val, top_name = [], []
			i = 0
			for k in sorted(pssm_dict, key=pssm_dict.get, reverse=True):
				i += 1
				if i > top_num: #take top_num AAs
					break
				top_val.append(pssm_dict[k])
				top_name.append(k)
			top_score.append(top_val)
			path_num *= (i-1)

		scores = []
		for t in itertools.product(*top_score):
			scores.append(sum(t))
		hla_path_score[hla] = sorted(scores, reverse=True)

		print "Finish", hla, hla_path_score[hla][0:10]

print len(hla_path_score)

for pep, hla_list in hla_binding.iteritems():
	'''
	score = []
	min_top5_score = []
	for hla in check_hla:
		s = 0
		for i in range(len(pep)):
			for j in range(len(aa_names)):
				if pep[i] == aa_names[j]:
					s += pssm_hla[hla][i][j]
					break
		score.append(s)
		min_s = 0
		for i in range(len(pep)):
			tmp = sorted(pssm_hla[hla][i], reverse=True)
			min_s += tmp[4]
		min_top5_score.append(min_s)
	print pep + "\t",
	for i in range(len(score)):
		print "\t" + str(score[i] - min_top5_score[i]),
	print
	'''
	'''
	for hla in hla_list:
		s = 0
		for i in range(len(pep)):
			for j in range(len(aa_names)):
				if pep[i] == aa_names[j]:
					s += pssm_hla[hla][i][j]
					break

		min_s = 0
		for i in range(len(pep)):
			tmp = sorted(pssm_hla[hla][i], reverse=True)
			min_s += tmp[4]

		max_s = 0
		for i in range(len(pep)):
			max_s += max(pssm_hla[hla][i])

		print pep + "\t" + hla + "\t" + seq_mass[pep] + "\t" + str(s) + "\t" + str(min_s) + "\t" + str(max_s)
	'''

	for hla in hla_list:
		print hla
		s = 0
		for i in range(len(pep)):
			for j in range(len(aa_names)):
				if pep[i] == aa_names[j]:
					s += pssm_hla[hla][i][j]
					break
		print s
		for score in hla_path_score[hla]:
			if score <= s:
				print pep + "\t" + hla + "\t" + seq_mass[pep] + "\t" + str(s) + "\t" + str(score)
				break