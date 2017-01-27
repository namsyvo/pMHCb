'''
Generating peptide sequences from mass specs data using dynamic programming
Input: the right most pick (m), which is likely the mass of whole peptide sequences,
20 amoino acids with thier theoretical mass (isotopics)
Output: all possible peptide sequences of which theoretical mass sum equal to m
@Nam Sy Vo 2016
'''

#paras: float expMass, float expPPM, string[] aaName, float[] aaMass
def generate_sequences(m, v, aaName):

	aaSeqNumb = [[] for _ in range(len(v) + 1)]
	aaSeqName = [[] for _ in range(len(v) + 1)]
	aaSeqMass = [[] for _ in range(len(v) + 1)]

	for i in range(0, len(v) + 1):
		aaSeqNumb[i] = [-1 for _ in range(m + 1)]
		aaSeqName[i] = [[] for _ in range(m + 1)]
		aaSeqMass[i] = [[] for _ in range(m + 1)]

	#if m=0 then return empty set to get the m (1 set)
	for i in range(0, len(v) + 1):
		aaSeqNumb[i][0] = 1
		aaSeqName[i][0].append("")
		aaSeqMass[i][0].append("")

	#if no aa given, 0 ways to get the m
	for j in range(1, m + 1):
		aaSeqNumb[0][j] = 0

	for i in range(1, len(v) + 1):
		for j in range(1, m + 1):
			if v[i - 1] <= j:
				aaSeqNumb[i][j] = aaSeqNumb[i - 1][j] + aaSeqNumb[i][j - v[i - 1]]
				if aaSeqNumb[i - 1][j] != 0 and aaSeqNumb[i][j - v[i - 1]] != 0:
					#get Names
					aa = []
					for n in aaSeqName[i][j - v[i-1]]:
						aa.append(n + aaName[i-1])
					aaSeqName[i][j] = aaSeqName[i - 1][j] + aa
					#get Mass
					mass = []
					for p in aaSeqMass[i][j - v[i-1]]:
						mass.append(p + str(v[i-1]) + ";")
					aaSeqMass[i][j] = aaSeqMass[i - 1][j] + mass;
				elif aaSeqNumb[i][j - v[i - 1]] != 0:
					#get Names
					aa = []
					for n in aaSeqName[i][j - v[i-1]]:
						aa.append(n + aaName[i-1])
					aaSeqName[i][j] = aa
					#get Mass
					mass = []
					for p in aaSeqMass[i][j - v[i-1]]:
						if j - v[i - 1] == 0:
							mass.append(str(v[i-1]) + ";")
						else:
							mass.append(p + str(v[i-1]) + ";")
					aaSeqMass[i][j] = mass;
				elif aaSeqNumb[i - 1][j] != 0:
					aaSeqName[i][j] = aaSeqName[i - 1][j]
					aaSeqMass[i][j] = aaSeqMass[i - 1][j]
			else:
				#just copy the value from the top
				aaSeqNumb[i][j] = aaSeqNumb[i - 1][j]
				if aaSeqNumb[i - 1][j] != 0:
					aaSeqName[i][j] = aaSeqName[i - 1][j]
					aaSeqMass[i][j] = aaSeqMass[i - 1][j]
	'''
	#print out the matrix
	print
	print "# of seqs for each mass:"
	print "\t\t\tMass"
	print "\t\t\t",
	for j in range(m + 1):
		print str(j) + "\t",
	print
	print "\t\t",
	for j in range(m + 1):
		print "-\t",
	print

	print "Pepties\t",
	print("0  |\t"),
	for j in range(m + 1):
		print str(aaSeqNumb[0][j]) + "\t",
	print
	for i in range(1, len(v) + 1):
		print "\t\t" + str(v[i-1]) + "  |\t",
		for j in range(m + 1):
			print str(aaSeqNumb[i][j]) + "\t",
		print
	'''

	return aaSeqName[len(v)][m], aaSeqMass[len(v)][m]

def generate_peptides(expMass, expPPM, aaName, aaMass):

	M = expMass - 15.994915 - 2*1.007825 - 1.007825 # exclusive H2O and H+
	m = int(expMass) - 19; #substract mass of H2O and H+
	v = [int(x) for x in aaMass]
	aaMassDict = {}
	for i in range(len(aaName)):
		aaMassDict[aaName[i]] = aaMass[i]

	err_rate = expPPM * pow(10, -6)
	pep_seq, pep_mass = {}, {}
	for md in [-1, 0, 1]:
		aa_seq, aa_mass = generate_sequences(m + md, v, aaName)
		for i in range(len(aa_seq)):
			for k in [8, 9, 10, 11]:
				if len(aa_seq[i]) == k:
					mass = 0.0
					for aa in aa_seq[i]:
						mass += aaMassDict[aa]
					if abs((mass/M - 1)) <= err_rate:
						if k not in pep_seq:
							pep_seq[k] = []
							pep_mass[k] = []
						pep_seq[k].append(aa_seq[i])
						pep_mass[k].append(aa_mass[i])
	return pep_seq, pep_mass

#Testing
import sys

if __name__ == "__main__":

	expMass = float(sys.argv[1])
	expPPM = float(sys.argv[2]) 

	aaName, aaMass = [], []
	f = open(sys.argv[3])
	for line in f:
		tokens = line.strip().split("\t")
		aaName.append(tokens[0])
		aaMass.append(float(tokens[1]))

	pep_seq, pep_mass = generate_peptides(expMass, expPPM, aaName, aaMass)

	#print out the peptides
	print "Peptides for Mass " + str(expMass) + " with Delta Mass " + str(expPPM) + " PPM (aa mass " + str(int(expMass) - 19) + "):"
	for k in pep_seq:
		print "There are totally " + str(len(pep_seq[k])) + " candidate " + str(k) + "-mer sequences."

	for k in pep_seq:
		print str(k) + "-mer sequences:"
		for i in range(len(pep_seq[k])):
			print str(i + 1) + " : " + str(pep_seq[k][i]) + " (" + str(pep_mass[k][i]) + ")"