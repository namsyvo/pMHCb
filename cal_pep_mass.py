import sys

aaMass = {}
f = open(sys.argv[1])
for line in f:
	tokens = line.strip().split("\t")
	aaMass[tokens[0]] = float(tokens[1])

pep = sys.argv[2]
m = 0
for p in pep:
	m += aaMass[p]
print "pep:", m
print "pep + H2O + H:", m + 19