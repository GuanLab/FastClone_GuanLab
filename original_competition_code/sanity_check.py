import sys,re

MUT=open(sys.argv[1],'r')
TRUTH=open(sys.argv[2],'r')
PRED=open(sys.argv[3],'r')

line=MUT.readline()
while (re.match('^##',line) is not None):
	line=MUT.readline()
	pass

for line in MUT:
	line=line.rstrip()	
	pred=PRED.readline()
	pred=pred.rstrip()
	if (re.search('True',line) is not None):
		line1=TRUTH.readline()
		line1=line1.rstrip()
		print (pred,line1,line)




