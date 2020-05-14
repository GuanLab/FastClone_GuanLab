import re
import os
import sys

def guess(FILE_MUT):
	result_file=open('1A.txt','w')
	result_file.write('0.85\n')
	result_file.close()

	result_file=open('1B.txt','w')
	result_file.write('2\n')
	result_file.close()

	
	result_file=open('2B.txt','w')
	result_file.close()

	
	result_file=open('3A.txt','w')
	result_file.close()

	result_file=open('3B.txt','w')
	result_file.close()

	
	
	MUT=open(FILE_MUT,'r')
	line=MUT.readline()
	while (re.search('^##',line) is not None):
		line=MUT.readline()
		pass
	
	i=0
	for line in MUT:
		i+=1
	MUT.close()

	result_file=open('1C.txt','w')
	false_n=int(i*0.05)
	cluster_1=int(i*0.5)
	cluster_2=i-false_n-cluster_1
	result_file.write("1\t%d\t0.85\n" % cluster_1)
	result_file.write("2\t%d\t0.42\n" % cluster_2)
	result_file.write("3\t%d\t0.00\n" % false_n)
	
	result_file.close()

	i_max=i
	j_max=i

	result_file=open('2B.txt','w')
	i=0
	while (i<i_max):
		j=0
		result_file.write('1')
		j+=1
		while (j<j_max):
			result_file.write('\t1')
			
			j+=1
		result_file.write('\n')
		
		i+=1
	result_file.close()
	

	result_file=open('2A.txt','w')
	i=0
	while (i<i_max):
		result_file.write('1\n')
		i+=1
	result_file.close()



	result_file=open('3A.txt','w')
	result_file.write('1\t0')
	result_file.close()
	result_file=open('3B.txt','w')
	i=0
	while (i<i_max):
		j=0
		result_file.write('0')
		j+=1
		while (j<j_max):
			result_file.write('\t0')
			
			j+=1
		result_file.write('\n')
		
		i+=1
	result_file.close()
	

#guess(sys.argv[2])
