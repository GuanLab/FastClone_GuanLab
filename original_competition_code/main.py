import sys,re
import het_utils as het
import numpy as np
import os

import sc1a as sc1a
import sc1b as sc1b
import sc1c as sc1c
import sc2ab as sc2ab
import sc3ab as sc3ab
import sc3ab_final as sc3ab_final
import guess as guess 

## calculate clonal frequency through battenburg
## usage: python3 __main__.py #battenburg file, ## cna file;

guess.guess(sys.argv[2])
#os.system("python3 guess.py "+sys.argv[1]+" "+sys.argv[2])

het.fillbat(sys.argv[1])

(sc1a_het,sc1a_het_length)=sc1a.cal_1A_BAT_single(sys.argv[1],sys.argv[2])
result_file=open('1A.txt','w')
result_file.write(str(sc1a_het)+'\n')
rho=sc1a_het
#print(sc1a_het_length)
result_file.close()


(sc1b_number,cluster_freq,sc1c_cut,peak,fall)=sc1b.cal_1A_MUT_single(sys.argv[2],'battenburg_filled.txt',rho)
result_file=open('1B.txt','w')
result_file.write('%d\n' % sc1b_number)
result_file.close()

print(cluster_freq,sc1c_cut)

(sc1c_result,frequency)=sc1c.cal_1C_MUT_top(sys.argv[2],'battenburg_filled.txt',cluster_freq,sc1c_cut,peak,rho)
print(sc1c_result,frequency)
result_file=open('1C.txt','w')
i=1
for number in sc1c_result:
	result_file.write('%d\t%d\t%.4f\n' % (i,number,frequency[i-1]))
	i=i+1
result_file.close()
	


#result_file.write(str(sc1b_number)+'\n')
#result_file.close()

(sc2a_result,result_length,label_matrix)=sc2ab.cal_2A_MUT_top(sys.argv[2],'battenburg_filled.txt',cluster_freq,sc1c_cut,peak,rho)
#for i in sc2a_result:
#	print(i)
#print(matrix)
#print(result_length)
result_file=open('2A.txt','w')
for number in sc2a_result:
	result_file.write('%d\n' % number)
result_file.close()

#(sc3a_result,sc3b_matrix)=sc3ab.cal_3_naive(cluster_freq,sc2a_result)
(tree)=sc3ab.cal_3_naive(cluster_freq,sc2a_result,label_matrix)
sc3ab_final.sc3ab_final()


#print(sc3a_result)
result_file=open('3A.txt','w')
l1=len(tree)
l2=len(tree[0])
i=0
while (i<l1):
	j=0
	result_file.write('%d' % tree[i][j])
	j+=1
	while (j<l2):
		result_file.write('\t%d'% tree[i][j])
		j+=1
	result_file.write('\n')
	i+=1


result_file.close()

os.system('rm 3B_tmp.txt')
os.system('rm 3B_t.txt')
