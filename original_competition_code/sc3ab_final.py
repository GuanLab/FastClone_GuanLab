import numpy as np

def sc3ab_final():
	result_3b=open('3B.txt','wb')
	
	sc3b_matrix=open("3B_tmp.txt",'r')
	sc3b_matrix_T=open("3B_t.txt",'r')
#	a = (sc2b_matrix + sc3b_matrix + sc3b_matrix.T - 1).min()
	
	with open('2B.txt') as f:
		for line in f:
			line=line.rstrip()
			table= [float(x) for x in line.split('\t')]
			table= np.asarray(table)
		
			line1=sc3b_matrix.readline()
			line1=line1.rstrip()
			table1= [float(x) for x in line1.split('\t')]
			table1= np.asarray(table1)
			
			line2=sc3b_matrix_T.readline()
			line2=line2.rstrip()
			table2= [float(x) for x in line2.split('\t')]
			table2= np.asarray(table2)

			b=(table+table1+table2-1).max()
			if (b>0):
				tablenew=table1-b
			else:
				tablenew=table1

			tablenew= np.asmatrix(tablenew.clip(min=0.0))
			
			np.savetxt(result_3b, tablenew, fmt="%0.4f", delimiter="\t")
			
	return()
