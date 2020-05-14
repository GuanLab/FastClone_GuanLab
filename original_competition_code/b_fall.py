## yuanfang
import het_utils as het
import re
import copy
import scipy
import scipy.stats as stats
import numpy as np


### calculate 1c

def cal_1B_MUT_multi(FILE_MUT,FILE_BAT,rho):
	
	(array)=het.record_allcorrect_cnv_region(FILE_BAT)
	pred_false=het.detect_false_snp_all(FILE_MUT)

	output=[]

	MUT=open(FILE_MUT,'r')
	line = MUT.readline()
	while (re.search('^##', line) is not None):
		line=MUT.readline()
		pass
	line=line.rstrip()
	title = line.split("\t")
	pred_i=0
	for line in MUT:
		pred_tmp=pred_false[pred_i]
		if (pred_tmp<=0.5):
			line=line.rstrip()
			table=line.split("\t")
			i=0
			dict={}
			for the_header in title:
				val=table[i]
				dict[the_header]=val
				i=i+1
			tumor=dict['tumor'].split(':')
			tumor_ratio=float(tumor[4])
			pos=int(dict['POS'])

			if ((dict['#CHROM']=='X') or(dict['#CHROM']=='Y')):
				# suppose no CNV on X
				ratio1=tumor_ratio/2
				ratio2=tumor_ratio/2
				ratio3=tumor_ratio/2
				output.append(ratio1)
			else:
				array_list=array[int(dict['#CHROM'])].split("\t")
				max_entropy=0;
				for member_list in array_list:

					tmp=[]
					if (member_list == ''):
						tmp.append(0)
						tmp.append(999999999)
						tmp.append(1)
						tmp.append(1)
						tmp.append(1)
						tmp.append('NA')
					else:
						tmp=member_list.split('_')
					if ((pos>int(tmp[0])) and (pos<int(tmp[1]))):
						## status when there is only one state:
						if (tmp[5] == 'NA'):
							## cnv is before the snp, this condition can either be the biggest clone or smaller clone
							up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
							ratio1=tumor_ratio*up/2	
							if (ratio1<0.6):
								print(ratio1)
								output.append(ratio1)
						
							## snp is before the cnv, this condition can only be the biggest clone 
							#(ratio2 and ratio3 comparing to parent clone only
							if (int(tmp[2])!=0):
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio2=tumor_ratio*up/(int(tmp[2]))/2
							else:
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio2=tumor_ratio*up/(int(tmp[3]))/2
							if (ratio2<0.6):
								output.append(ratio2)

							if (int(tmp[3])!=0):
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio3=tumor_ratio*up/(int(tmp[3]))/2
							else:
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio3=tumor_ratio*up/(int(tmp[2]))/2
							if (ratio3<0.6):
								output.append(ratio3)

							prior1=float(tmp[4]);

							# for cnv before snp compared to all possible clones
						### two state condition: 
						else:
							#three condistions: 1) cluster frequency > cnv proportion*rho; 2) cluster frequency <cnv proportion*rho; 3) (cnv frequency*rho+cluster frequency)<rho (double branch situation)
							tmp[2]=int(tmp[2])
							tmp[3]=int(tmp[3])
							tmp[5]=int(tmp[5])
							tmp[6]=int(tmp[6])
							if ((tmp[2]==1)and(tmp[3]==1)):
								prior=float(tmp[7])
								if (tmp[5]>tmp[6]):
									major=tmp[5]
									minor=tmp[6]
								else:
									major=tmp[6]
									minor=tmp[5]

							elif ((tmp[5]==1)and(tmp[6]==1)):
								prior=float(tmp[4])
								if (tmp[2]>tmp[3]):
									major=tmp[2]
									minor=tmp[3]
								else:
									major=tmp[3]
									minor=tmp[2]
						
							elif ((tmp[2]+tmp[3])==2):
								prior=float(tmp[7])
								if (tmp[5]>tmp[6]):
									major=tmp[5]
									minor=tmp[6]
								else:
									major=tmp[6]
									minor=tmp[5]
								
							elif ((tmp[5]+tmp[6])==2):
								prior=float(tmp[4])
								if (tmp[2]>tmp[3]):
									major=tmp[2]
									minor=tmp[3]
								else:
									major=tmp[3]
									minor=tmp[2]

							else:
								prior=float(tmp[7])
								if (tmp[5]>tmp[6]):
									major=tmp[5]
									minor=tmp[6]
								else:
									major=tmp[6]
									minor=tmp[5]
							

							cr=rho*prior
							
							major=int(major)
							minor=int(minor)
							up=2*(1-rho)+2*rho*(1-prior)+rho*prior*(major+minor);

							#1 only need to compare to the large frequencies
							# find cr with maj
							ratio1=(tumor_ratio*up-cr*major+cr)/2
							if (ratio1<0.6):
								output.append(ratio1)
							# find cr with min
							ratio2=(tumor_ratio*up-cr*minor+cr)/2

							if (ratio2<0.6):
								output.append(ratio2)
							## compare to the smaller frequencies
							ratio=(tumor_ratio*up)/2
							if (ratio<0.6):
								output.append(ratio)

	
		pred_i+=1
	
	return(output)

	
	


