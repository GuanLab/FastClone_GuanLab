## yuanfang test
import het_utils as het
import re
import numpy as np
import b_fall as b_fall
from scipy.stats.kde import gaussian_kde
from scipy import signal
import copy
import random


## remove false snps and find number of peaks from the the single cnv regions
def cal_1A_MUT_single(FILE_MUT,FILE_BAT,rho): #MUT: input mutect file, BAT: input battenburg file
	(array)=het.record_single_cnv_region(FILE_BAT)
	pred_false=het.detect_false_snp_all(FILE_MUT)
	# test if pred_false length equals total length of the MUT file
	MUT=open(FILE_MUT,'r')
	line=MUT.readline()
	while (re.search('^##',line) is not None):
		line=MUT.readline()
		pass
	i=0
	for line in MUT:
		i+=1
	MUT.close()

	if (len(pred_false)!=i):
		die

	MUT=open(FILE_MUT,'r')
	line = MUT.readline()
	while (re.search('^##', line) is not None):
		line=MUT.readline()
		pass
	line=line.rstrip()
	title = line.split("\t")
	all_single_val=[]
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
				pass

			tumor=dict['tumor'].split(':')
			tumor_ratio=float(tumor[4])
			pos=int(dict['POS'])
			
			## only use 1:1 on regular chromosomes only
			if ((dict['#CHROM'] =='X') or (dict['#CHROM']=='Y')):
				tumor_ratio=tumor_ratio/2
				all_single_val.append(tumor_ratio)
			else:
				array_list=array[int(dict['#CHROM'])].split("\t")
				for member_list in array_list:
					if (member_list ==''):
						pass
					else:
						tmp=member_list.split('_')
						if ((pos>int(tmp[0])) and (pos<int(tmp[1])) and (int(tmp[2])==1) and (int(tmp[3])==1)):
							all_single_val.append(tumor_ratio)
							pass
						pass
		pred_i+=1

		
	length=len(all_single_val)

	if (length>100):
		fall=0
		pass
	else:
		all_single_val=b_fall.cal_1B_MUT_multi(FILE_MUT,FILE_BAT,rho) #MUT: input mutect file, BAT: input battenburg file
		fall=1
		pass


	nparam_density = gaussian_kde(all_single_val)
	ind = np.linspace(0,1,1000)
	kdepdf = nparam_density.evaluate(ind)
	cluster_num=0
	cluster=[]
	cut_all=[]
	peak=[]
	i_tmp=0
	old_val=kdepdf[i_tmp]
	up=0
	i_tmp+=1
	cut=0
	
	while (i_tmp<500):
		val=kdepdf[i_tmp]
		if (up==0):
			if (val>=old_val):
				old_val=val
			else:
				up=1
				cluster_num+=1
			#	print(i_tmp,cluster_num)
				freq=i_tmp/1000
				cluster.append(freq)
				cut_all.append(cut)
				peak.append(val)
				old_val=val
		else:
			if (up==1): ## down the hill
				if (val<=old_val):
					old_val=val
				else:
					up=0
					cut=i_tmp/1000
					old_val=val

		i_tmp+=1
	count=[]
	i_tmp=0;
	for i in cluster:
		count.append(0)
		i_tmp=i_tmp+1
	
	i_tmp_max=i_tmp
	cut_all_tmp=copy.copy(cut_all)
	cut_all_tmp.append(0.51)

	for val in all_single_val:
		i_tmp=0
		while (i_tmp<i_tmp_max):
			if ((val>cut_all_tmp[i_tmp]) and (val<cut_all_tmp[i_tmp+1])):
				count[i_tmp]=count[i_tmp]+1
			i_tmp=i_tmp+1

	cluster_num_final=0
	cluster_final=[]
	cut_all_final=[]
	peak_final=[]
	i_tmp=0
	for val in count:
		all_count=sum(count)
		r=val/all_count;
		if ((r>0.02) and (val>3)):
			cluster_num_final=cluster_num_final+1
			cluster_final.append(cluster[i_tmp])
			cut_all_final.append(cut_all[i_tmp])
			peak_final.append(peak[i_tmp])
		i_tmp=i_tmp+1
			
	if (len(cluster_final)<1):
		cluster_num_final=cluster_num
		cluster_final=copy.copy(cluster)
		cut_all_final=copy.copy(cut_all)
		peak_final=copy.copy(peak)

		
	if (rho==0.85):
		pass
	else:
		diff=(rho/2-cluster_final[len(cluster_final)-1])
	
		if (diff>0.05):
			cluster_num_final=cluster_num_final+1
		
		
		while ((diff<-0.05) and (len(cluster_final)>1)):
			cluster_num_final=cluster_num_final-1
			del cluster_final[-1]
			del cut_all_final[-1]
			del peak_final[-1]
			diff=(rho/2-cluster_final[len(cluster_final)-1])

	if (cluster_num_final<1):
		cluster_num_final=1

	print(cluster_final,cut_all_final,peak_final)
	return(cluster_num_final,cluster_final,cut_all_final,peak_final,fall)		
	
		
		

							
