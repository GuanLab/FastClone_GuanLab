## yuanfang
import het_utils as het
import re
import copy
import scipy
import scipy.stats as stats
import numpy as np


### calculate 1c

def cal_1C_MUT_top(FILE_MUT,FILE_BAT,cluster_freq,cut,peak,rho):
	sd=het.sd_single_largestc(FILE_MUT,FILE_BAT,cluster_freq,cut,rho)
	#print(sd)
	#print(cut, cluster_freq)
	#print (rho)
	(array)=het.record_allcorrect_cnv_region(FILE_BAT)
#	print(array)
	pred_false=het.detect_false_snp_all(FILE_MUT)
	result_file=open('bbb.txt','w')
	for number in pred_false:
		result_file.write('%f\n' % number)
	result_file.close()

	# test if pred_false length equals total lenth
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
	
	## assign each snp to a cluster
	n_total=len(cluster_freq)
	snp_count=[0] * (n_total+1) ## the last one being the false snp

	## find all false snps >0.5
	snp_count[n_total]=len([1 for i in pred_false if i > 0.5])	

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
				a=[ratio1,ratio2,ratio3]
				max_entropy=0
				for the_ratio in a:
					i=0
					while (i<n_total):
						value=float(stats.norm.cdf(the_ratio,cluster_freq[i],sd))
						if (value>0.5):
							value=1-value
						value=value*peak[i]

						if (value>max_entropy):
							max_entropy=value
							cluster_id_tmp=i
						i+=1

				#cluster_id_final=n_total-cluster_id_tmp
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
						
					if (n_total<3):
						if (tmp[5] == 'NA'):
							pass
						else:
							tmp[5]='NA'
							tmp[4]=1
							tmp[3]=1
							tmp[2]=1
					if ((pos>int(tmp[0])) and (pos<int(tmp[1]))):
						## status when there is only one state:
						if (tmp[5] == 'NA'):
							## cnv is before the snp, this condition can either be the biggest clone or smaller clone
							up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
							ratio1=tumor_ratio*up/2	
						
							## snp is before the cnv, this condition can only be the biggest clone 
							#(ratio2 and ratio3 comparing to parent clone only
							if (int(tmp[2])!=0):
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio2=tumor_ratio*up/(int(tmp[2]))/2
							else:
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio2=tumor_ratio*up/(int(tmp[3]))/2

							if (int(tmp[3])!=0):
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio3=tumor_ratio*up/(int(tmp[3]))/2
							else:
								up=2*(1-rho)+rho*(int(tmp[2])+int(tmp[3]));
								ratio3=tumor_ratio*up/(int(tmp[2]))/2
							prior1=float(tmp[4]);

							# for cnv before snp compared to all possible clones
							a=[ratio1]
							for the_ratio in a:
								i=0
								while (i<n_total):
									value=float(stats.norm.cdf(the_ratio,cluster_freq[i],sd))
									if (value>0.5):
										value=1-value

									value=value*peak[i]
									if (value>max_entropy):
										max_entropy=value
										cluster_id_tmp=i
									i+=1
							# for cnv after snp, compared to only the parent clone
							a=[ratio1,ratio2,ratio3]
							for the_ratio in a:
								i=n_total-1
								value=float(stats.norm.cdf(the_ratio,cluster_freq[i],sd))
								if (value>0.5):
									value=1-value
								value=value*peak[i]

								if (value>max_entropy):
									max_entropy=value
									cluster_id_tmp=i

				#			cluster_id_final=n_total-cluster_id_tmp

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
							# find cr with min
							ratio2=(tumor_ratio*up-cr*minor+cr)/2
							a=[ratio1,ratio2]
							for the_ratio in a:
								if (the_ratio>cr):
									i=0
									while (i<n_total):
										if (cluster_freq[i]>cr):
											value=float(stats.norm.cdf(the_ratio,cluster_freq[i],sd))
											if (value>0.5):
												value=1-value

											if (value>max_entropy):
												max_entropy=value
												cluster_id_tmp=i
										i+=1

		
							## compare to the smaller frequencies
							ratio=(tumor_ratio*up)/2
							the_ratio=ratio
							## plus the condition that cr and snv are on different branches
							if ((the_ratio<cr) or ((the_ratio+cr)<rho)):
								i=0
								while (i<n_total):
									if (cluster_freq[i]<cr):
										value=float(stats.norm.cdf(the_ratio,cluster_freq[i],sd))
										if (value>0.5):
											value=1-value
										value=value*peak[i]

										if (value>max_entropy):
											max_entropy=value
											cluster_id_tmp=i
									i+=1

			snp_count[cluster_id_tmp]+=1
				
	
		pred_i+=1
	
	last=snp_count.pop()
	snp_count.reverse()
	snp_count_final=snp_count
#	print(snp_count)
#	print(snp_count_final)
	snp_count_final.append(last)

	cluster_new=copy.copy(cluster_freq)
	cluster_new.reverse()

	cluster_freq_final=cluster_new
	if (last>0):
		cluster_freq_final.append(0.000)
	else:
		del snp_count_final[-1]
		
	
	i_tmp=0
	for val in cluster_freq_final:
		val=val*2
		cluster_freq_final[i_tmp]=val
		i_tmp+=1
	
	return(snp_count_final,cluster_freq_final)
