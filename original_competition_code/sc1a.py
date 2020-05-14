
import sc1b as sc1b

## the functions to derive 1A predictions from BATTENBURG only
def cal_1A_BAT_single(FILE_BAT,FILE_MUT): #BAT: input file;
	BAT=open(FILE_BAT,'r')
	line = BAT.readline()
	line=line.rstrip()
	title = line.split("\t")
	rho_total=0
	rho_count=0
	rho_total_2=0
	rho_count_2=0
	rho_total_back=0
	rho_count_back=0
	rho_total_length=0
	rho_count_length=0
	rho_count_alt=0
	for line in BAT:
		line=line.rstrip()
		table=line.split("\t");
		i=0
		dict={}
		for the_header in title:
			val=table[i]
			dict[the_header]=val;
			#print(the_header,val,dict[the_header])
			i=i+1
			pass
		dict['BAF']=float(dict['BAF'])
		if (dict['BAF']<0.5):
			dict['BAF']=1-dict['BAF']

		if (int(dict['nMaj1_A'])==int(dict['nMin1_A'])):
			if (dict['nMaj1_B'] =='NA'):
				if (int(dict['nMaj1_A'])==1):
					length_tmp=(int(dict['endpos'])-int(dict['startpos']));
					rho_count_alt+=1
			pass
#			print "equal"
		else:
			if (dict['nMaj1_B'] =='NA'):
				rho_tmp_up=1-2*float(dict['BAF'])
				rho_tmp_down=1-float(dict['nMaj1_A'])+(float(dict['nMin1_A'])+float(dict['nMaj1_A'])-2)*float(dict['BAF'])
				rho_tmp=rho_tmp_up/rho_tmp_down
				length_tmp=(int(dict['endpos'])-int(dict['startpos']));
				rho_total+=rho_tmp
				rho_count+=1
			
			#	print(dict['nMin1_A'],dict['nMaj1_A'],dict['BAF'],rho_tmp_up,rho_tmp_down,rho_tmp)
				pass
			else:
				rho_tmp_up=1-2*float(dict['BAF'])
				
				a_tmp=(float(dict['nMaj1_A'])+float(dict['nMin1_A']))*float(dict['frac1_A'])+(float(dict['nMaj2_A'])+float(dict['nMin2_A']))*float(dict['frac2_A'])-2
				b_tmp=float(dict['nMaj1_A'])*float(dict['frac1_A'])+float(dict['nMaj2_A'])*float(dict['frac2_A'])-1
				rho_tmp=rho_tmp_up/(a_tmp*float(dict['BAF'])-b_tmp)
				rho_total_2+=rho_tmp
				rho_count_2+=1
	#			print (dict['chr'],dict['endpos'])
			
	#			print (rho_tmp)

			pass
				
		pass
	if (rho_count_2>0):
		rho_final=rho_total_2/rho_count_2
	else:
		if ((rho_count)>3):
			rho_final=(rho_total+rho_total_2)/(rho_count+rho_count_2)
		else:
		## use cluster values
			(sc1b_number,cluster_freq,sc1c_cut,peak,fall)=sc1b.cal_1A_MUT_single(FILE_MUT,'battenburg_filled.txt',0.85)
			rho_final=cluster_freq[len(cluster_freq)-1]*2

	if (rho_count_alt>2):
		(sc1b_number,cluster_freq,sc1c_cut,peak,fall)=sc1b.cal_1A_MUT_single(FILE_MUT,'battenburg_filled.txt',0.85)
		#rho_final=cluster_freq[len(cluster_freq)-1]*2
		rho_final_tmp=cluster_freq[len(cluster_freq)-1]*2
		if ((rho_count_2>0) or (rho_count>3)):
			diff=abs(rho_final_tmp-rho_final)
			if (diff<0.08):
				if (rho_final_tmp==0.85):
					pass
				else:
					if (fall==1):
						pass
					else:			
						rho_final=rho_final_tmp
		else:
			rho_final=rho_final_tmp
				

	print('rho_count,rho_count_alt')
	print(rho_count,rho_count_alt)
	rho_final_length=1
	if (rho_final>1):
		rho_final=1
	if (rho_final<0):
		rho_final=0
	if (rho_final_length>1):
		rho_final_length=1
	if (rho_final_length<0):
		rho_final_length=0
	
	return (rho_final,rho_final_length)
	pass
