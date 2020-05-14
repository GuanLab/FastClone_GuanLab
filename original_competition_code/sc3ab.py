### yuanfang mar-29-2015
import copy
import numpy as np

def cal_3_naive(cluster_freq,sc2a_result,all_result_by_cluster):
	
	freq_3a=copy.copy(cluster_freq)
	freq_3a.reverse()
	frequency_dict={}
	i=1
	for val in freq_3a:
		frequency_dict[i]=val*2
		i+=1

	tree=[]

	f=freq_3a[0]
	freq_3a.pop(0)
	f=f*2
	tree.append([])
	tree[0].append(1)
	tree[0].append(0)
	mydict={}
	mydict[1]=f
	mychild_sum={}
	mychild_sum[1]=0
	myfather={}
	
	myfather[1]='0'
	
	id=2

	while (len(freq_3a)>0):
		print(freq_3a)
		print(tree)
		print(id)
		f=freq_3a[0]
		freq_3a.pop(0)
		f=f*2
		if (len(mydict)==1):
			tree.append([])
			tree[1].append(2)
			tree[1].append(1)
			mychild_sum[1]+=frequency_dict[2]
			mydict[2]=f
			mychild_sum[2]=0
			myfather[2]='1'
			id+=1
		else:
			
			## tree try to grow from the largest node
			sorted_keys=sorted(mydict, reverse=True, key=mydict.get)
			print(sorted_keys)
			for root_tmp in sorted_keys:
				if ((f+mychild_sum[root_tmp])<frequency_dict[root_tmp]):
					tree.append([])
					tree[id-1].append(id)
					tree[id-1].append(root_tmp)
					mychild_sum[root_tmp]+=f
					mychild_sum[id]=0
					if (id in myfather):
						myfather[id]+='\t'+str(root_tmp)
						allfather=myfather[root_tmp].split('\t')
						print(allfather)
						for root_ttt in allfather:
							myfather[id]+='\t'+str(root_ttt)
			
					else:
						myfather[id]=str(root_tmp)
						allfather=myfather[root_tmp].split('\t')
						for root_ttt in allfather:
							myfather[id]+='\t'+str(root_ttt)
					break

				## values always sorted so no need to compare	
			mydict[id]=f

			id+=1
	
	tree_matrix=[]
	i=1
	while (i<id):
		tree_matrix.append([])
		j=1
		while (j<id):
			ch=0
			if i in myfather:
				
				allfather=myfather[i].split('\t')
				for root_ttt in allfather:
					if int(root_ttt)==j:
						ch=1
			tree_matrix[i-1].append(ch)
			j+=1
		i+=1
	
	all_result_by_cluster=np.array(all_result_by_cluster,dtype='f4')
	all_result_by_cluster=all_result_by_cluster[:,::-1]
	tree_matrix=np.array(tree_matrix)
#	print(all_result_by_cluster[:5,:])
	#sc2b_matrix=all_result_by_cluster.dot(all_result_by_cluster.T)
	sc3b_matrix=all_result_by_cluster.dot(tree_matrix.T).dot(all_result_by_cluster.T)
	np.fill_diagonal(sc3b_matrix, 0, wrap=False)
	np.savetxt('3B_tmp.txt',sc3b_matrix,'%.4f',delimiter='\t')
	sc3b_matrix=sc3b_matrix.T
	np.savetxt('3B_t.txt',sc3b_matrix,'%.4f',delimiter='\t')
	print('finished multiply')
	return tree
			
			
	
