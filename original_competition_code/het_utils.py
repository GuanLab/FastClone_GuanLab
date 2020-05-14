## yuanfang test
import re, glob, os
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import b_fall as b_fall
import statistics

def record_single_cnv_region(FILE_BAT): #BAT: input file;
    BAT=open(FILE_BAT,'r')
    line = BAT.readline()
    line=line.rstrip()
    title = line.split("\t")
    #array_maj=np.zeros((24,300000000), dtype=np.int)
    #array_min=np.zeros((24,300000000), dtype=np.int)
    array = ['' for i in range(24)]
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
        if (dict['nMaj1_B'] =='NA'):
            chr=dict['chr']
            if (chr =='X'):
                chr=23
            else:
                chr=int(chr)
            i_start=dict['startpos']
            i_end=dict['endpos']
            #print(chr)
            if (array[chr]==''):
                array[chr]+=i_start+'_'+i_end+'_'+dict['nMin1_A']+'_'+dict['nMaj1_A']
            else:
                array[chr]+='\t'+i_start+'_'+i_end+'_'+dict['nMin1_A']+'_'+dict['nMaj1_A']

#        print(array[chr])
    return (array)




def record_allcorrect_cnv_region(FILE_BAT): #BAT: input file;
    BAT=open(FILE_BAT,'r')
    line = BAT.readline()
    line=line.rstrip()
    title = line.split("\t")
    #array_maj=np.zeros((24,300000000), dtype=np.int)
    #array_min=np.zeros((24,300000000), dtype=np.int)
    array = ['' for i in range(25)]
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
        chr=dict['chr']
        if (chr =='X'):
            chr=23
        elif (chr=='Y'):
            chr=24
        else:
            chr=int(chr)
        i_start=dict['startpos']
        i_end=dict['endpos']
        #print(chr)


        if (dict['nMaj2_A'] == 'NA'):
            attached=i_start+'_'+i_end+'_'+dict['nMin1_A']+'_'+dict['nMaj1_A']+'_'+dict['frac1_A']+'_'+dict['nMaj2_A']+'_'+dict['nMin2_A']+'_'+dict['frac2_A']
            if (array[chr]==''):
                array[chr]=attached
            else:
                array[chr]+="\t"+attached
        else:
            label=''

            if (dict['nMaj1_F']=='NA'):
                pass
            else:
                if ((int(dict['nMaj1_F'])==1)and(int(dict['nMin1_F'])==1)):
                    label='F'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
                if ((int(dict['nMaj2_F'])==1)and(int(dict['nMin2_F'])==1)):
                    label='F'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached

            if (dict['nMaj1_E']=='NA'):
                pass
            else:
                if ((int(dict['nMaj1_E'])==1)and(int(dict['nMin1_E'])==1)):
                    label='E'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
                if ((int(dict['nMaj2_E'])==1)and(int(dict['nMin2_E'])==1)):
                    label='E'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached

            if (dict['nMaj1_D']=='NA'):
                pass
            else:
                if ((int(dict['nMaj1_D'])==1)and(int(dict['nMin1_D'])==1)):
                    label='D'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
                if ((int(dict['nMaj2_D'])==1)and(int(dict['nMin2_D'])==1)):
                    label='D'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached

            if (dict['nMaj1_C']=='NA'):
                pass
            else:
                if ((int(dict['nMaj1_C'])==1)and(int(dict['nMin1_C'])==1)):
                    label='C'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
                if ((int(dict['nMaj2_C'])==1)and(int(dict['nMin2_C'])==1)):
                    label='C'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
            if (dict['nMaj1_B']=='NA'):
                pass
            else:
                if ((int(dict['nMaj1_B'])==1)and(int(dict['nMin1_B'])==1)):
                    label='B'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
                if ((int(dict['nMaj2_B'])==1)and(int(dict['nMin2_B'])==1)):
                    label='B'
                    attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                    if (array[chr]==''):
                        array[chr]=attached
                    else:
                        array[chr]+="\t"+attached
            if (label == ''):
                label='A'
                attached=i_start+'_'+i_end+'_'+dict['nMin1_'+label]+'_'+dict['nMaj1_'+label]+'_'+dict['frac1_'+label]+'_'+dict['nMaj2_'+label]+'_'+dict['nMin2_'+label]+'_'+dict['frac2_'+label]
                if (array[chr]==''):
                    array[chr]=attached
                else:
                    array[chr]+="\t"+attached



    #    print(array[chr])
    return (array)


def detect_false_snp_all(FILE_MUT):
    X_train=[]
    Y_train=[]
    X_test=[]
    false={}

    path = os.path.dirname(os.path.realpath(__file__)) + "/tmp_train_old"
    files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".truth.scoring_vcf.vcf")]
    for TRAIN in files:



        MUT=open(TRAIN,'r')
        line=MUT.readline()
        while (re.search('^##',line) is not None):
            line=MUT.readline()
            pass
        line=line.rstrip()
        title = line.split("\t")
        for line in MUT:
            line=line.rstrip()
            table=line.split("\t")
            i=0;
            dict={}
            for the_header in title:
                val=table[i]
                dict[the_header]=val
                i=i+1
                pass

            label=table[i];
            if (label=='False'):
                location=dict['#CHROM']+'_'+dict['POS']
                false[location]=0

    path = os.path.dirname(os.path.realpath(__file__)) + "/tmp_train_new"
    files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".truth.scoring_vcf.vcf")]
    for TRAIN in files:
        print(TRAIN)
        MUT=open(TRAIN,'r')
        line=MUT.readline()
        while (re.search('^##',line) is not None):
            line=MUT.readline()
            pass
        line=line.rstrip()
        title = line.split("\t")
        for line in MUT:
            line=line.rstrip()
            table=line.split("\t")
            i=0;
            dict={}
            for the_header in title:
                val=table[i]
                dict[the_header]=val
                i=i+1
                pass

            label=table[i];
            if (label=='False'):
                location=dict['#CHROM']+'_'+dict['POS']
                false[location]=0



    path = os.path.dirname(os.path.realpath(__file__)) + "/tmp_train"
    files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".truth.scoring_vcf.vcf")]
    for TRAIN in files:
        MUT=open(TRAIN,'r')
        line=MUT.readline()
        while (re.search('^##',line) is not None):
            line=MUT.readline()
            pass
        line=line.rstrip()
        title = line.split("\t")
        for line in MUT:
            line=line.rstrip()
            table=line.split("\t")
            i=0;
            dict={}
            for the_header in title:
                val=table[i]
                dict[the_header]=val
                i=i+1
                pass

            label=table[i];
            if (label=='False'):
                location=dict['#CHROM']+'_'+dict['POS']
                false[location]=0

    path = os.path.dirname(os.path.realpath(__file__)) + "/tmp_train_new"
    files = [os.path.join(path, d, f)
             for d in os.listdir(path)
             for f in os.listdir(os.path.join(path, d))
             if f.endswith(".truth.scoring_vcf.vcf")]
    for TRAIN in files:
        BAT=TRAIN
        BAT=re.sub('truth.scoring_vcf.vcf', 'battenberg.txt', BAT)
        ## need to return and detect false snps
#        (array)=record_top_cnv_region(BAT)
        MUT=open(TRAIN,'r')
        line=MUT.readline()
        while (re.search('^##',line) is not None):
            line=MUT.readline()
            pass
        line=line.rstrip()
        title = line.split("\t")
        for line in MUT:
            line=line.rstrip()
            table=line.split("\t")
            i=0;
            dict={}
            for the_header in title:
                val=table[i]
                dict[the_header]=val
                i=i+1
                pass

            label=table[i];
            if (label=='False'):
                location=dict['#CHROM']+'_'+dict['POS']
                false[location]=0
                label=1
            else:
                label=0
            Y_train.append(label)
            if (re.search('rs', dict['ID']) is not None):
                rs=1
            else:
                rs=0
                pass

            tumor=dict['tumor'].split(':')
            tumor_ratio=float(tumor[4])

            normal=dict['normal'].split(':')
            normal_ratio=float(normal[4])

            if ((dict['#CHROM']=='X') or (dict['#CHROM'] =='Y')):
                tumor_ratio=tumor_ratio/2
                normal_ratio=normal_ratio/2

            X_train.append([rs, tumor_ratio, normal_ratio])


    MUT=open(FILE_MUT,'r')
    line=MUT.readline()
    while (re.search('^##',line) is not None):
        line=MUT.readline()
        pass
    line=line.rstrip()
    title = line.split("\t")
    answer=[]
    j=0
    for line in MUT:
        line=line.rstrip()
        table=line.split("\t")
        i=0;
        dict={}
        for the_header in title:
            val=table[i]
            dict[the_header]=val
            i=i+1
            pass

        if (re.search('rs', dict['ID']) is not None):
            rs=1
        else:
            rs=0
            pass
        location=dict['#CHROM']+'_'+dict['POS']
        if location in false:
            answer.append(1)
        else:
            answer.append(0)
        j+=1

        tumor=dict['tumor'].split(':')
        tumor_ratio=float(tumor[4])

        normal=dict['normal'].split(':')
        normal_ratio=float(normal[4])
        X_test.append([rs, tumor_ratio, normal_ratio])
    forest = RandomForestRegressor(n_estimators = 300)
    forest = forest.fit(X_train,Y_train)
    pred=forest.predict(X_test)

    j=0
    for ppp in pred:
        if (answer[j] ==1):
            pred[j]=1
        j+=1

    return(pred)


def sd_single_largestc(FILE_MUT,FILE_BAT,cluster_freq,cut,rho):
    (array_single)=record_single_cnv_region(FILE_BAT)
#    print (array_single)
    pred_false=detect_false_snp_all(FILE_MUT)

    all_single_val=[]
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
        pred_i+=1
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
                    ## only use 1:1 on regular chromosomes only
            if ((dict['#CHROM'] =='X') or (dict['#CHROM']=='Y')):
                ratio1=tumor_ratio/2
                i=n_total-1
                if (tumor_ratio>cut[i]):
                    all_single_val.append(ratio1)
                pass
            else:
                array_list=array_single[int(dict['#CHROM'])].split("\t")
                for member_list in array_list:
                    if (member_list ==''):
                        pass
                    else:
                        tmp=member_list.split('_')
                        if ((pos>int(tmp[0])) and (pos<int(tmp[1])) and (int(tmp[2])==1) and (int(tmp[3])==1)):

                            i=n_total-1
                            if (tumor_ratio>cut[i]):
                                all_single_val.append(tumor_ratio)
                            pass
                        pass
                    pass
                pass
            pass
        pass

    length=len(all_single_val)
    if (length>100):
        pass
    else:
        all_single_val=b_fall.cal_1B_MUT_multi(FILE_MUT,FILE_BAT,rho) #MUT: input mu

    sd=statistics.stdev(all_single_val)
    MUT.close()
    return(sd)




def fillbat(FILE_BAT):

    length={}
    length['1']=249250621
    length['2']=243199373
    length['3']=198022430
    length['4']=191154276
    length['5']=180915260
    length['6']=180915260
    length['7']=159138663
    length['8']=146364022
    length['9']=141213431
    length['10']=135534747
    length['11']=135006516
    length['12']=133851895
    length['13']=115169878
    length['14']=107349540
    length['15']=102531392
    length['16']=90354753
    length['17']=81195210
    length['18']=78077248
    length['19']=59128983
    length['20']=63025520
    length['21']=48129895
    length['22']=51304566
    length['X']=155270560
    length['Y']=59373566


    BAT=open(FILE_BAT,'r')

    line=BAT.readline()
    line=line.rstrip()

    mydict={}
    for line in BAT:
        line=line.rstrip()
        table=line.split('\t')
        chrom=table[0]
        start=table[1]
        end=table[2]
        attach=start+'_'+end
        if chrom in mydict:
            mydict[chrom]+='\t'+attach
        else:
            mydict[chrom]=attach
    BAT.close()

    BAT=open(FILE_BAT,'r')
    NEW=open('battenburg_filled.txt','w')
    line=BAT.readline()
    line=line.rstrip()
    NEW.write(line+'\n')

    count={}
    count_end={}
    for line in BAT:
        line=line.rstrip()
        table=line.split('\t')

        all_pairs=mydict[table[0]].split("\t")
        number=len(all_pairs)
        if (number==1):
            NEW.write(table[0])
            NEW.write('\t1\t%d' % length[table[0]])
            i=3
            while (i<len(table)):
                NEW.write('\t'+table[i])
                i=i+1

        else:
            NEW.write(table[0])
            ## determine start point
            if table[0] in count:
                tmp1=all_pairs[count[table[0]]-1].split("_")
                tmp2=all_pairs[count[table[0]]].split("_")
                val=(int(int(tmp1[1])+int(tmp2[0]))/2)
                NEW.write('\t%d' % val)
                count[table[0]]+=1
            else:
                NEW.write('\t1')
                count[table[0]]=1
            ## determin the end point
            if table[0] in count_end:
                if (count_end[table[0]]<(number-1)):
                    tmp1=all_pairs[count_end[table[0]]].split("_")
                    tmp2=all_pairs[count_end[table[0]]+1].split("_")
                    val=(int(int(tmp1[1])+int(tmp2[0]))/2)
                    NEW.write('\t%d' % val)
                    count_end[table[0]]+=1
                else:
                    NEW.write('\t%d' % length[table[0]])
            else:
                tmp1=all_pairs[0].split("_")

                tmp2=all_pairs[1].split("_")
                val=(int(int(tmp1[1])+int(tmp2[0]))/2)
                NEW.write('\t%d' % val)
                count_end[table[0]]=1
            i=3
            while (i<len(table)):
                NEW.write('\t')
                NEW.write(table[i])
                i+=1
        NEW.write('\n')
    NEW.close()
    BAT.close()
