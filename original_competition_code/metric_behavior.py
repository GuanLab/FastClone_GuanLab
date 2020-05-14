from SMCScoring import *
import numpy as np
import csv

tsv_dir = './scoring_metric_data/text_files/' # directory to save tsv's to

def scoring1A_behavior(method='abs'):
    guesses = np.array(range(0,101))/100.0
    truths = np.array(range(10,110,10))/100.0
    res = []
    for i in range(len(truths)):
        for j in range(len(guesses)):
            res.append( [truths[i],guesses[j],calculate1A(guesses[j],truths[i], method)] )
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1A_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()


def scoring1B_behavior(method='normalized'):
    guesses = np.array(range(1,11))
    truths = np.array(range(1,6))
    res = [] 
    for i in range(len(truths)):
        for j in range(len(guesses)):
            res.append( [truths[i],guesses[j],calculate1B(guesses[j],truths[i], method)] )
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1B_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

def scoring1C_behavior(method='abs'):
    # Baseline Truth
    t_phis = np.array([.85,.5,.3])
    t_nssms = np.array([200,200,200])
    t_entry = zip(t_nssms,t_phis)
    n_iter = 100

    # Zero-mean noise in phi
    res = []
    concentration = [2, 5, 10, 100, 1000]
    for c in concentration:
        for i in range(n_iter):
            phis = []
            for p in t_phis:
                phis.append(np.random.beta(p*c,(1-p)*c))
            res.append([c,phis,calculate1C(t_entry,zip(t_nssms,phis), method)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1C_phi_ZM_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    # Systematic over/underestimation of phi
    sys_errors = np.array([.01,.03,.05,0.075,.10])
    sys_errors = np.concatenate((sys_errors,-sys_errors))
    res = []
    for sys_error in sys_errors:
        res.append([sys_error, calculate1C(t_entry,zip(t_nssms,t_phis+sys_error), method)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1C_phi_sys_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    # zero mean error in both phi and the number of SSMs
    res = []
    concentration = [2, 5, 10, 100, 1000]
    for c in concentration:
        for i in range(n_iter):
            rd = np.random.dirichlet([c,c,c]) * sum(t_nssms)
            rd = map(round,rd)
            remainder = sum(t_nssms) - sum(rd)
            #Randomly assign remainder
            rd[np.random.randint(0,len(rd))] += remainder
            rd = map(int,rd)
            res.append([c,rd,calculate1C(t_entry,zip(rd,t_phis), method)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1C_nssm_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    concentration = [2, 5, 10, 100, 1000, None]
    res = [[''] + ['beta_err_conc_nssm=' + str(x) for x in concentration]] # matrix of metric scores - phi error is in the rows, # SSMs error in the columns
    for c_phi in concentration:
        res_row = []
        for c_nssm in concentration:
            temp = []
            for i in range(n_iter):
                if not c_nssm is None: # if error should be added to the number of SSMs assigned to each subclone do it
                    # Determine the proportion of SSMs to assign to each cluster
                    rd_nssm = np.random.dirichlet([c_nssm,c_nssm,c_nssm]) * sum(t_nssms)
                    rd_nssm = map(round,rd_nssm)
                    remainder = sum(t_nssms) - sum(rd_nssm)
                    # Randomly assign remainder
                    rd_nssm[np.random.randint(0,len(rd_nssm))] += remainder
                    rd_nssm = map(int,rd_nssm)
                else:
                    rd_nssm = t_nssms

                if not c_phi is None: # similar for the phis
                    phis = [np.random.beta(p*c_phi,(1-p)*c_phi) for p in t_phis]
                else:
                    phis  = t_phis


                temp.append(calculate1C(t_entry,zip(rd_nssm,phis), method))
            res_row.append(np.mean(temp))
        res.append(res_row)
    res = [map(str,x) for x in res]
    res = [res[0]] + [['beta_err_conc_phi=' + str(concentration[i])] + res[i+1] for i in range(len(concentration))]
    res = ['\t'.join(x) for x in res]

    f = open(tsv_dir + 'scoring1C_interaction_behavior_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()            



    res = []
    # Collapse first two clusters
    phis = [(.85+.5)/2.0, .3]
    nssms = [400,200]
    entry = zip(nssms,phis)
    res.append(["Collapse12", calculate1C(t_entry,entry, method)])
    
    # Collapse last two clusters
    phis = [.85, (.5+.3)/2.0]
    nssms = [200,400]
    entry = zip(nssms,phis)
    res.append(["Collapse23", calculate1C(t_entry,entry, method)])

    # Collapse all clusters
    phis = [.55]
    nssms=[600]
    entry = zip(nssms,phis)
    res.append(["Collapse123", calculate1C(t_entry,entry, method)])

    # Assume all SSMs are clonal
    phis = [.85]
    nssms=[600]
    entry = zip(nssms,phis)
    res.append(["All_Clonal", calculate1C(t_entry,entry, method)])

    # For splits, phis were calculated as +/- 0.05 from center.  
    # 0.05 was obtained empirically by calculating the mean of the top and bottom half of a binomial sample with depth = 50. e.g.:
    # s = np.random.binomial(50,.85,(1000))
    # np.mean(s[s<np.median(s)]) / 50.0
    # 0.90
    # np.mean(s[s>np.median(s)]) / 50.0
    # 0.80

    # Split cluster 1
    phis = [.9,.8, .5, .3]
    nssms = [100,100,200,200]
    entry = zip(nssms,phis)
    res.append(["Split1", calculate1C(t_entry,entry, method)])

    # Split cluster 2
    phis = [.85, .55,.45, .3]
    nssms = [200,100,100,200]
    entry = zip(nssms,phis)
    res.append(["Split2", calculate1C(t_entry,entry, method)])
    # Split cluster 3
    phis = [.85, .5,.35,.25]
    nssms = [200,200,100,100]
    entry = zip(nssms,phis)
    res.append(["Split3", calculate1C(t_entry,entry, method)])
    
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring1C_cases_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()        

def scoring2A_behavior(tst_big_mat=True, tst_rand_reassign=True, tst_closest_reassign=True, method='default', verbose=False):
    '''Test the scoring behaviour of different metrics for evaluating Sub-Challenge 2, under various conditions.

    :param tst_big_mat: boolean for whether to test all the mistake scenarios with a larger number of clusters
    :param tst_rand_reassign: boolean for whether to test reassigning a portion of mutations to a random new cluster
    :param tst_closest_reassign: boolean for whether to test reassigning a portion of the mutations to the nearest cluster
    :param method: scoring metric to use
    :param verbose: boolean for whether to print output on the status of the function
    '''

    # Test the scoring behavior when the predicted CCM comes from one of the presepeficied 'mistake scenarios'
    # i.e. mistakes in the co-clustering assignment that we expect people to make
    if verbose:
        print('Testing scoring for SC2 using the ' + method + ' scoring metric:')
        print('Scoring behavior for mistake scenarios with 3 clusters...')
    size_clusters = 200 # true size of each cluster
    n_clusters = 3 # true number of clusters
    big_extra_num = 33 # number of mutations to add in the BigExtra case
    # True CCM:
    t_ccm, t_clusters = get_ccm('Truth',size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)
    scenarios = ["SplitClusterBot", "MergeClusterBot", "OneCluster", "NCluster", "SmallExtra", "BigExtra"] # mistake scenarios to be tested

    # Cases:
    res = []
    for sc in scenarios:
        # calculate the CCM
        ccm = get_ccm(sc ,t_ccm=t_ccm, t_clusters=t_clusters, size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)
        # calculate the score for the given scenario
        res.append([sc ,calculate2(ccm, t_ccm, method=method)])

    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open(tsv_dir + 'scoring2A_cases_' + method + '.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    if tst_big_mat:
        # Same thing but with more groups (all of the same size)
        # True CCM:
        n_clusters = 10 # true number of clusters
        big_extra_num = 10
        if verbose:
            print('Scoring behavior for mistake scenarios with ' + str(n_clusters) + ' clusters...')

        t_ccm, t_clusters = get_ccm('Truth',size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)

        # Cases:
        res_more_cl = []

        for sc in scenarios:
            # calculate the CCM
            ccm = get_ccm(sc ,t_ccm=t_ccm, t_clusters=t_clusters, size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)
            # calculate the score for the given scenario
            res_more_cl.append([sc ,calculate2(ccm, t_ccm, method=method)])

        res_more_cl = [map(str,x) for x in res_more_cl]
        res_more_cl = ['\t'.join(x) for x in res_more_cl]
        f = open(tsv_dir + 'scoring2A_big_cases_' + method + '.tsv', 'w')
        f.write('\n'.join(res_more_cl))
        f.close()

    if tst_rand_reassign or tst_closest_reassign:
        # Random re-assignment to arbitrary cluster and/or closest cluster with p=0.01,.03,.05,.1,.15,.25,.5 x 100
        if verbose:
            print('Scoring behavior of random reassignments to different clusters...')
            print('    Types of reassignments to perform:')
        res = {}
        if tst_rand_reassign:
            print('       - to random cluster')
            res['rand'] = []
        if tst_closest_reassign:
            print('        - to closest cluster')
            res['closest'] = []

        # Testing parameters
        p_errors = [0.01,0.03,0.05,.1,.15,.25,.5] # probability that a mutation is reassigned
        n_iter = 5 # number of iterations to perform for each p_err value
        size_clusters = 200 # true size of each cluster
        n_clusters = 5 # true number of clusters

        t_ccm, t_clusters = get_ccm('Truth',size_clusters=size_clusters, n_clusters=n_clusters)

        for p_err in p_errors:
            if verbose:
                print('     Using a probability of reassignment of ' + str(p_err) + ':')
            for i in range(n_iter):
                if verbose:
                    print('          Iteration ' + str(i+1) + ' of ' + str(n_iter) + '...')

                clusters = {}
                if tst_rand_reassign:
                    clusters['rand'] = np.copy(t_clusters)
                if tst_closest_reassign:
                    clusters['closest'] = np.copy(t_clusters)

                for j in range(t_ccm.shape[0]):
                    if np.random.random() < p_err:
                        if tst_rand_reassign:
                            cluster = np.argmax(clusters['rand'][j,:])
                            if np.random.random() < 0.5:
                                cluster -= 1
                            else:
                                cluster += 1
                            cluster = cluster % n_clusters

                            clusters['rand'][j,:] = 0
                            clusters['rand'][j,cluster] = 1

                        if tst_closest_reassign:
                            cluster = np.argmax(clusters['closest'][j,:])
                            if cluster == (n_clusters - 1):
                                cluster = (n_clusters - 2)
                            elif cluster == 0:
                                cluster = 1
                            else:
                                if np.random.random() < 0.5:
                                    cluster = (cluster + 1)
                                else:
                                    cluster = (cluster - 1)
                            clusters['closest'][j,:] = 0
                            clusters['closest'][j,cluster] = 1
                if tst_rand_reassign:
                    ccm = np.dot(clusters['rand'],clusters['rand'].T)
                    res['rand'].append([p_err,calculate2(ccm,t_ccm, method=method)])
                if tst_closest_reassign:
                    ccm = np.dot(clusters['closest'],clusters['closest'].T)
                    res['closest'].append([p_err,calculate2(ccm,t_ccm, method=method)])

        # Output the results
        if tst_rand_reassign:
            res['rand'] = [map(str,x) for x in res['rand']]
            res['rand'] = ['\t'.join(x) for x in res['rand']]
            f = open(tsv_dir + 'scoring2A_random_reassignment_' + method + '.tsv', 'w')
            f.write('\n'.join(res['rand']))
            f.close()
        if tst_closest_reassign:
            res['closest'] = [map(str,x) for x in res['closest']]
            res['closest'] = ['\t'.join(x) for x in res['closest']]
            f = open(tsv_dir + 'scoring2A_closest_reassignment_' + method + '.tsv', 'w')
            f.write('\n'.join(res['closest']))
            f.close()

    
def scoring2B_behavior(tst_betas=True, tst_prob_mod=True, tst_prob_mod_err=True, method='pseudoV', verbose=True):
    if tst_betas:
        if verbose:
            print 'Testing adding beta error to the true CCM:'
        t_ccm, t_clusters = get_ccm('Truth',size_clusters=200, n_clusters=3, big_extra_num=33)

        n_uniq = len(np.triu_indices(t_ccm.shape[0],k=1)[0])
        res = []
        concentrations = [1000,100,50,25,10,5,3,1]
        n_iter = 5 # number of iterations to perform for each concentration value

        for c in concentrations:
            if verbose:
                print ('       Beta Error with concentration param = ' + str(c) + '...')
            for i in range(n_iter):
                if verbose:
                    print('           Iteration ' + str(i+1) + ' of ' + str(n_iter) + '...')

                ccm = np.copy(t_ccm)
                ccm[np.triu_indices(t_ccm.shape[0],k=1)] -= np.random.beta(1,c,n_uniq) # subtract beta error from the upper triangular part of the true CCM
                ccm[np.tril_indices(t_ccm.shape[0],k=-1)] = 0 # ensure the matrix is symmetrical
                ccm = ccm + ccm.T
                np.fill_diagonal(ccm,1) # ensure the matrix has 1's along the diagonal
                ccm = np.abs(ccm) # ensure the matrix has values between 0 and 1
                res.append([c,calculate2(ccm,t_ccm, method=method)])
        res = [map(str,x) for x in res]
        res = ['\t'.join(x) for x in res]
        f = open(tsv_dir + 'scoring2B_beta_' + method + '.tsv', 'w')
        f.write('\n'.join(res))
        f.close()

    if tst_prob_mod:
        size_clusters = 200
        n_clusters = 3
        big_extra_num = 33

        t_ccm, t_clusters = get_ccm('Truth', size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)

        TwoBscenarios = ['SmallExtra', "BigExtra", "OneCluster", "NCluster", "SplitClusterBot", "MergeClusterBot"]
        scoring_data = {}
        if tst_prob_mod_err:
            error_data = {}
        for sc in TwoBscenarios:
            if verbose:
                print('  Scenario: ' + sc)
            ccm = get_ccm(sc,t_ccm=t_ccm, t_clusters=t_clusters, size_clusters=size_clusters, n_clusters=n_clusters, big_extra_num=big_extra_num)
            ccm_ones = (ccm == 1)
            ccms = [ccm]
            if verbose:
                print '       Calculating SC2 score with certainty of 1...'
            data = [calculate2(ccm, t_ccm, method=method)]

            diag = np.diag_indices(ccm.shape[0])

            if tst_prob_mod_err:
                stds = [0.01,0.03,0.05,0.1,0.15,0.2]
                errs = {} # error matrix for each std value - keep error matrix the same for all probabalistic ccm's
                data_err = {} # ccm's for each certainty level with including small errors
                for std in stds:
                    err = np.random.normal(0,std,ccm.shape) # add some random error to each entry in the matrix
                    err = np.triu(err, 1) # make sure the error matrix is symmetrical and the diagonal is all 0's
                    err = err + np.transpose(err)
                    errs[std] = err

                for std in stds:
                    ccm_err = (ccm - errs[std])

                    ccm_err[ccm_err > 1] = 1 # make sure all values are between 0 and 1
                    ccm_err[ccm_err < 0] = 0

                    ccms_err = [ccm_err]
                    if verbose:
                        print '       Calculating SC2 score with certainty of 1 with Errors using std dev of ' + str(std) + '...'
                    data_err[std] = [calculate2(ccm_err, t_ccm, method=method)]

            probs = [0.95,0.9,0.85,0.8,0.75,0.7]
            for prob in probs: # prob is your certainty of your results
                ccm_prob = np.copy(ccm) + (1-prob) # change 0's to 1-prob
                ccm_prob[ccm_ones] -= 2*(1-prob) # change 1's to prob
                ccm_prob[diag] = 1 # diagonal should always be ones
                ccms.append(ccm_prob)

                if verbose:
                    print '       Calculating SC2 score with certainty of ' + str(prob) + '...'
                data.append(calculate2(ccm_prob, t_ccm, method=method))

                if tst_prob_mod_err:
                    for std in stds:
                        ccm_prob_err = (ccm_prob - errs[std])

                        ccm_prob_err[ccm_prob_err > 1] = 1 # make sure all values are between 0 and 1
                        ccm_prob_err[ccm_prob_err < 0] = 0

                        ccms_err.append(ccm_prob_err)
                        if verbose:
                            print '       Calculating SC2 score with certainty of ' + str(prob) + ' with Errors using std dev of ' + str(std) + '...'
                        data_err[std].append(calculate2(ccm_prob_err, t_ccm, method=method))

            scoring_data[sc] = data
            if tst_prob_mod_err:
                error_data[sc] = data_err

        with open(tsv_dir + '2B_prob_scoring_' + method + '.tsv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Scenario', 1]+probs)
            for key, val in scoring_data.iteritems():
                writer.writerow([key] + val)

        if tst_prob_mod_err:
            for std in stds:
                with open(tsv_dir + '2B_prob_scoring_with_err_' + str(std) + '_' + method + '.tsv', 'w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Scenario', 1]+probs)
                    for key, val in error_data.iteritems():
                        writer.writerow([key] + val[std])





def scoring3A_behavior(method="orig", verbose=False, weights=None, save=True, pc_amount='more', full_matrix=True, in_mat=2):
    '''Scoring behaviour of subchallenge 3 metrics

    Attributes:
    :param method: method or list of methods to use when evaluating each subchallenge scenario
    :param verbose: boolean for whether to output details of the scoring metrics
    :param weights: weights to pass into the scoring function for calculating the weighted average of multiple metrics
    :param save: boolean for whether or not to save the results of the scoring behaviour to a tsv file
    :param pc_amount: amount of pseudo counts to include for each matrix, 'more', 'less' or 'none'
    :param in_mat: number representing which matrices to use in calculating the SC3 scoring metric
        Options:
            1 - use all input matrics i.e. CCM, ADM, ADM^T and CM
            2 - use all except co-clustering matrix (CCM)
            3 - use all except ancestor descendant matrix (ADM)
            4 - use all except ADM^T
            5 - use all except cousin matrix (CM)
    '''
    # True values for each attribute
    t_ccm, t_clusters = get_ccm("Truth")
    t_ad = get_ad("Truth")

    res = list() # results of using the given method to score each scenario

    if pc_amount is 'more':
        n_pc = np.sqrt(t_ccm.shape[0]) # number of pseudo counts to include
        pc_ext = "_more_pc" # file name extension for the amount of pseudo counts to include
    elif pc_amount is 'less':
        n_pc = np.log(t_ccm.shape[0])
        pc_ext = ""
    elif pc_amount is 'none':
        n_pc = 0
        pc_ext = "_no_pc"

    in_mat_ext = {1:'_all',
                  2:'_nc',
                  3:'_na',
                  4:'_nat',
                  5:'_ncous'}[in_mat]

    for scenario in scenarios:
        if verbose:
            print '\nScenario %s' % scenario

        if scenario == 'Truth':
            res.append(['Truth',
                        calculate3(t_ccm,t_ad,t_ccm,t_ad,
                                   method=method, verbose=verbose, weights=weights, pseudo_counts=n_pc, full_matrix=full_matrix, in_mat=in_mat)])
        else:
            ccm = get_ccm(scenario, t_ccm=t_ccm)
            ad = get_ad(scenario, t_ad=t_ad)
            res.append([scenario,
                        calculate3(ccm,ad,t_ccm,t_ad,
                                   method=method, verbose=verbose, weights=weights, pseudo_counts=n_pc, full_matrix=full_matrix, in_mat=in_mat)])


    if save:
        if full_matrix:
            tri_ext = "_full"
        else:
            tri_ext = "_triu"
        if isinstance(method, list):
            f = open(tsv_dir + 'scoring3A_all_cases_' + '_'.join([m + in_mat_ext for m in method]) + pc_ext + tri_ext + '.tsv', 'w')
        else:
            f = open(tsv_dir + 'scoring3A_all_cases_' + method + in_mat_ext +  pc_ext + tri_ext + '.tsv', 'w')
        out_res = [map(str,x) for x in res]
        out_res = ['\t'.join(x) for x in out_res]
        f.write('\n'.join(out_res))
        f.close()

    return res

def scoring3A_behavior_all(verbose=True):
    for method in ['pseudoV',
                    'orig',
                    'mcc',
                    'pearson',
                    'spearman',
                    'aupr',
                    'sqrt',
                    'sym_pseudoV',
                   ['pseudoV', 'mcc', 'pearson'],
                   ['pseudoV', 'pearson', 'sym_pseudoV'],
                   ['aupr', 'sqrt', 'sym_pseudoV'],
                   ['aupr', 'sqrt', 'sym_pseudoV', 'pearson']]:
        for fm in [True, False]:
            for pc in ['none', 'less', 'more']:
                for input in range(5):
                    print 'Starting %s - Pseudo Counts: %s - Full Matrix: %s, Input Matrices Index: %s' % (method,pc,fm, input)

                    scoring3A_behavior(method=method, verbose=verbose,pc_amount=pc, full_matrix=fm, in_mat=input+1)
                    print 'Done %s - Pseudo Counts: %s - Full Matrix: %s, Input Matrices Index: %s' % (method,pc,fm, input)

def scoring3A_weight_behavior(methods=["pseudoV", "pearson", "sym_pseudoV"], verbose=False, res=None, in_mat=2):
    '''Create the data on how the weights used in subchallenge 3 affect the score using the given scoring methods

    Attributes:
    :param methods: list of methods to use when evaluating each subchallenge scenario
    :param verbose: boolean for whether to output details of the scoring metrics
    '''
    # True values for each attribute
    if res is None:
        res = np.transpose(np.asarray([[row[1] for row in scoring3A_behavior(method, verbose=verbose)] for method in methods]))
    print res

    wght_res = {'Case':scenarios}
    n_method = len(methods)
    weights = [0,0.5,1]
    weight_list = get_weights(n_method, weights)
    print weight_list
    for wght in weight_list:
        norm_wght = wght / float(sum(wght))
        scores = np.sum(norm_wght * res, 1)
        if set(wght) == {0,weights[-1]} or set(wght) == {weights[-1]}:
            key = '+'.join([methods[i] for i in range(n_method) if wght[i] == weights[-1]])
            wght_res[key] = scores
        else:
            wght_res[str(wght)] = scores


    in_mat_ext = {1:'_all',
                  2:'_nc',
                  3:'_na',
                  4:'_nat',
                  5:'_ncous'}[in_mat]

    with open(tsv_dir + 'weights3A_all_cases_' + '_'.join([m + in_mat_ext for m in methods]) + '.tsv', 'wb') as f:
        fields = sorted(wght_res.keys())
        print fields
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        for i in range(len(scenarios)):
            w.writerow({field:wght_res[field][i] for field in fields})

    return wght_res, res

def get_weights(n_objects, weights, cur_weight=None, list_weights= list()):
    """Generate a list of all possible unique combinations of weight assignments from the given list
    of weights. Does not make assumptions about the contents of the list, so some weight combinations
    may be duplicates if the list of weights contains elements that are multiples of each other.

    :param n_objects: number of objects that need to be assigned weights
    :param weights: list of possible weights
    :param cur_weight: current iteration of weight assignments
    :param list_weights: list of all weight assignments generated so far
    :return: list of unique weight assignments
    """
    if n_objects is 0:
        if len(set(cur_weight)) > 1 or cur_weight[0] is weights[-1]: # only include the case where all the weights are the same once
            if not (len(set(cur_weight)) <= 2 and 0 in set(cur_weight)) or weights[-1] in set(cur_weight): # only include the cases where one or more of the weights
                                                                                                        # are zero and all the other weights are the same once
                wght = np.copy(cur_weight)
                list_weights.append(wght)
    else:
        if cur_weight is None:
            cur_weight = [0] * n_objects
        for i in weights:
            cur_weight[n_objects-1] = i
            get_weights(n_objects-1, weights, cur_weight, list_weights)

    return np.asarray(list_weights)


def overall_score(p_cell, t_cell, p_ncluster, t_ncluster, p_1c, t_1c,  p_ccm, t_ccm, p_ad, t_ad, weights = [2000, 2000, 3000, 6000, 7000], verbose=False):
    '''Get the overall score for a prediction

    Attributes:
    :params p_cell, t_cell: predicted and true cellularity
    :params p_ncluster, t_ncluster: predicted and true number of clusters
    :params p_1c, t_1c: zippped variable with the size of each clusters and the cellullar frequency of each cluster
    :params p_ccm, t_ccm: co-clustering matrix for the mutations
    :params p_ad, t_ad: ancestry-descendant matrix for the mutations
    :param weights: weights for scores from each challenge when combining to give an overall score

    '''
    s1A = calculate1A(p_cell, t_cell)
    s1B = calculate1B(p_ncluster, t_ncluster)
    s1C = calculate1C(p_1c, t_1c)
    s2 = calculate2(p_ccm, t_ccm)
    s3 = calculate3(p_ccm, p_ad, t_ccm, t_ad)
    scores = [s1A,s1B,s1C,s2,s3]

    soverall = np.sum(np.multiply(scores, weights))

    scores.append(soverall)
    if verbose:
        print('Scores:\nSC1 - Part A: %s, Part B: %s, Part C: %s\nSC2 - %s\nSC3 - %s\nOverall: %s' % tuple(scores))

    return soverall

def score_all(l_p_cell, t_cell,
              l_p_ncluster, t_ncluster,
              l_p_1c, t_1c,
              l_p_ccm, t_ccm,
              l_p_ad, t_ad,
              weights = [2000, 2000, 3000, 6000, 7000],
              verbose = False):
    """Calculate the overall score across all sub-challenges for multiple scenarios.

    :param l_p_cell: list of predicted cellularity for each scenario
    :param t_cell: true cellularity
    :param l_p_ncluster: list of predicted number of mutation clusters for each scenario
    :param t_ncluster:
    :param l_p_1c: list of zip objects containing the predicted number of mutations in each cluster and the predicted
        cellular frequency of each cluster. one zip object for each scenario.
    :param t_1c: zip object containing the true number of mutations in each cluster and the true cellular frequency
        of each cluster
    :param l_p_ccm: list of predicted co-clustering matrices for each scenario
    :param t_ccm: true co-clustering matrix
    :param l_p_ad: list of predicted ancestor-descendant matrices for each scenario
    :param t_ad: true ancestor-descendant matrix
    :param weights: weights for each part of each sub-challenge tp be used when calculating the overall score
    :param verbose: boolean for whether to print information about the function execution
    :return: list of overall scores for each of the given scenarios
    """
    scores = dict()
    for i in range(len(l_p_cell)):
        if verbose:
            print 'Scenario: %s' % scenarios[i]
        score = overall_score(
            l_p_cell[i], t_cell,
              l_p_ncluster[i], t_ncluster,
              l_p_1c[i], t_1c,
              l_p_ccm[i], t_ccm,
              l_p_ad[i], t_ad,
              weights,
              verbose)
        scores[scenarios[i]] = score

    if verbose:
        print ("Overall Scores")
        print scores

    return scores

def rank(scores):
    '''Calculate the rank order of a list of scores (either from a single sub-challenge or overall scores)

    :param scores: list of scores
    :return: list of relative ranking of each score
    '''
    scores = np.asarray(scores)
    order = scores.argsort()
    rank = order.argsort()
    return rank + 1

def get_ccm(scenario, t_ccm=None, t_clusters=None, size_clusters=100, n_clusters=6, big_extra_num=15, nssms=None):
    '''Find the co-clustering matrix for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :params t_ccm, t_clusters: - optional value for the true co-clustering matrix and the true cluster assignments,
            to avoid computing it multiple times
    :param return: (if scenario is 'Truth') the true co-clustering matrix and the clusters used to generate this matrix
                    (otherwise) the co-clustering matrix fro the given scenario
    '''

    if t_clusters is None:
        t_clusters = np.zeros((n_clusters*size_clusters,n_clusters))
        for i in range(n_clusters):
            t_clusters[i*size_clusters:(i+1)*size_clusters,i] = 1 #assign each cluster
    if t_ccm is None and nssms is None:
        t_ccm = np.dot(t_clusters,t_clusters.T)

    if scenario in "Truth":
        return t_ccm, t_clusters
    elif "ParentIs" in scenario:
        return t_ccm
    elif scenario is "OneCluster":
        return np.ones((nssms,nssms))
    elif "NCluster" in scenario:
        return np.identity(nssms)
    elif "SplitCluster" in scenario:
        clusters = np.zeros((n_clusters*size_clusters,n_clusters+1))
        clusters[:,:-1] = np.copy(t_clusters)
        if "SplitClusterBot" in scenario:
            clusters[(n_clusters-0.5) * size_clusters:n_clusters*size_clusters,n_clusters-1] = 0
            clusters[(n_clusters-0.5)*size_clusters:n_clusters*size_clusters,n_clusters] = 1
            return np.dot(clusters, clusters.T)
            return np.dot(clusters, clusters.T)
        elif scenario is "SplitClusterMidOneChild":
            clusters[2.5*size_clusters:3*size_clusters,2] = 0
            clusters[2.5*size_clusters:3*size_clusters,n_clusters] = 1
            return np.dot(clusters, clusters.T)
        elif scenario is "SplitClusterMidMultiChild":
            clusters[1.5*size_clusters:2*size_clusters,1] = 0
            clusters[1.5*size_clusters:2*size_clusters,n_clusters] = 1
            return np.dot(clusters, clusters.T)

    elif scenario is "MergeClusterBot":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[(n_clusters-1)*size_clusters:n_clusters*size_clusters,n_clusters-2] = 1 #fix cluster 5 (originally cluster 6)
        clusters[(n_clusters-2)*size_clusters:(n_clusters-1)*size_clusters,n_clusters-2] = 0 #merge clusters 4 and 5 (from true phylogeny)
        clusters[(n_clusters-3)*size_clusters:(n_clusters-2)*size_clusters,n_clusters-3] = 1
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterMid&BotOneChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[5*size_clusters:6*size_clusters, 2] = 1 #merge clusters 3 and 6
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterMid&BotMultiChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[5*size_clusters:6*size_clusters, 4] = 1 #fix cluster 5 (originally cluster 6)
        clusters[4*size_clusters:5*size_clusters, 1] = 1 #merge clusters 2 and 5 (from true phylogeny)
        clusters[4*size_clusters:5*size_clusters, 4] = 0
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterTop&Mid":
        clusters = np.zeros((6*size_clusters,5))
        clusters[0:2*size_clusters, 0] = 1 #merged cluster from clusters 1 and 2 in true phylogeny
        clusters[2*size_clusters:3*size_clusters, 1] = 1
        clusters[3*size_clusters:4*size_clusters, 2] = 1
        clusters[4*size_clusters:5*size_clusters, 3] = 1
        clusters[5*size_clusters:6*size_clusters, 4] = 1
        return np.dot(clusters, clusters.T)
    if 'Extra' in scenario:
        clusters = np.zeros((n_clusters*size_clusters,n_clusters+1))
        clusters[:,:-1] = np.copy(t_clusters)
        if "SmallExtra" in scenario:
            num_extra = 1
        elif "BigExtra" in scenario:
            num_extra = big_extra_num

        for i in range(n_clusters):
            clusters[size_clusters*i:size_clusters*i+num_extra,i] = 0
            clusters[size_clusters*i:size_clusters*i+num_extra,n_clusters] = 1
        return np.dot(clusters,clusters.T)
    else:
        raise LookupError("Invalid scenario")

def get_ad(scenario, t_ad=None, size_clusters=100, nssms=None):
    '''Find the ancestry-descendant matrix for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param t_ad: optional value for the true AD matrix, to avoid comupting it multiple times
    '''
    if t_ad is None and nssms is None:
        t_ad = np.zeros((6*size_clusters,6*size_clusters))
        t_ad[0:size_clusters, size_clusters:] = 1
        t_ad[size_clusters:2*size_clusters, 3*size_clusters:5*size_clusters] = 1
        t_ad[2*size_clusters:3*size_clusters, 5*size_clusters:6*size_clusters] = 1

    if scenario in ["Truth", "SplitClusterBotSame", "MergeClusterBot"]:
        return t_ad
    elif scenario is "SplitClusterBotDiff":
        ad = np.copy(t_ad)
        ad[5*size_clusters:5.5*size_clusters,5.5*size_clusters:6*size_clusters] = 1.
        return ad
    elif scenario is "SplitClusterMidOneChild":
        ad = np.copy(t_ad)
        ad[2.5*size_clusters:3*size_clusters,5*size_clusters:6*size_clusters] = 0
        return ad
    elif scenario is "SplitClusterMidMultiChild":
        ad = np.copy(t_ad)
        ad[1.5*size_clusters:2*size_clusters,3*size_clusters:5*size_clusters] = 0
        return ad
    elif scenario == "MergeClusterMid&BotOneChild":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, 5*size_clusters:6*size_clusters] = 0
        return ad
    elif scenario == "MergeClusterMid&BotMultiChild":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters, 4*size_clusters:5*size_clusters] = 0
        ad[4*size_clusters:5*size_clusters, 3*size_clusters:4*size_clusters] = 1
        return ad
    elif scenario == "MergeClusterTop&Mid":
        ad = np.zeros((6*size_clusters,6*size_clusters))
        ad[0:2*size_clusters, 2*size_clusters:] = 1
        ad[2*size_clusters:3*size_clusters, 5*size_clusters:] = 1
        return ad
    elif scenario is "ParentIsSibling":
        ad = np.copy(t_ad)
        ad[3*size_clusters:4*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsGrandparent":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        return ad
    elif scenario is "ParentIsAunt":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsCousin":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        ad[5*size_clusters:6*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsSiblingWithChildren":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        return ad
    elif scenario is "ParentIsNieceWithChildren":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        ad[5*size_clusters:6*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 6's ancestry
        return ad
    elif scenario is "OneCluster":
        return np.zeros((nssms,nssms))
    elif scenario is "NClusterOneLineage":
        return np.triu(np.ones((nssms,nssms)), k=1)
    elif scenario is "NClusterTwoLineages":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[2:3*size_clusters+2,3*size_clusters+2:] = 0
        return ad
    elif scenario is "NClusterCorrectLineage":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[size_clusters:2*size_clusters,range(2*size_clusters,3*size_clusters)+range(5*size_clusters,6*size_clusters)] = 0 # equivalent of cluster 2 from true AD matrix
        ad[2*size_clusters:3*size_clusters,3*size_clusters:5*size_clusters] = 0 # cluster 3 from true AD matrix
        ad[3*size_clusters:4*size_clusters,4*size_clusters:] = 0 # cluster 4 from true AD matrix
        ad[4*size_clusters:5*size_clusters,5*size_clusters:6*size_clusters] = 0 # cluster 5 from true AD matrix
        return ad
    elif scenario is "SmallExtraNewBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters)+range(size_clusters+1,2*size_clusters)+range(3*size_clusters+1,4*size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraCurBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters)+range(2*size_clusters+1,3*size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraMid":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraTop":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[size_clusters*i,
               range(1,size_clusters)+
               range(size_clusters+1,2*size_clusters)+
               range(2*size_clusters+1,3*size_clusters)+
               range(3*size_clusters+1,4*size_clusters)+range(4*size_clusters+1,5*size_clusters)+range(5*size_clusters+1,6*size_clusters)] = 1
        return ad
    elif scenario is "BigExtraNewBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters)+
               range(size_clusters+15,2*size_clusters)+
               range(3*size_clusters+15,4*size_clusters),
            size_clusters*i:size_clusters*i+15] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraCurBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters)+range(2*size_clusters+15,3*size_clusters),size_clusters*i:size_clusters*i+15] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraMid":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters),(size_clusters*i):(size_clusters*i+15)] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraTop":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[size_clusters*i:size_clusters*i+15,
            range(15,size_clusters)+
            range(1*size_clusters+15,2*size_clusters)+
            range(2*size_clusters+15,3*size_clusters)+range(3*size_clusters+15,4*size_clusters)+
            range(4*size_clusters+15,5*size_clusters)+range(5*size_clusters+15,6*size_clusters)] = 1
        return ad
    else:
        raise LookupError("Invalid scenario")

def get_cluster_size(scenario, t_size=100):
    '''Find the size of each cluster for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param t_size: True size of each cluster. It is assumed the the true clusters are all the same size.
    '''
    t_n = 6 # true number of clusters
    if scenario is "Truth" or "ParentIs" in scenario:
        return np.repeat(t_size, t_n)
    elif scenario is "OneCluster":
        return np.array([t_n * t_size])
    elif "NCluster" in scenario:
        return np.repeat(1,t_n * t_size)
    elif "Split" in scenario:
        return np.concatenate((np.repeat(t_size,t_n-1),np.array([t_size / 2, t_size / 2])))
    elif "Merge" in scenario:
        return np.concatenate((np.array([2 * t_size]),np.repeat(t_size,t_n - 2)))
    elif "SmallExtra" in scenario:
        return np.concatenate((np.repeat(t_size - 1,t_n), np.array([t_n])))
    elif "BigExtra" in scenario:
        return np.concatenate((np.repeat(t_size - 15,t_n), np.array([15*t_n])))
    else:
        raise LookupError("Invalid scenario")

def get_cf(scenario, size_clusters):
    '''Find the cellular frequency for each cluster for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param size_cluster: the size of eah cluster in the given scenario
    '''
    return size_clusters / float(sum(size_clusters)) # default is frequency relative to size of cluster


def get_cellularity(scenario):
    '''Find the cellularity for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    '''
    t_cell = 0.7 # true cellularity
    if scenario is "Truth":
        return t_cell
    else:
        return t_cell + np.random.normal(scale=0.05)

# list of scenarios for testing the metric behaviour for sub-challenge 2 and 3 scores and the overall score
scenarios = ["Truth", "OneCluster", "NClusterOneLineage", "NClusterTwoLineages", "NClusterCorrectLineage",
                "ParentIsNieceWithChildren", "ParentIsSiblingWithChildren", "ParentIsCousin","ParentIsAunt", "ParentIsGrandparent", "ParentIsSibling",
                "BigExtraTop", "BigExtraMid", "BigExtraCurBot", "BigExtraNewBot",
                "SmallExtraTop", "SmallExtraMid", "SmallExtraNewBot", 'SmallExtraCurBot',
                "SplitClusterMidMultiChild", "SplitClusterMidOneChild", "SplitClusterBotSame", "SplitClusterBotDiff",
                "MergeClusterTop&Mid", "MergeClusterMid&BotMultiChild", "MergeClusterMid&BotOneChild","MergeClusterBot"]

def scoringtotal_behavior(verbose=False):
    '''Scoring behaviour of all subchallenge metrics together

    Attributes:
    :param verbose: boolean for whether to output details of the scoring metrics
    '''
    # True values for each attribute
    t_cell = get_cellularity("Truth")
    t_cluster_size = get_cluster_size("Truth")
    t_ncluster = len(t_cluster_size)
    t_cf = get_cf("Truth", t_cluster_size)
    t_1C = zip(t_cluster_size, t_cf)

    t_ccm = get_ccm("Truth")
    t_ad = get_ad("Truth")

    # list of predicted value for each scenario
    l_p_cell = list()
    l_p_ncluster = list()
    l_p_1C = list()
    l_p_ccm = list()
    l_p_ad = list()

    for scenario in scenarios:
        l_p_cell.append(get_cellularity(scenario))
        p_cluster_size = get_cluster_size(scenario)
        p_cf = get_cf(scenario, p_cluster_size)
        l_p_ncluster.append(len(p_cluster_size))
        l_p_1C.append(zip(p_cluster_size, p_cf))

        l_p_ccm.append(get_ccm(scenario))
        l_p_ad.append(get_ad(scenario))

    scores =  score_all(l_p_cell, t_cell,
              l_p_ncluster, t_ncluster,
              l_p_1C, t_1C,
              l_p_ccm, t_ccm,
              l_p_ad, t_ad, verbose=verbose)

    f = open(tsv_dir + "all_scores.csv", 'wb')
    wr = csv.writer(f)
    wr.writerow(scores.keys())
    wr.writerow(scores.values())
    f.close()

    return scores




if __name__ == '__main__':
    methods ={
        '1A':['abs', 'sqr'],
        '1B':['orig', 'normalized'],
        '1C':['abs', 'sqr'],
        '2':["orig",
            "sqrt",
            "pseudoV",
            "sym_pseudoV",
            "spearman",
            "pearson",
            "aupr",
            "mcc",
            "default"],
        '3': ['pseudoV',
              'sym_pseudoV',
              'pearson',
              'aupr',
              'sqrt',
              'mcc',
              'orig']
    }
    for m in methods['1A']:
        print 'Scoring 1A Behavior with method ' + m + '...'
        scoring1A_behavior(m)

    for m in methods['1B']:
        print 'Scoring 1B Behavior with method ' + m + '...'
        scoring1B_behavior(m)

    for m in methods['1C']:
        print 'Scoring 1C Behavior with method ' + m + '...'
        scoring1C_behavior(m)
    for m in methods['2']:
        print 'Scoring 2B Behavior with method ' + m + '...'
        scoring2B_behavior(method=m, verbose=True, tst_betas=True, tst_prob_mod_err=False, tst_prob_mod=False)
        scoring2A_behavior(method=m, verbose=True, tst_closest_reassign=False, tst_big_mat=False)


    print 'Scoring 3A Behavior...'
    scoring3A_behavior_all(verbose=True)

    print 'Scoring 3A Behavior using multiple metrics with different weights...'
    scoring3A_weight_behavior(verbose=True)
