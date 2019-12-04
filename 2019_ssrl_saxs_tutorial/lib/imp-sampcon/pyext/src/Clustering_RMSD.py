import os,sys,random,numpy,math
import pickle

import scipy as sp
from scipy import spatial
import scipy.stats
import time

def get_sample_identity(idfile_A, idfile_B):
    # whether a run belongs to run1 or run2
    sampleA_models=[]
    sampleB_models=[]

    with open(idfile_A, 'r') as f:
        for line in f:
            mod = line.split("/")[-1]
            sampleA_models.append(int(mod.strip("\n").split(" ")[1]))
    f.close()
    
    with open(idfile_B, 'r') as f:
        for line in f:
            mod = line.split("/")[-1]
            sampleB_models.append(int(mod.strip("\n").split(" ")[1]))
    return sampleA_models, sampleB_models


def get_cutoffs_list(distmat, gridSize):

    mindist = distmat.min()
    maxdist = distmat.max()

    print("Minimum and maximum pairwise model distances:",mindist, maxdist)
    cutoffs=numpy.arange(mindist+gridSize,maxdist,gridSize)
    return cutoffs


def precision_cluster_test(distmat,numModels,rmsd_cutoff):
    #STEP 2. Populate the neighbors ofa given model
    print("The new way")
    stime = time.time()
    neighbors=[]
    for count in range(numModels):
        neighbors.append([count])  # model is a neighbor of itself
    inds = numpy.argwhere(distmat <= rmsd_cutoff) # set of i,j indices that pass
    print("There are",len(inds),"neighbors")
    for x in inds:
        i = x[0]
        j = x[1]
        if i>j:
            neighbors[i].append(j)
            neighbors[j].append(i)
    print("Done the new way in",time.time()-stime)

    print("Now the old way")
    stime = time.time()
    neighbors2=[]
    for count in range(numModels):
        neighbors2.append([count])  # model is a neighbor of itself

    for i in range(numModels-1):
        if i%1000 == 0:
            print("On model", i, "of", numModels)
        for j in range(i+1,numModels):
            if distmat[i][j]<=rmsd_cutoff:
                neighbors2[i].append(j)
                neighbors2[j].append(i)

    print("Done the old way in",time.time()-stime)

    print("Now check similarity")
    for i in range(numModels):
         if set(neighbors[i]) != set(neighbors2[i]):
             print("Uh oh!!", i, neighbors[i], neighbors2[i])
             print("New way dists", [distmat[i][j] for j in neighbors[i]])
             print("Old way dists", [distmat[i][j] for j in neighbors2[i]])
             print("Clustering cutoff", rmsd_cutoff)
             exit()
    exit()

def precision_cluster(distmat,numModels,rmsd_cutoff):
    #STEP 2. Populate the neighbors ofa given model
    
    # neighbors is a list of lists of integers, where each integer is the index of the model
    # in all_models (A + B)
    neighbors=[]
    for count in range(numModels):
        neighbors.append([count])  # model is a neighbor of itself

    inds = numpy.argwhere(distmat <= rmsd_cutoff) # set of i,j indices that pass
    
    for x in inds:
        i = x[0]
        j = x[1]
        # Only count each contribution once
        if i>j:
            neighbors[i].append(j)
            neighbors[j].append(i)

    #STEP 3. Get the weightiest cluster, and iterate
    unclustered=[]
    boolUnclustered=[]
    for i in range(numModels):
        unclustered.append(i)
        boolUnclustered.append(True)

    cluster_members=[] # list of lists : one list per cluster
    cluster_centers=[]

    while len(unclustered)>0:
        # get cluster with maximum weight
        max_neighbors=0
        currcenter=-1
        for eachu in unclustered:  # if multiple clusters have same maxweight this tie is broken arbitrarily! 
            if len(neighbors[eachu])>max_neighbors:
                max_neighbors=len(neighbors[eachu])
                currcenter=eachu   
   
        #form a new cluster with u and its neighbors
        cluster_centers.append(currcenter)
        cluster_members.append([n for n in neighbors[currcenter]]) 

        #update neighbors 
        for n in neighbors[currcenter]:
            #removes the neighbor from the pool
            unclustered.remove(n) #first occurence of n is removed. 
            boolUnclustered[n]=False # clustered

        for n in neighbors[currcenter]:
            for unn in neighbors[n]: #unclustered neighbor
                if not boolUnclustered[unn]:
                    continue
                neighbors[unn].remove(n)
    
    return cluster_centers, cluster_members

def get_contingency_table(num_clusters,cluster_members,n1_models,n2_models):
    full_ctable=numpy.zeros((num_clusters,2))
   
    # enumerate over all clusters
    # ic is cluster number, cluster is the cluster_memer idices for models
    # in that cluster
    for ic,cluster in enumerate(cluster_members):

	# Ok, get this index from all_models
        for member in cluster:
            if member < n1_models:
                #print("run1", ic, cluster, full_ctable[ic][0])
                full_ctable[ic][0]+=1.0
            else:
                #print("run2", ic, cluster, full_ctable[ic][1])
                full_ctable[ic][1]+=1.0

    ## now normalize by number of models in each run
    numModelsRun1 = float(numpy.sum(full_ctable,axis=0)[0])
    numModelsRun2 = float(numpy.sum(full_ctable,axis=0)[1])

    reduced_ctable=[]
    retained_clusters=[]
    
    for i in range(num_clusters):
        #print("  ",i, [full_ctable[i][0],full_ctable[i][1]])
        if full_ctable[i][0]<=10.0 or full_ctable[i][1]<=10.0:
            #if full_ctable[i][0]<=0.10*numModelsRun1 and full_ctable[i][1] <= 0.10*numModelsRun2:
            continue
        reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
        retained_clusters.append(i)

    return numpy.array(reduced_ctable),retained_clusters

def test_sampling_convergence(contingency_table,total_num_models):
    if len(contingency_table)==0:
        return (0.0,1.0)

    ct = numpy.transpose(contingency_table)

    [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)
    if dof==0.0:
        cramersv=0.0
    else:
        cramersv=math.sqrt(chisquare/float(total_num_models))
        
    return(pvalue,cramersv)

def percent_ensemble_explained(ctable,total_num_models):
    if len(ctable)==0:
        return 0.0
    percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)
    return percent_clustered

def get_clusters(cutoffs_list, distmat_full, total_num_models, n1_models, n2_models, sysname):
    #Do Clustering on a Grid
    pvals=[]
    cvs=[]
    percents=[]
    f1=open("%s.ChiSquare_Grid_Stats.txt" % sysname, 'w+')
    for c in cutoffs_list:
        cluster_centers,cluster_members=precision_cluster(distmat_full,
                                                          total_num_models, c)
        #print("C, CC", c, cluster_centers)
        ctable,retained_clusters=get_contingency_table(len(cluster_centers), cluster_members, n1_models, n2_models)

        (pval,cramersv)=test_sampling_convergence(ctable, total_num_models)

        percent_explained= percent_ensemble_explained(ctable, total_num_models)

        pvals.append(pval)
        cvs.append(cramersv)
        percents.append(percent_explained)
        
        print(str(c)+" "+str(pval)+" "+str(cramersv)+" "+str(percent_explained), file=f1)

    return pvals, cvs, percents

def get_sampling_precision(cutoffs_list, pvals, cvs, percents):
    sampling_precision=max(cutoffs_list)
    pval_converged=0.0
    cramersv_converged=1.0
    percent_converged=0.0

    for i in range(len(cutoffs_list)):
        if percents[i]>80.0:
            if pvals[i]>0.05 or cvs[i]<0.10:
                if sampling_precision>cutoffs_list[i]:
                    sampling_precision=cutoffs_list[i]
                    pval_converged=pvals[i]
                    cramersv_converged=cvs[i]
                    percent_converged=percents[i]
        else:
            sampling_precision=max(cutoffs_list)

    return sampling_precision,pval_converged,cramersv_converged,percent_converged
