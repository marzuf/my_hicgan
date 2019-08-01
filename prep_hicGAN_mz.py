	
import pandas as pd
import numpy as np

def prep_hicGAN_data(inFile, keepRatio, hicFormat, chrSize = None, binSize=None, nColToSkip=3, fieldSep="\t"):
    if hicFormat == "dixon":
        hic_DT = pd.read_csv(inFile, index_col=list(range(nColToSkip)), sep=fieldSep, header=None)
        hic_mat = hic_DT.values
       
    if hicFormat == "rao":
        #mat_dim = int(math.ceil(chrSize*1.0/binSize))     
        mat_dim = int(math.floor(chrSize*1.0/binSize) + 1)   
        hic_mat = np.zeros((mat_dim,mat_dim))
        for line in open(inFile).readlines():
            idx1, idx2, value = int(line.strip().split('\t')[0]),\
                                    int(line.strip().split('\t')[1]),\
                                        float(line.strip().split('\t')[2])                           
            hic_mat[int(idx1/binSize)][int(idx2/binSize)] = value
            hic_mat[int(idx2/binSize)][int(idx1/binSize)] = value
            # or to make it symmetric
            # hic_mat += hic_mat .T - np.diag(hic_mat .diagonal())
        
    assert hic_mat.shape[0] == hic_mat.shape[1]
             
    hic_mat_init = hic_mat.copy()             
             
    totSum = hic_mat.sum()
    countsToRemove = int((totSum - keepRatio*totSum)/2.0)
    # Return random integers from low (inclusive) to high (exclusive).
    rm_idxs = np.random.randint(hic_mat.shape[0], size=(2, countsToRemove))
    
    print(" > START REMOVING COUNTS")
    
        
    for k in range(countsToRemove):
        i = rm_idxs[0][k]
        j = rm_idxs[1][k]    
        hic_mat[i,j] = hic_mat[i,j] - 1
        hic_mat[j,i] = hic_mat[j,i] - 1

    print("totSum before: " + str(totSum))
    print("keepRatio * totSum: " + str(totSum*keepRatio))
    
    totSum_after = hic_mat.sum()
    
    print("totSum after: " + str(totSum_after))    	
        

    return hic_mat_init, hic_mat
        
