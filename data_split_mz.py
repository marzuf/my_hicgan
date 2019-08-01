import os, time, pickle, sys, math,random
import numpy as np
import hickle as hkl
# hickle is hdf5 based clone of pickle
import pandas as pd
from prep_hicGAN_mz import *

# python data_split_mz.py output_datasplit

myargs = sys.argv

print(myargs)

#sys.exit()

if len(myargs) == 3:
    data_path = myargs[1].rstrip('/')
    save_dir = myargs[2].rstrip('/')
    if not os.path.exists(data_path):
        print('Data path wrong!')
        sys.exit()
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
elif len(myargs) == 2:
    save_dir = myargs[1].rstrip('/')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
else :
    print("ERROR # of arguments")
    sys.exit(1)

binSize = 50 * 1000
setDir = "/media/electron"
setDir = ""
dataFold = setDir + "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/KARPAS/DMSO/NORM"
filePref = "KARPAS_DMSO_"
fileSuf = "_" + str(binSize) + "_aggregCounts_ICE_matrix_pos.txt"

# MZ:
# lr = low resolution data (the downsampled data)
# hr = high resolution data 

#def hic_matrix_extraction(DPATH,res=10000,norm_method='NONE'):
def hic_matrix_extraction(mainFold, filePrefix, fileSuffix):

    hr_contacts_dict={}
    lr_contacts_dict={}

    for nchr in range(22):        
        chromo = "chr" + str(nchr+1)        
        print("> START " + chromo)
        
        hicFile = os.path.join(mainFold, filePrefix + chromo + fileSuffix)
        assert os.path.exists(hicFile )
        sizeFile = os.path.join(mainFold, chromo + ".size")
        assert os.path.exists(sizeFile )

        all_chrom_len = {item.split()[0]:int(item.strip().split()[1]) for item in open(sizeFile).readlines()}

        hr_contacts_dict[chromo], lr_contacts_dict[chromo] = prep_hicGAN_data(inFile = hicFile,keepRatio = 0.5,hicFormat = "dixon")

    nb_hr_contacts={item:sum(sum(hr_contacts_dict[item])) for item in hr_contacts_dict.keys()}
    nb_lr_contacts={item:sum(sum(lr_contacts_dict[item])) for item in lr_contacts_dict.keys()}
    
    return hr_contacts_dict,lr_contacts_dict,nb_hr_contacts,nb_lr_contacts

#hr_contacts_dict,lr_contacts_dict,nb_hr_contacts,nb_lr_contacts = hic_matrix_extraction(data_path)
hr_contacts_dict,lr_contacts_dict,nb_hr_contacts,nb_lr_contacts = hic_matrix_extraction(mainFold = dataFold, filePrefix = filePref, fileSuffix=fileSuf)

max_hr_contact = max([nb_hr_contacts[item] for item in nb_hr_contacts.keys()])
max_lr_contact = max([nb_lr_contacts[item] for item in nb_lr_contacts.keys()])

hr_contacts_norm_dict = {item:np.log2(hr_contacts_dict[item]*max_hr_contact/sum(sum(hr_contacts_dict[item]))+1) for item in hr_contacts_dict.keys()}
lr_contacts_norm_dict = {item:np.log2(lr_contacts_dict[item]*max_lr_contact/sum(sum(lr_contacts_dict[item]))+1) for item in lr_contacts_dict.keys()}

max_hr_contact_norm={item:hr_contacts_norm_dict[item].max() for item in hr_contacts_dict.keys()}
max_lr_contact_norm={item:lr_contacts_norm_dict[item].max() for item in lr_contacts_dict.keys()}

hkl.dump(nb_hr_contacts,'%s/nb_hr_contacts.hkl'%save_dir)
hkl.dump(nb_lr_contacts,'%s/nb_lr_contacts.hkl'%save_dir)

hkl.dump(max_hr_contact_norm,'%s/max_hr_contact_norm.hkl'%save_dir)
hkl.dump(max_lr_contact_norm,'%s/max_lr_contact_norm.hkl'%save_dir)



# MZ: if I understood correctly, function used to split the genome-wide Hi-C data into intra-chromo Hi-C
# "only kept intrachromosomal interactions in a reasonable genomic distance and then cropped the matrices into non-overlapping patches with size 400 x 400 kb^2"
# "each patch treated as an indiviudal sample which contains 40x40=1600 pixels"

# "only kept patches of which genomic distances between 2 loci is less than 2 Mb, thus filtering too far off-diagonal patches"

# MZ: size -> gives the size of the patch in # bins

# Returns a tuple of arrays, one for each dimension of a, containing the indices of the non-zero elements in that dimension. 
# The values in *a* are always tested and returned in row-major, C-style order. The corresponding non-zero values can be obtained with: a[nonzero(a)]

def crop_hic_matrix_by_chrom(chrom,norm_type,size=40 ,thred=200):

# MZ: the *thred* in "crop_hic_matrix_by_chrom" => abs(idx1-idx2) < thred => take only bin pair separated by < 2 Mb 
# MZ: the *thred* in "quality_control" => the ratio of nonzero values (take only if a certain ratio of interactions are nonzero)

    #thred=2M/resolution
    #norm_type=0-->raw count
    #norm_type=1-->log transformation
    #norm_type=2-->scaled to[-1,1]after log transformation
    #norm_type=3-->scaled to[0,1]after log transformation
    distance=[]
    crop_mats_hr=[]
    crop_mats_lr=[]    
    row,col = hr_contacts_norm_dict[chrom].shape
    if row<=thred or col<=thred:
        print('HiC matrix size wrong!')
        sys.exit()
    def quality_control(mat,thred=0.05):
        if len(mat.nonzero()[0])<thred*mat.shape[0]*mat.shape[1]: # MZ: len(mat.nonzero()[0]) => will give the # of values that are non zero; thred = the desired ratio of nonzero 
            return False
        else:
            return True
        
    for idx1 in range(0,row-size,size):
        for idx2 in range(0,col-size,size):
            if abs(idx1-idx2)<thred:    # MZ: if understood correctly, check that not separated by > 2 Mb (200 bins of 10 kb)
                if quality_control(lr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]):  # MZ: if understood correctly, check at least a certain ratio of interactions are non-zero
                    distance.append([idx1-idx2,chrom])
                    if norm_type==0:
                        lr_contact = lr_contacts_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        hr_contact = hr_contacts_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                    elif norm_type==1:
                        lr_contact = lr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        hr_contact = hr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                    elif norm_type==2:
                        lr_contact_norm = lr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        hr_contact_norm = hr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        lr_contact = lr_contact_norm*2.0/max_lr_contact_norm[chrom]-1
                        hr_contact = hr_contact_norm*2.0/max_hr_contact_norm[chrom]-1
                    elif norm_type==3:
                        lr_contact_norm = lr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        hr_contact_norm = hr_contacts_norm_dict[chrom][idx1:idx1+size,idx2:idx2+size]
                        lr_contact = lr_contact_norm*1.0/max_lr_contact_norm[chrom]
                        hr_contact = hr_contact_norm*1.0/max_hr_contact_norm[chrom]
                    else:
                        print('Normalization wrong!')
                        sys.exit()
                    
                    crop_mats_lr.append(lr_contact)
                    crop_mats_hr.append(hr_contact)
    crop_mats_hr = np.concatenate([item[np.newaxis,:] for item in crop_mats_hr],axis=0)
    crop_mats_lr = np.concatenate([item[np.newaxis,:] for item in crop_mats_lr],axis=0)
    return crop_mats_hr,crop_mats_lr,distance


def data_split(chrom_list,norm_type):
    random.seed(100)
    distance_all=[]
    assert len(chrom_list)>0
    hr_mats,lr_mats=[],[]
    for chrom in chrom_list:
        crop_mats_hr,crop_mats_lr,distance = crop_hic_matrix_by_chrom(chrom,norm_type,size=40 ,thred=200)
        distance_all+=distance
        hr_mats.append(crop_mats_hr)
        lr_mats.append(crop_mats_lr)
    hr_mats = np.concatenate(hr_mats,axis=0)
    lr_mats = np.concatenate(lr_mats,axis=0)
    hr_mats=hr_mats[:,np.newaxis]
    lr_mats=lr_mats[:,np.newaxis]
    hr_mats=hr_mats.transpose((0,2,3,1))
    lr_mats=lr_mats.transpose((0,2,3,1))
    return hr_mats,lr_mats,distance_all


hr_mats_train,lr_mats_train,distance_train = data_split(['chr%d'%idx for idx in list(range(1,18))],norm_type=0) # MZ: train on the chromo 1-17
hr_mats_test,lr_mats_test,distance_test = data_split(['chr%d'%idx for idx in list(range(18,23))],norm_type=0) # MZ: test on the chromo 18-23
hkl.dump([lr_mats_train,hr_mats_train,distance_train],'%s/train_data_raw_count.hkl'%save_dir)
hkl.dump([lr_mats_test,hr_mats_test,distance_test],'%s/test_data_raw_count.hkl'%save_dir)

hr_mats_train,lr_mats_train,distance_train = data_split(['chr%d'%idx for idx in list(range(1,18))],norm_type=1)
hr_mats_test,lr_mats_test,distance_test = data_split(['chr%d'%idx for idx in list(range(18,23))],norm_type=1)
hkl.dump([lr_mats_train,hr_mats_train,distance_train],'%s/train_data_log_trans.hkl'%save_dir)
hkl.dump([lr_mats_test,hr_mats_test,distance_test],'%s/test_data_log_trans.hkl'%save_dir)

hr_mats_train,lr_mats_train,distance_train = data_split(['chr%d'%idx for idx in list(range(1,18))],norm_type=2)
hr_mats_test,lr_mats_test,distance_test = data_split(['chr%d'%idx for idx in list(range(18,23))],norm_type=2)
hkl.dump([lr_mats_train,hr_mats_train,distance_train],'%s/train_data_log_trans_scaled_sym.hkl'%save_dir)
hkl.dump([lr_mats_test,hr_mats_test,distance_test],'%s/test_data_log_trans_scaled_sym.hkl'%save_dir)

hr_mats_train,lr_mats_train,distance_train = data_split(['chr%d'%idx for idx in list(range(1,18))],norm_type=3)
hr_mats_test,lr_mats_test,distance_test = data_split(['chr%d'%idx for idx in list(range(18,23))],norm_type=3)
hkl.dump([lr_mats_train,hr_mats_train,distance_train],'%s/train_data_log_trans_scaled_asym.hkl'%save_dir)
hkl.dump([lr_mats_test,hr_mats_test,distance_test],'%s/test_data_log_trans_scaled_asym.hkl'%save_dir)


























































