import os, time, sys, math,random
import numpy as np
import hickle as hkl


# python data_split.py <INPUT_DIR> <FILE_PREFIX> <FILE_SUFFIX HIGH-RESOL> <FILE_SUFFIX LOW-RESOL> <OUT_DIR>
# python data_split.py INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt downsample16_merged_agg.txt INPUT_MATS

if len(sys.argv) != 6:
    print("ERROR: invalid number of arguments!")
    sys.exit(1)

input_dir = sys.argv[1]
file_prefix = sys.argv[2]
file_suffix_noDS = sys.argv[3]
file_suffix_ds = sys.argv[4]
output_dir = sys.argv[5]

out_dir = os.path.join(output_dir, file_prefix)
os.makedirs(out_dir, exist_ok=True)

if not os.path.exists(input_dir):
    print("ERROR: data path wrong !")
    sys.exit(1)

### HARD-CODED STUFF:
chromo_list = list(range(1,23))   #chr1-chr22
train_chromo_list = list(range(1,18))
test_chromo_list = list(range(18,23))

size_file = "all_chr_length.txt"
assert os.path.exists(size_file)
chrSizeDict = {item.split()[0]:int(item.strip().split()[1]) for item in open(size_file).readlines()}

bin_size = 10000
# set to work for 10kb data -> 200*10*1000 = 2 MB
maxDistBin = 200
# set to work for 10kb data -> 40*10kb = 400 kb
patchSizeBin = 40

print("!!! HARD-CODED: ")
print("chromo_list:")
print(chromo_list)
print("train_chromo_list:")
print(train_chromo_list)
print("test_chromo_list:")
print(test_chromo_list)
print("bin_size = " + str(bin_size))
print("maxDistBin = " + str(maxDistBin))
print("patchSizeBin = " + str(patchSizeBin))
#################################################################################################################################################################################
################################################################################################################################################################################# def hic_matrix_extraction
################################################################################################################################################################################# 

def hic_matrix_extraction(inDir,filePrefix, hr_fileSuffix, lr_fileSuffix, chrSizeDict, binSize, idxInBin=True):

    # PREPARE THE HIGH-RESOLUTION DATA
    print("*** PREPARING HIGH-RESOLUTION DATA")
    hr_contacts_dict={}
    for chrom in chromo_list:
        chromo = "chr" + str(chrom)
        chr_size = chrSizeDict[chromo]
        mat_dim = int(math.ceil(chr_size*1.0/binSize))
        print("> START " + chromo)
        print("... size = " + str(chr_size) + " ( " + str(mat_dim) + " bins)")
        hr_hic_file = '%s/%s_%s_%s'%(inDir, filePrefix, chromo, hr_fileSuffix)
        print(hr_hic_file)
        assert os.path.exists(hr_hic_file)
        print("... build HR matrix from:\t" + hr_hic_file)
        hr_contact_matrix = np.zeros((mat_dim,mat_dim))
        for line in open(hr_hic_file).readlines():
            idx1, idx2, value = int(line.strip().split('\t')[0]),int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])        
            #print("idx1 = " + str(idx1) + " - idx2 = " + str(idx2) + " - value = " + str(value))
            if not idxInBin:
                idx1 /= binSize
                idx2 /= binSize
            assert idx1 <= idx2
            assert idx1 <= mat_dim
            assert idx2 <= mat_dim
            hr_contact_matrix[idx1][idx2] = value
            
        # end-for iterating over interactions (lines in file)
        hr_contact_matrix+= hr_contact_matrix.T - np.diag(hr_contact_matrix.diagonal())
        hr_contacts_dict[chromo] = hr_contact_matrix
    # end-for iterating over chromosomes

    # PREPARE THE LOW-RESOLUTION DATA
    print("*** PREPARING LOW-RESOLUTION DATA")
    lr_contacts_dict={}
    for chrom in chromo_list:
        chromo = "chr" + str(chrom)
        chr_size = chrSizeDict[chromo]
        mat_dim = int(math.ceil(chr_size*1.0/binSize))
        print("> START " + chromo)
        print("... size = " + str(chr_size) + " ( " + str(mat_dim) + " bins)")

        lr_hic_file = '%s/%s_%s_%s'%(inDir, filePrefix, chromo, lr_fileSuffix)
        print(lr_hic_file)
        assert os.path.exists(lr_hic_file)
        print("... build LR matrix from:\t" + lr_hic_file)
        lr_contact_matrix = np.zeros((mat_dim,mat_dim))
        for line in open(lr_hic_file).readlines():
            idx1, idx2, value = int(line.strip().split('\t')[0]),int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])        
            #print("idx1 = " + str(idx1) + " - idx2 = " + str(idx2) + " - value = " + str(value))
            if not idxInBin:
                idx1 /= int(binSize)
                idx2 /= int(binSize)
            assert idx1 <= idx2
            assert idx1 <= mat_dim
            assert idx2 <= mat_dim
            lr_contact_matrix[idx1][idx2] = value

        # end-for iterating over iterations (lines in file)
        lr_contact_matrix+= lr_contact_matrix.T - np.diag(lr_contact_matrix.diagonal())
        lr_contacts_dict[chromo] = lr_contact_matrix

    # end-for iterating over chromosomes

    nb_hr_contacts={item:sum(sum(hr_contacts_dict[item])) for item in hr_contacts_dict.keys()}
    nb_lr_contacts={item:sum(sum(lr_contacts_dict[item])) for item in lr_contacts_dict.keys()}
    
    return hr_contacts_dict,lr_contacts_dict,nb_hr_contacts,nb_lr_contacts

#################################################################################################################################################################################
################################################################################################################################################################################# def crop_hic_matrix_by_chrom
################################################################################################################################################################################# 

def crop_hic_matrix_by_chrom(chrom,norm_type,patch_size,thresh_dist):
    #thred=2M/binSizeolution
    #norm_type=0-->raw count
    #norm_type=1-->log transformation
    #norm_type=2-->scaled to[-1,1]after log transformation, default
    #norm_type=3-->scaled to[0,1]after log transformation
    distance=[]
    crop_mats_hr=[]
    crop_mats_lr=[]    
    row,col = hr_contacts_norm_dict[chrom].shape
    if row<=thresh_dist or col<=thresh_dist:
        print('HiC matrix size wrong!')
        sys.exit(1)

    def quality_control(mat,quality_thresh=0.05):
        if len(mat.nonzero()[0]) < quality_thresh*mat.shape[0]*mat.shape[1]:
            return False
        else:
            return True
        
    for idx1 in range(0,row-patch_size,patch_size):
        for idx2 in range(0,col-patch_size,patch_size):
            if abs(idx1-idx2)<thresh_dist:
                if quality_control(lr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]):
                    distance.append([idx1-idx2,chrom])
                    if norm_type==0:
                        lr_contact = lr_contacts_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                        hr_contact = hr_contacts_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                    elif norm_type==1:
                        lr_contact = lr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                        hr_contact = hr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                    elif norm_type==2:
                        lr_contact_norm = lr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                        hr_contact_norm = hr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                        lr_contact = lr_contact_norm*2.0/max_lr_contact_norm[chrom]-1
                        hr_contact = hr_contact_norm*2.0/max_hr_contact_norm[chrom]-1
                    elif norm_type==3:
                        lr_contact_norm = lr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
                        hr_contact_norm = hr_contacts_norm_dict[chrom][idx1:idx1+patch_size,idx2:idx2+patch_size]
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

#################################################################################################################################################################################
################################################################################################################################################################################# def data_split
################################################################################################################################################################################# 



def data_split(chromo_list,norm_type):
    random.seed(100)
    distance_all=[]
    assert len(chromo_list) > 0
    hr_mats,lr_mats=[],[]
    for chrom in chromo_list:
        crop_mats_hr,crop_mats_lr,distance = crop_hic_matrix_by_chrom(chrom=chrom, norm_type=norm_type, patch_size=patchSizeBin ,thresh_dist=maxDistBin)
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


#################################################################################################################################################################################


print("> extract Hi-C data and prepare matrices... ")

hr_contacts_dict,lr_contacts_dict,nb_hr_contacts,nb_lr_contacts = hic_matrix_extraction(inDir=input_dir ,filePrefix=file_prefix, hr_fileSuffix=file_suffix_noDS, lr_fileSuffix=file_suffix_ds, chrSizeDict=chrSizeDict, binSize = bin_size)

max_hr_contact = max([nb_hr_contacts[item] for item in nb_hr_contacts.keys()])
max_lr_contact = max([nb_lr_contacts[item] for item in nb_lr_contacts.keys()])

#normalization

print("> extract normalized Hi-C data... ")

hr_contacts_norm_dict = {item:np.log2(hr_contacts_dict[item]*max_hr_contact/sum(sum(hr_contacts_dict[item]))+1) for item in hr_contacts_dict.keys()}
lr_contacts_norm_dict = {item:np.log2(lr_contacts_dict[item]*max_lr_contact/sum(sum(lr_contacts_dict[item]))+1) for item in lr_contacts_dict.keys()}

max_hr_contact_norm={item:hr_contacts_norm_dict[item].max() for item in hr_contacts_dict.keys()}
max_lr_contact_norm={item:lr_contacts_norm_dict[item].max() for item in lr_contacts_dict.keys()}


# WRITE NB CONTACT FILES
nb_hr_contactsFile = os.path.join(out_dir, out_dir + "_nb_hr_contacts.hkl")
hkl.dump(nb_hr_contacts, nb_hr_contactsFile)
print("... written: " + nb_hr_contactsFile)

nb_lr_contactsFile = os.path.join(out_dir, out_dir + "_nb_lr_contacts.hkl")
hkl.dump(nb_lr_contacts,nb_lr_contactsFile)
print("... written: " + nb_lr_contactsFile)


# WRITE MAX CONTACT FILES
max_hr_contact_normFile = os.path.join(out_dir, out_dir + "_max_hr_contact_norm.hkl")
hkl.dump(max_hr_contact_norm,max_hr_contact_normFile)
print("... written: " + max_hr_contact_normFile)

max_lr_contact_normFile = os.path.join(out_dir, out_dir + "_max_lr_contact_norm.hkl")
hkl.dump(max_lr_contact_norm,max_lr_contact_normFile)
print("... written: " + max_lr_contact_normFile)

# SPLIT AND WRITE TRAIN AND TEST DATA

print("> split data for training and testing... ")

hr_mats_train,lr_mats_train,distance_train = data_split(['chr%d'%idx for idx in train_chromo_list],norm_type=2)
hr_mats_test,lr_mats_test,distance_test = data_split(['chr%d'%idx for idx in test_chromo_list],norm_type=2)


train_data_file = os.path.join(out_dir, out_dir + "_train_data.hkl")
hkl.dump([lr_mats_train,hr_mats_train,distance_train], train_data_file)
print("... written: " + train_data_file)

test_data_file = os.path.join(out_dir, out_dir + "_test_data.hkl")
hkl.dump([lr_mats_test,hr_mats_test,distance_test], test_data_file)
print("... written: " + test_data_file)


#uncomment to save the raw readscount data
#hr_mats_train,lr_mats_train,distance_train = data_split(['chrom%d'%idx for idx in list(range(1,18))],norm_type=0)
#hr_mats_test,lr_mats_test,distance_test = data_split(['chrom%d'%idx for idx in list(range(18,23))],norm_type=0)
#hkl.dump([lr_mats_train,hr_mats_train,distance_train],'data/%s/train_data_raw_count.hkl'%cell)
#hkl.dump([lr_mats_test,hr_mats_test,distance_test],'data/%s/test_data_raw_count.hkl'%cell)
