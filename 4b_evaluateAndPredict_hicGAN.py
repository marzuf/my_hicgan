import os, time, pickle, random, time, sys, math
from datetime import datetime
import numpy as np
from time import localtime, strftime
import logging, scipy
import hickle as hkl
import tensorflow as tf
import tensorlayer as tl
from tensorlayer.layers import *
import matplotlib.pyplot as plt
from skimage.measure import compare_mse
from skimage.measure import compare_ssim

import datetime
import re

startTime = datetime.datetime.now()

# python 4b_evaluateAndPredict_hicGAN.py <TEST_DATA_FILE> <OUTPUT_FOLDER>
# python 4b_evaluateAndPredict_hicGAN.py INPUT_HKL_chr1-4/KARPAS_DMSO/KARPAS_DMSO_test_data.hkl TRAIN_hicGAN_chr1-4/KARPAS_DMSO/checkpoint EVAL_PRED_HICGAN_chr1-4/KARPAS_DMSO

# CHANGE MZ
#os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[1]
#model_path=sys.argv[2].rstrip('/')
#cell=sys.argv[3]

#modeldir = "/media/electron/mnt/etemp/marie/hicGAN/hicgan_virtual/KARPAS_data/TRAIN_hicGAN_chr1-4/KARPAS_DMSO/checkpoint"#g_hicgan_720.npz"

assert len(sys.argv) == 4

test_data_file=sys.argv[1].rstrip('/')
model_dir = sys.argv[2]  # the checkpoint folder
output_folder=sys.argv[3]

if not os.path.exists(test_data_file):
    print("ERROR: invalid path input test_data_file !")
    sys.exit(1)


if not os.path.exists(model_dir):
    print("ERROR: invalid path to model_dir !")
    sys.exit(1)

os.makedirs(output_folder, exist_ok=True)

# TRAIN_hicGAN_chr1-4/KARPAS_DMSO/checkpoint/g_hicgan_2765.npz
all_check_files = [x for x in os.listdir(model_dir) if x.endswith(".npz") and x.startswith("g_hicgan")]

if len(all_check_files) <= 0:
    print("ERROR: no checkpoint files found")
    sys.exit(1)

# HARD-CODED SETTINGS
pred_batchSize = 64
matching_pattern = 'g_hicgan_(\d+).npz'
p = re.compile(matching_pattern)
minIdx = 100
idx_step = 100

modelSelectFile = os.path.join(output_folder, "model_selection_logFile.txt")
print("... write model selection in:\t" + modelSelectFile)

logFile = os.path.join(output_folder, "eval_and_pred_hicgan_logFile.txt")
print("... write logs in:\t" + logFile)

mylog = open(logFile,"a+") 
mylog.write("! HARD-CODED SETTINGS !")
mylog.write("> pred_batchSize\t=\t" + str(pred_batchSize))
mylog.write("> matching_pattern\t=\t" + str(matching_pattern))
mylog.write("> minIdx\t=\t" + str(minIdx))
mylog.write("> idx_step\t=\t" + str(idx_step))
mylog.close() 

def calculate_psnr(mat1,mat2):
    data_range=np.max(mat1)-np.min(mat1)
    err=compare_mse(mat1,mat2)
    return 10 * np.log10((data_range ** 2) / err)

def calculate_ssim(mat1,mat2):
    data_range=np.max(mat1)-np.min(mat1)
    return compare_ssim(mat1,mat2,data_range=data_range)

def hicGAN_g(t_image, is_train=False, reuse=False):
    w_init = tf.random_normal_initializer(stddev=0.02)
    b_init = None  # tf.constant_initializer(value=0.0)
    g_init = tf.random_normal_initializer(1., 0.02)
    with tf.variable_scope("hicGAN_g", reuse=reuse) as vs:
        n = InputLayer(t_image, name='in')
        n = Conv2d(n, 64, (3, 3), (1, 1), act=tf.nn.relu, padding='SAME', W_init=w_init, name='n64s1/c')
        temp = n
        # B residual blocks
        for i in range(5):
            nn = Conv2d(n, 64, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n64s1/c1/%s' % i)
            nn = BatchNormLayer(nn, act=tf.nn.relu, is_train=is_train, gamma_init=g_init, name='n64s1/b1/%s' % i)
            nn = Conv2d(nn, 64, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n64s1/c2/%s' % i)
            nn = BatchNormLayer(nn, is_train=is_train, gamma_init=g_init, name='n64s1/b2/%s' % i)
            nn = ElementwiseLayer([n, nn], tf.add, name='b_residual_add/%s' % i)
            n = nn
        n = Conv2d(n, 64, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, b_init=b_init, name='n64s1/c/m')
        n = BatchNormLayer(n, is_train=is_train, gamma_init=g_init, name='n64s1/b/m')
        n = ElementwiseLayer([n, temp], tf.add, name='add3')
        # B residual blacks end. output shape: (None,w,h,64)
        n = Conv2d(n, 128, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n128s1/1')  # changed MZ: still written n256 in original file
        #n = Conv2d(n, 256, (3, 3), (1, 1), act=None, padding='SAME', W_init=w_init, name='n256s1/2')
        n = Conv2d(n, 1, (1, 1), (1, 1), act=tf.nn.tanh, padding='SAME', W_init=w_init, name='out')
        return n
    
t_image = tf.placeholder('float32', [None, None, None, 1], name='image_input')
net_g = hicGAN_g(t_image, is_train=False, reuse=False)   

def hicGAN_predict(batch):
    sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=False))
    tl.layers.initialize_global_variables(sess)
    tl.files.load_and_assign_npz(sess=sess, name=model_dir, network=net_g)
    out = np.zeros(lr_mats_test.shape)
    for i in range(int(out.shape[0]/batch)):
        out[batch*i:batch*(i+1)] = sess.run(net_g.outputs, {t_image: lr_mats_test[batch*i:batch*(i+1)]})
    out[batch*(i+1):] = sess.run(net_g.outputs, {t_image: lr_mats_test[batch*(i+1):]})
    return out

def evaluate(gan_idx,batch=64):
    sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=False))
    tl.layers.initialize_global_variables(sess)
    model_file = '%s/g_hicgan_%d.npz'%(model_dir,gan_idx)
    print("... start evaluating: " + model_file)
    tl.files.load_and_assign_npz(sess=sess, name=model_file, network=net_g)
    out = np.zeros(lr_mats_test.shape)
    for i in range(int(out.shape[0]/batch)):
        out[batch*i:batch*(i+1)] = sess.run(net_g.outputs, {t_image: lr_mats_test[batch*i:batch*(i+1)]})
    out[batch*(i+1):] = sess.run(net_g.outputs, {t_image: lr_mats_test[batch*(i+1):]})
    return out


matching_st = [int(p.search(s).group(1)) for s in  all_check_files]
max_model_idx = max(matching_st)
min_model_idx = min([idx for idx in matching_st if idx >= minIdx])

# INPUT_HKL_chr1-4/KARPAS_DMSO/KARPAS_DMSO_test_data.hkl
# lr_mats_test,hr_mats_test,_ = hkl.load("INPUT_HKL_chr1-4/KARPAS_DMSO/KARPAS_DMSO_test_data.hkl")
lr_mats_test,hr_mats_test,_ = hkl.load(test_data_file)

f_out = open(modelSelectFile,'w')
mse_median_list = []
#model search index
for i in range(min_model_idx,max_model_idx,idx_step):
    out = evaluate(i)
    hkl.dump(out, "out_tmp.hkl")
    #mse = np.median(map(compare_mse,out[:,:,:,0],hr_mats_test[:,:,:,0]))
    mse = np.median(list(map(compare_mse,out[:,:,:,0],hr_mats_test[:,:,:,0])))
    mse_median_list.append(mse)
    f_out.write('%d\t%.5f\t%s\n'%(min_model_idx+idx_step*i,mse,model_dir))
best_idx = min_model_idx+idx_step*(mse_median_list.index(min(mse_median_list)))
f_out.close()

sr_mats_pre = evaluate(best_idx)

txt = "The best model is gan_%d" % best_idx
print(txt)
mylog = open(logFile,"a+") 
mylog.write(txt)
mylog.close() 


#uncomment to test with orginal readscount data
#def trans_norm2readscount_hr(chrom,mat):
#    max_hr_contact = max([nb_hr_contacts[item] for item in nb_hr_contacts.keys()])
#    mat_hat=(mat+1)*0.5*max_hr_contact_norm[chrom]
#    return (np.exp2(mat_hat)-1)*nb_hr_contacts[chrom]/max_hr_contact
#    
#max_hr_contact_norm = hkl.load('data/%s/max_hr_contact_norm.hkl')
#max_lr_contact_norm = hkl.load('data/%s/max_lr_contact_norm.hkl')
#nb_hr_contacts = hkl.load('data/%s/nb_hr_contacts.hkl')
#nb_lr_contacts = hkl.load('data/%s/nb_lr_contacts.hkl')
#chrom_list = [item[1] for item in distance_all]
#sr_mats_pre_readscount = map(trans_norm2readscount_hr,chrom_list,sr_mats_pre[:,:,:,0])
#_,hr_mats_test_readscount,_ = hkl.load('data/GM12878/test_data_raw_count.hkl')
#mse_hicGAN=map(compare_mse,hr_mats_test_readscount[:,:,:,0],sr_mats_pre_readscount)
#psnr_hicGAN=map(calculate_psnr,hr_mats_test_readscount[:,:,:,0],sr_mats_pre_readscount)
#ssim_hicGAN=map(calculate_ssim,hr_mats_test_readscount[:,:,:,0],sr_mats_pre_readscount)
#print 'mse_hicGAN:%.5f'%np.median(mse_hicGAN)
#print 'psnr_hicGAN:%.5f'%np.median(psnr_hicGAN)
#print 'ssim_hicGAN:%.5f'%np.median(ssim_hicGAN)

outfile_pred = output_folder + "hicGAN_predicted.npz"
sr_mats_pre = hicGAN_predict(batch=pred_batchSize)
np.savez(outfile_pred)
print("... hicGAN predictions written in " + outfile_pred)
    
mse_hicGAN_norm = map(compare_mse, hr_mats_test[:,:,:,0],sr_mats_pre[:,:,:,0])
psnr_hicGAN_norm = map(calculate_psnr, hr_mats_test[:,:,:,0],sr_mats_pre[:,:,:,0])
ssim_hicGAN_norm = map(calculate_ssim, hr_mats_test[:,:,:,0],sr_mats_pre[:,:,:,0])
#txt1 = 'mse_hicGAN_norm:%.5f' % np.median(mse_hicGAN_norm)
#txt2 = 'psnr_hicGAN_norm:%.5f' % np.median(psnr_hicGAN_norm)
#txt3 = 'ssim_hicGAN_norm:%.5f' % np.median(ssim_hicGAN_norm)
txt1 = 'mse_hicGAN_norm:%.5f' % np.median(list(mse_hicGAN_norm))
txt2 = 'psnr_hicGAN_norm:%.5f' % np.median(list(psnr_hicGAN_norm))
txt3 = 'ssim_hicGAN_norm:%.5f' % np.median(list(ssim_hicGAN_norm))
print(txt1)
print(txt2)
print(txt3)
mylog = open(logFile,"a+") 
mylog.write(txt1)
mylog.write(txt2)
mylog.write(txt3)
mylog.close() 

print("... written: " + logFile)




################################################################################################
################################################################################################ DONE
################################################################################################
endTime = datetime.datetime.now()
print("*** DONE")
print(str(startTime) + " - " + str(endTime))

mylog = open(logFile,"a+") 
mylog.write(str(startTime) + "\n" + str(endTime))
mylog.close() 

print("... logs written in: " + logFile)


  

