#!/usr/bin/env python
import os,sys,time
import math
#import matplotlib.pyplot as plt

import gdal,gdalconst
from scipy.stats import ks_2samp
import numpy as np
from numpy import *
import scipy as Sci
import scipy.linalg
import scipy.io as sio
from  scipy  import ndimage
import matplotlib.pyplot as plt
import struct,shutil
import gdal,gdalconst
from gdalconst import *
from scipy import stats
from scipy import signal
import multiprocessing
import re
import scipy.linalg
import scipy.io as sio
import linecache

global stBurst
global endBurst
global stswath
global endswath
global Range_look
global Azimuth_look
stBurst=1
endBurst=27
stswath=2
endswath=2
Range_look=20
Azimuth_look=4

outputWinFirstLine=1
outputWinLastLine=24100
outputWinFirstPix=1
outputWinLastPix=121000

#******************************************************************************
def Linux_cmd(python_instruction):
    #python_instruction='source '+work_dir+'/'+Configfile
    print 'python_instruction=',python_instruction
    with os.popen(python_instruction,"r") as p:
        cmd_profile_copy=p.read()
    #print 'cmd_profile_copy=',cmd_profile_copy       
    return cmd_profile_copy

#******************************************************************************
#******************************************************************************************#
def get_parameter(First_param,file_name,format_flag=1,Second_param=None,Third_param=None):
    Read_contine_flag=0
    orbit_info=""
    value=None
   
    for line in open(file_name):
        if format_flag==1:
            if not (line.find(First_param)):
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                return value 
        
        if format_flag==2:
            if  not (line.find(Second_param)):
                Read_contine_flag=1             
            if (Read_contine_flag==1) and (not (line.find(First_param))):  #Be careful
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))                             
            if  (not (line.find(Third_param))):  #Be careful                      
                Read_contine_flag=0              
                return value
               

        if format_flag==3:            
            if not (line.find(First_param)):
                index=line.find(':')
                pixel_time=(line[(index+1):].strip(' \n\t')).split(' ')[1].split(':')
                return pixel_time  


        if format_flag==4:   # for orbit

            
            if not (line.find(First_param)):
                index=line.find(':')
                
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1;
                #print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+1)):
                    return orbit_info  
        if format_flag==5:   # for corgistration

            
            if not (line.find(First_param)):
                index=line.find(':')
                
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1;
                print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+2)):
                    return orbit_info 
        if format_flag==6:   # for cpm

            
            if not (line.find(First_param)):
             
                Read_contine_flag=1;
                #print line
                continue                 
            if (Read_contine_flag>=1) :
               
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                nofCoef= 0.5*((int(Second_param)+1)**2+int(Second_param)+1)
                if (Read_contine_flag==int(nofCoef+1)):
                    return orbit_info 
         

#---------------------------------------------

#*********************************************
##################################################################################
def function_concatenate(filename_list):
    #print 'filename_list=',filename_list
    work_dir=os.getcwd()
    

    for cint_file in filename_list:
        print  cint_file.strip()
        file_fold_name=cint_file.strip()
 
        os.chdir(file_fold_name)
        print "file=",os.path.splitdrive(file_fold_name)
        print "file1=",os.path.splitext(file_fold_name)
        print "file2=",os.path.split(file_fold_name)[1]
        master_time=os.path.split(file_fold_name)[1].split('_')[0]
        slave_time=os.path.split(file_fold_name)[1].split('_')[1]
        print 'master_time=',master_time
        print 'slave_time=',slave_time
        #########################
        Output_fold    =work_dir+'/'+'Ras_copy'+'_'+str(master_time.strip())
        if not os.path.exists(Output_fold):
            os.makedirs(Output_fold)
        #########################
        python_instruction='rm -fr *.jpg'
        print 'python_instruction=',python_instruction
        
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################
        #########################
        python_instruction='concatenate_master.py  '+str(stBurst)+':'+str(endBurst)+' '+'master.res'+' '+\
        str(master_time.strip())+'_iw_'+str(stswath)+'.raw'
        print 'python_instruction=',python_instruction
        
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################
        python_instruction='sar_image_multilook.py'+' '+str(master_time.strip())+'_iw_'+str(stswath)+'.res'+' '+\
                           str(master_time.strip())+'_iw_'+str(stswath)+'.res'+' '+\
                           str(master_time.strip())+'_iw_'+str(stswath)+'.raw'+' '+'crop'+' '+\
                           str(Range_look)+' '+ str(Azimuth_look)
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################
        #########################

        python_instruction='concatenate_slave.py  '+str(stBurst)+':'+str(endBurst)+' '+'resample'+' '+\
        'slave_rsmp.raw'+' '+'False False False False'
        print 'python_instruction=',python_instruction
        
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################

        python_instruction='sar_image_multilook.py'+' '+str(master_time.strip())+'_iw_'+str(stswath)+'.res'+' '+\
                           'slave_rsmp.res'+' '+'slave_rsmp.raw'+' '+'resample'+' '+\
                           str(Range_look)+' '+ str(Azimuth_look)
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################
        #########################

        python_instruction='concatenate_ifg_res.py'+' '+str(stBurst)+':'+str(endBurst)+' '+'subtr_refdem'+' '+\
        'mergedIfgBursts_srd.raw'+' '+'False False False False'
        print 'python_instruction=',python_instruction
        
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################

        python_instruction='sar_image_multilook.py'+' '+str(master_time.strip())+'_iw_'+str(stswath)+'.res'+' '+\
                           'mergedIfgBursts_srd.res'+' '+'mergedIfgBursts_srd.raw'+' '+'subtr_refdem'+' '+\
                            str(Range_look)+' '+ str(Azimuth_look)
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
        #########################

        master_resfile=os.path.abspath(str(master_time.strip())+'_iw_'+str(stswath)+'.res')
        print 'master_resfile=',master_resfile

        First_line=get_parameter('First_line (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        Last_line =get_parameter('Last_line (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        First_pixel=get_parameter('First_pixel (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        Last_pixel=get_parameter('Last_pixel (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL')
 
    
        formatData1=get_parameter('Data_output_format',master_resfile,2,'*_Start_resample','* End_resample:_NORMAL')

        Image_width=int(Last_pixel)-int(First_pixel)+1
        print  'Image_width=',Image_width
        python_instruction='cpxfiddle -w '+str(Image_width)+' -e 0.2 -s 1.8 -q mixed -o sunraster  -c hsv -M 40/8 -f cr4 -l1 -p1 '+\
        ' -P'+str(Image_width)+' '+'mergedIfgBursts_srd.raw  > interferogram_mix_srd.ras'
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy

        jpg_name=str(master_time.strip())+'_'+str(slave_time.strip())+'.jpg'
        python_instruction='convert interferogram_mix_srd.ras '+' '+jpg_name
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy


        master_resfile=str(master_time.strip())+'_iw_'+str(stswath)+'_mul.res'
        print master_resfile

        First_line=get_parameter('First_line (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        Last_line =get_parameter('Last_line (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        First_pixel=get_parameter('First_pixel (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL') 
        Last_pixel=get_parameter('Last_pixel (w.r.t. original_image)',master_resfile,2,'*_Start_crop','* End_crop:_NORMAL')
 

        Crop_Image_width=int(Last_pixel)-int(First_pixel)+1
        print  'Crop_Image_width=',Crop_Image_width
        python_instruction='cpxfiddle -w '+str(Crop_Image_width)+' -e 0.2 -s 1.8 -q mixed -o sunraster  -c hsv -M 2/2 -f cr4 -l1 -p1 '+\
        ' -P'+str(Crop_Image_width)+' '+'mergedIfgBursts_srd_mul.raw  > interferogram_mix_mul.ras'
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy

        jpg_name=str(master_time.strip())+'_'+str(slave_time.strip())+'_iw_'+str(stswath)+'_mul.jpg'
        python_instruction='convert interferogram_mix_mul.ras '+' '+jpg_name
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy


        
        python_instruction="cp *.jpg"+' '+Output_fold
        cmd_profile_copy=Linux_cmd(python_instruction)
        print 'cmd_profile_copy=',cmd_profile_copy
##############
def main():
     
    cores=multiprocessing.cpu_count()
    print 'Total cores is ', cores
    #pool=multiprocessing.Pool(processes=(cores-1))
      
    #stBurst=1
    #endBurst=2
    #stswath=1
    #endswath=1
    ##############################################################
    ##Parameter of filter Windows
    #############################################################
   ##############################################################
    ##Parameter of filter Windows
    #############################################################
   
    work_dir=os.getcwd();
    #SBAS_work_dir=work_dir+'/SMALL_BASELINES'
    #print 'SBAS_work_dir=',SBAS_work_dir

    if os.path.exists(work_dir):    
        file_dir=os.listdir(work_dir)
    else:
        print "This path is wrong !!!"
        sys.exit(1)


   

    filename_list=[]
    for filename in file_dir:    
        #m_index=re.search('\d\d\d\d_\d\d\d\d',filename)
        m_index=re.search('20200616_20\d\d\d\d\d\d',filename)
        #M_index=re.search('\d_\d',filename)
        if m_index is not None:
        #print 'filename of Interferometry folder =',filename
            filename_list.append(work_dir+'/'+filename+'/'+filename+'_iw'+str(stswath)+'_fine')

    filename_list.sort()
     
    print 'filename_list=',filename_list
    fileinfo = open('Interfergrom_combination.list','w')
    for i in range(len(filename_list)):                      
        fileinfo.write(filename_list[i].strip()+'\n')
               
    fileinfo.close()
    #function_concatenate(filename_list)

   

    """
    filename_list=[]
    f=open('patch.list')
    patch_file_list=f.readlines()
    for patch_file in  patch_file_list:   
        print patch_file.strip()
        filename_list.append(patch_file.strip())
    #filename_list.sort()
    Num_patch =len(patch_file_list)
    #print 'Num_SLC =',Num_patch
    #pool.map(homofilter_KS,cint_minrefdem_file)
        #homofilter_KS(patch_Name)
    homofilter_KS()
    """   
        
    
if __name__== '__main__':
    main()
