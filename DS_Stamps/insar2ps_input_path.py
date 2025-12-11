#!/usr/bin/env python
import multiprocessing
import numpy as np
from scipy.optimize import minimize
import os,sys,time
import numpy as np
from scipy.stats import ks_2samp
import math
import numpy as np
from numpy import *
import numpy
#import matplotlib.pyplot as plt
import scipy
import gdal,gdalconst
from gdalconst import *
from struct import *
import datetime
import re
import shutil
#version 2021 01 16
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
def Linux_cmd(python_instruction):
    #python_instruction='source '+work_dir+'/'+Configfile
    print 'python_instruction=',python_instruction
    with os.popen(python_instruction,"r") as p:
        cmd_profile_copy=p.read()
    #print 'cmd_profile_copy=',cmd_profile_copy       
    return cmd_profile_copy

#*********************************************
Amplitude_dispersion="0.49"
#######################################################

Configfile='StaMPS_CONFIG.bash'

##N1datadir="/home/wuwenhao/data_disk/kunming/data_processing"
master_date="20191002"
Geo_date="20191008"
swath_num='2'
N1datadir="/home/wuwenhao/disk_external/Guangxi/data_processing"
work_dir=os.getcwd();
print 'work_dir=' ,work_dir
##################################################################

#####provide convenience for users
##################################################################
python_instruction='source '+work_dir+'/'+Configfile
cmd_turn=Linux_cmd(python_instruction)
print 'cmd_turn=',cmd_turn

os.chdir(work_dir)

#SAR_Interfer_Data_path=N1datadir
#print 'SAR_Interfer_Data_path=',SAR_Interfer_Data_path

#if os.path.exists(SAR_Interfer_Data_path):    
#    file_dir=os.listdir(SAR_Interfer_Data_path)
#else:
#    print "This path is wrong !!!"
#    sys.exit(1)
filename_list=[]
large_mini_stack=[]
original_large_mini_stack=[]
Interfergrom_combination_info=[]

for filename in open('Int_Data_merge.list'):
    print 'filename=',filename
    swath_Number=os.path.dirname(filename.strip().split('iw')[1])
    master_date=os.path.basename(filename.strip()).split('_')[0].strip()
    slave_date=os.path.basename(filename.strip()).split('_')[1].strip()
    filename_list.append(slave_date.strip()) 
    Interfergrom_combination_info.append(filename.strip())
filename_list.append(master_date)
print 'filename_list=',filename_list
#time.sleep(10000000000)
###############################################################
Geo_date=filename_list[0].strip()
for temp_file in filename_list:    
    if abs(int(master_date)-int(temp_file)) > abs(int(master_date)-int(Geo_date)):
        Geo_date=temp_file.strip()
print 'Geo_date=',Geo_date 


Master_Date_fold='INSAR_'+master_date.strip()
if not os.path.exists(Master_Date_fold):
    os.makedirs(Master_Date_fold)

#mkdir $workdir/SLC
#cp $N1datadir/timetableN1.txt $workdir/SLC/make_slcs.list
Insar_master_check=0
#time.sleep(10000000000)
##################################################################
for insar_file_name in Interfergrom_combination_info:
    print '\n'
    print 'insar_file_name=',insar_file_name
    file_time    =os.path.basename(insar_file_name).split('_')[1].strip()
    print 'file_time=',file_time
    if int(file_time.strip())!=int(master_date.strip()):
        original_fold=insar_file_name.strip()

        print 'original_fold=',original_fold
        target_fold  =work_dir+'/'+Master_Date_fold+'/'+file_time
        print 'work_dir=',work_dir
        print 'file_time=',file_time
        print 'target_fold=',target_fold
        if not os.path.exists(target_fold):
            os.makedirs(target_fold)
        else:
            screen_print=target_fold+' is alread here'
            print screen_print 
            print '\n'

        [dirname,filename]=os.path.split(original_fold.strip())
        print(dirname,"\n",filename)
        
        inner_original_fold=''
        #shutil.copy()
        #time.sleep(100000000000)
        print 'original_fold=',original_fold
        print 'target_fold=',target_fold
        [dirname,filename]=os.path.split(original_fold.strip())
        print(dirname,"\n",filename)
        #inner_original_fold=filename+'_iw'+swath_num+'_fine'
        inner_original_fold=''
        original_source=dirname+'/'+file_time+'/'+'Crop'+'/'+'resample_crop.res' 
        target_source  =target_fold+'/'+'slave.res'
        shutil.copy(original_source,target_source)
        #time.sleep(100000000000)
        #original_source_res=original_fold+'/'+inner_original_fold+'/'+'slave_rsmp_crop.res' 
        #target_source_res  =target_fold+'/'+'slave.res'
        #shutil.copy(original_source_res,target_source_res)
        #time.sleep(100000000000)
        slave_rsmp_raw=get_parameter('Data_output_file',original_source,2,'*_Start_resample','* End_resample:_NORMAL') 
        original_slave_data=get_parameter('Data_output_file',original_source,2,'*_Start_crop','* End_crop:_NORMAL')
        original_source=dirname+'/'+file_time+'/'+'Crop'+'/'+'resample_crop.raw' 
        print 'original_source=',original_source
        target_source  =target_fold+'/'+'resample_crop.raw'
        print 'target_source =',target_source 
        target_source_alias  =target_fold+'/'+'slave_res.slc'
        print 'target_source_alias=',target_source_alias
        #shutil.copy(original_source,target_source)
        if not (os.path.isfile(target_source)):        
            os.symlink(original_source,target_source)
        else:
            print  target_source+' is already here!!!\n'

        if os.path.isfile(original_source):
            if not (os.path.isfile(target_source_alias)):        
                os.symlink(original_source,target_source_alias)
            else:
                print  target_source_alias+' is already here!!!\n'

        ########################################################################
        #slave_rsmp_raw=get_parameter('Data_output_file',original_source_res,2,'*_Start_resample','* End_resample:_NORMAL') 
        """
        original_source=original_fold+'/'+inner_original_fold+'/'+'merged_h2ph_srd_crop.raw' 
        target_source  =target_fold+'/'+'merged_h2ph_srd.raw'
        
        #shutil.copy(original_source,target_source)
        if not (os.path.isfile(target_source)):        
            os.symlink(original_source,target_source)
        else:
            print  target_source+' is already here!!!\n'
        """
        
        ########################################################################
        
        print 'original_slave_data=',original_slave_data
        original_source=dirname+'/'+file_time+'/'+'Burst'+'/'+original_slave_data
      
        target_source  =target_fold+'/'+original_slave_data

        print 'target_source =',target_source 
        if os.path.isfile(original_source):
            if not (os.path.isfile(target_source)) :        
                os.symlink(original_source,target_source)
                print "link finished!!!\n\n\n"
            else:
                print  target_source+' is already here!!!\n'
    
        #time.sleep(100000000000)
        """
        original_source_res=original_fold+'/'+inner_original_fold+'/'+master_date+'_iw_'+swath_num+'_crop.res'
        target_source_res  =target_fold+'/'+'master.res'
        if not (os.path.isfile(target_source_res)):        
            shutil.copy(original_source_res,target_source_res)
            print '\n'
        else:
            print  target_source+' is already here!!!\n'  
 
        """

original_source_res=dirname+'/'+master_date.strip()+'/'+'Crop'+'/'+master_date.strip()+'_'+'iw'+swath_Number+'_crop.res'
target_source_res  =work_dir+'/'+Master_Date_fold+'/'+'master.res'
if not (os.path.isfile(target_source_res)):        
    shutil.copy(original_source_res,target_source_res)
    print '\n'
else:
    print  target_source+' is already here!!!\n'
    
        
sourceText   = get_parameter('Data_output_file',target_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
replaceText  = work_dir+'/'+Master_Date_fold+'/'+master_date+'_crop.slc'
print 'sourceText   = ',sourceText
print 'replaceText  =',replaceText
inputStream  = open(target_source_res,'r')
textStream   = inputStream.read()
inputStream.close()
outputStream = open(target_source_res,"w")
outputStream.write(textStream.replace(sourceText, replaceText))
outputStream.close()




Insar_master_check=0
if Insar_master_check==0:    
    #original_master_data=get_parameter('Data_output_file',target_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
    #print 'original_master_data=',original_master_data
    original_source     =dirname+'/'+master_date.strip()+'/'+'Crop'+'/'+master_date.strip()+'_'+'iw'+swath_Number+'_crop.raw'
    #original_source     =original_fold+'/'+inner_original_fold+'/'+os.path.basename(original_master_data.strip())
    #original_source=original_master_data
    print 'original_source=',original_source 
    target_source       =work_dir+'/'+Master_Date_fold+'/'+master_date+'_crop.slc'
    print 'target_source=',target_source

if Insar_master_check==1:
    original_master_data    =os.path.dirname(insar_file_name)+'/'+master_date+'/'+'resample_ministack.raw'
    print 'original_master_data=',original_master_data
    original_source   =original_master_data
    original_file_name=get_parameter('Data_output_file',original_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
    target_source     =work_dir+'/'+Master_Date_fold+'/'+os.path.basename(original_file_name)

    sourceText   = get_parameter('Data_output_format',target_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
    replaceText  = 'complex_real4'
    print 'sourceText   = ',sourceText
    print 'replaceText  =',replaceText
    inputStream  = open(target_source_res,'r')
    textStream   = inputStream.read()
    inputStream.close()
    outputStream = open(target_source_res,"w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()

if not (os.path.isfile(target_source)):        
    os.symlink(original_source,target_source)
    print '\n\n\n'
    print 'original_source=',original_source
    print 'target_source=',target_source 
    print '\n\n\n' 
else:
    print  target_source+' is already here!!!\n'
    #os.symlink(original_source,target_source)
    
    # replace crop tag in result file



print 'target_source_res=',target_source_res
python_instruction='python master_int2float.py '+target_source_res
#print 'python_instruction=',python_instruction
cmd_profile_copy=Linux_cmd(python_instruction)
#with os.popen(python_instruction,"r") as p:
#    cmd_profile_copy=p.read()
print 'cmd_profile_copy=',cmd_profile_copy  




python_file=open('timetableN1.txt','w')
for insar_file_name in filename_list:
    #print 'insar_file_name=',insar_file_name
    file_time    =insar_file_name.split('_')[-1].strip()
    print 'file_time=',file_time
    if file_time.strip()!=master_date.strip():
        original_fold_res=work_dir+'/'+Master_Date_fold+'/'+'master.res'
        target_fold_res  =work_dir+'/'+Master_Date_fold+'/'+file_time+'/'+'master.res'
        shutil.copy(original_fold_res,target_fold_res )
        python_file.writelines(file_time+'\n')
python_file.close()



os.chdir(Master_Date_fold)

#print os.path.split(original_master_data)
#time_tables=os.path.split(original_master_data)[0]+'/slcs.list'
#print time_tables
time_tables=work_dir+'/timetableN1.txt'
slc_list   = work_dir+'/'+Master_Date_fold+'/'+'slcs.list'
day_in   = work_dir+'/'+Master_Date_fold+'/'+'day.1.in'
shutil.copy(time_tables,slc_list)
shutil.copy(time_tables,day_in )

#
dem_input_card  = work_dir+'/'+'dem.dorisin'
target_dem_card = work_dir+'/'+Master_Date_fold+'/'+'dem.dorisin'
shutil.copy(dem_input_card,target_dem_card)

python_script_card   = work_dir+'/'+'slave_fileName_replace.py'
target_python_script = work_dir+'/'+Master_Date_fold+'/'+'slave_fileName_replace.py'
shutil.copy(python_script_card,target_python_script)


python_instruction='python slave_fileName_replace.py '
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 

python_instruction='make_coarse'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 

python_instruction='make_coreg_simple'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 


python_instruction='make_dems'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 

python_instruction='make_ifgs'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 
#make_coarse
#make_coreg_simple
#make_dems
#make_ifgs
######################################################
for insar_file_name in Interfergrom_combination_info:
    print '\n'
    print 'insar_file_name=',insar_file_name
    file_time    =os.path.basename(insar_file_name).split('_')[1].strip()
    print 'file_time=',file_time
    if file_time.strip()!=master_date.strip():
        original_fold=insar_file_name.strip()

        print 'original_fold=',original_fold
        target_fold  =work_dir+'/'+Master_Date_fold+'/'+file_time
        print 'work_dir=',work_dir
        print 'file_time=',file_time
        print 'target_fold=',target_fold
        [dirname,filename]=os.path.split(original_fold.strip())
        print(dirname,"\n",filename)
        
        original_source     =dirname+'/'+file_time+'/'+'CInt_Crop'+'/'+'srd_crop.raw' 
        target_source  =target_fold+'/'+'cint.minrefdem.raw'
        shutil.copy(original_source,target_source)
#########################################################################################




geo_path=work_dir+'/'+Master_Date_fold+'/'+Geo_date
os.chdir(geo_path)
python_instruction='step_geo'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 

master_path=work_dir+'/'+Master_Date_fold
os.chdir(master_path)
python_instruction='mt_prep '+Amplitude_dispersion+' '+'4 1 '+' '+' 50 50'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 


#    original_source=original_fold+'/'+'cint.raw'
#    target_source  =target_fold+'/'+'cint.raw'
#    if not (os.path.isfile(target_source)):        
#        os.symlink(original_source,target_source)
#    else:
#        print  target_source+' is already here!!!\n'
#cp $N1datadir/$N1file/burst4/slave.res       $workdir/$Master_Date/$TimeN1/slave.res
#ln -s $N1datadir/$N1file/burst4/${master_date}_iw_2_burst_4.raw       $workdir/$Master_Date/$TimeN1/.
#ln -s $N1datadir/$N1file/burst4/slave_rsmp.raw  $workdir/$Master_Date/$TimeN1/slave_res.slc
#ln -s $N1datadir/$N1file/burst4/slave_rsmp.raw  $workdir/$Master_Date/$TimeN1/.
#ln -s $N1datadir/$N1file/burst4/${TimeN1}_iw_2_burst_*.raw  $workdir/$Master_Date/$TimeN1/.
#echo "${TimeN1}_iw_3_burst_*.raw"
#ln -s $N1datadir/$N1file/burst4/cint.raw   $workdir/$Master_Date/$TimeN1/cint.raw


