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
master_date="20190420"
Geo_date="00000000"
swath_num='2'
#N1datadir="/home/wuwenhao/disk_external/kunming/data_processing/"
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
#    print "This path Int_Data_merge.listis wrong !!!"
#    sys.exit(1)
filename_list=[]
large_mini_stack=[]
original_large_mini_stack=[]
Interfergrom_combination_info=[]
for filename in open('Int_Data_merge.list'):
    #print 'filename=',filename.strip()
    master_date=os.path.basename(filename.strip()).split('_')[0].strip()
    slave_date=os.path.basename(filename.strip()).split('_')[1].strip()
    filename_list.append(slave_date) 
    Interfergrom_combination_info.append(filename.strip())    
filename_list.append(master_date)
#print 'filename_list=',filename_list
#print 'large_mini_stack=',large_mini_stack
Interfergrom_combination_info =list(set(Interfergrom_combination_info))
Interfergrom_combination_info =sorted(Interfergrom_combination_info)


Geo_date=filename_list[0].strip()
for temp_file in filename_list:    
    if abs(int(master_date)-int(temp_file)) > abs(int(master_date)-int(Geo_date)):
        Geo_date=temp_file.strip()
print 'Geo_date=',Geo_date 


for mini_stack_filename in open('mini_stack_combine.txt'):
    original_large_mini_stack.append('stamps_'+mini_stack_filename.split('\t')[5].strip())

large_mini_stack =list(set(original_large_mini_stack))
large_mini_stack=sorted(large_mini_stack)
print 'large_mini_stack=',large_mini_stack
#******************************************************************************#
file_mini_stack= open('mini_stack_combine.txt')
mini_info=file_mini_stack.readlines()
file_mini_stack.close()

#******************************************************************************
geo_initialization_index=0

for ks_stack_num in large_mini_stack:
    large_mini_stack_fold_name = ks_stack_num
    if  os.path.exists(large_mini_stack_fold_name):
        python_instruction='rm -fr '+os.path.abspath(large_mini_stack_fold_name)
        os.system(python_instruction)
    if not os.path.exists(large_mini_stack_fold_name):
        os.makedirs(large_mini_stack_fold_name) 
 
    Master_Date_fold=large_mini_stack_fold_name+'/'+'INSAR_'+master_date
    if not os.path.exists(Master_Date_fold):
       os.makedirs(Master_Date_fold)

    geo_initialization_index=0
    Insar_master_check=0
    ks_stack_index=ks_stack_num.split('_')[-1].strip()
##################################################################
    temp_index=0
    #print  mini_info
    for insar_file_name in mini_info:
        
        #print '\n\n\n'
        print 'insar_file_name=',insar_file_name
        #time.sleep(1000000000)
        
        file_time    =insar_file_name.split('\t')[3].strip()       
        master_date  =insar_file_name.split('\t')[2].strip()

        
        if (file_time.strip()!=master_date.strip()) and\
            int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):
            temp_index=temp_index+1
            print 'temp_index=',temp_index
            #print "int(insar_file_name.split('\t')[5].strip())=",int(insar_file_name.split('\t')[5].strip())
            print "Interfergrom_combination_info[temp_index-1]=",Interfergrom_combination_info[temp_index-1]
            [dirname,filename]  =os.path.split(Interfergrom_combination_info[temp_index-1].strip())
            print(dirname,"\n",filename) 
            temp_salve_file_time =os.path.basename(Interfergrom_combination_info[temp_index-1]).split('_')[1].strip()
            original_fold             =dirname+'/'+temp_salve_file_time+'/'+'Crop'
            original_master_fold      =dirname+'/'+master_date.strip()+'/'+'Crop'
            original_slave_Burst_fold =dirname+'/'+temp_salve_file_time+'/'+'Burst'
            print 'original_fold=',original_fold
            target_fold  =work_dir+'/'+Master_Date_fold+'/'+temp_salve_file_time
    
            if not os.path.exists(target_fold):
                os.makedirs(target_fold)
            else:
                screen_print=target_fold+' is alread here'
                print screen_print 
                print '\n'
            

            [dirname,filename]=os.path.split(original_fold.strip())
            print(dirname,"\n",filename)
            swath_num=dirname.split('iw')[1][0].strip()
            #inner_original_fold='/'+filename+'_iw'+swath_num+'_fine'
            inner_original_fold=''
            

            original_source_res=original_fold+'/'+inner_original_fold+'/'+'resample_crop.res' 
            target_source_res  =target_fold+'/'+'slave.res'
            shutil.copy(original_source_res,target_source_res)

            slave_rsmp_raw=get_parameter('Data_output_file',original_source_res,2,'*_Start_resample','* End_resample:_NORMAL') 

            original_source=original_fold+'/'+inner_original_fold+'/'+'resample_crop.raw' 
            target_source  =slave_rsmp_raw
            target_source_alias  =target_fold+'/'+'slave_res.slc'
            #shutil.copy(original_source,target_source)
            print 'original_source=',original_source
            print 'target_source=',target_source
            print 'target_source_alias  =',target_source_alias 
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
          

           
            ########################################################################
            original_slave_data=get_parameter('Data_output_file',original_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
            print 'original_slave_data=', original_slave_data
            master_original_source_res=original_fold+'/'+inner_original_fold+'/'+master_date+'_iw_'+swath_num+'_crop.res'
            print 'original_fold=',original_fold 
            original_source=original_slave_Burst_fold+'/'+original_slave_data.strip()
            print 'original_source=',original_source
            target_source  =target_fold+'/'+original_slave_data
            
            if os.path.isfile(original_source):
                if not (os.path.isfile(target_source)) :        
                    os.symlink(original_source,target_source)
                    print "link finished!!!\n\n\n"
                else:
                    print  target_source+' is already here!!!\n'
            #time.sleep(100000000000)

            original_source_res=original_master_fold+'/'+inner_original_fold+'/'+master_date+'_iw'+swath_num+'_crop.res'
            print 'original_source_res=',original_source_res
            target_source_res  =target_fold+'/'+'master.res'
            if not (os.path.isfile(target_source_res)):        
                shutil.copy(original_source_res,target_source_res)
                print '\n'
            else:
                print  target_source+' is already here!!!\n'  
 

            original_source_res=original_master_fold+'/'+inner_original_fold+'/'+master_date+'_iw'+swath_num+'_crop.res'
            target_source_res  =work_dir+'/'+Master_Date_fold+'/'+'master.res'
            if not (os.path.isfile(target_source_res)):        
                shutil.copy(original_source_res,target_source_res)
                print '\n'
            else:
                print  target_source+' is already here!!!\n'

    
    target_source_res  =work_dir+'/'+Master_Date_fold+'/'+'master.res'
    sourceText   = get_parameter('Data_output_file',target_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
    replaceText  = work_dir+'/'+Master_Date_fold+'/'+master_date+'_crop.slc'
    #print 'sourceText   = ',sourceText
    #print 'replaceText  =',replaceText
    inputStream  = open(target_source_res,'r')
    textStream   = inputStream.read()
    inputStream.close()
    outputStream = open(target_source_res,"w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()


    Insar_master_check=0
    if Insar_master_check==0:    
        original_master_data=get_parameter('Data_output_file',target_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
        print 'original_master_data=',original_master_data
        original_source     =original_master_fold+'/'+inner_original_fold+'/'+master_date+'_iw'+swath_num+'_crop.raw'
        
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
    python_instruction='master_int2float.py '+target_source_res
    print 'python_instruction=',python_instruction
    cmd_profile_copy=Linux_cmd(python_instruction) 
    print 'cmd_profile_copy=',cmd_profile_copy 

    
    python_file=open('timetableN1.txt','w')
    for insar_file_name in mini_info:
        
        print '\n\n\n'
        #insar_file_name=',insar_file_name
        file_time    =insar_file_name.split('\t')[3].strip()
        print 'file_time=',file_time
        master_date=insar_file_name.split('\t')[2].strip()
        
        if file_time.strip()!=master_date.strip() and\
        insar_file_name.split('\t')[5].strip()==ks_stack_index:
            
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
    """
    python_script_card   = work_dir+'/'+'slave_fileName_replace.py'
    target_python_script = work_dir+'/'+Master_Date_fold+'/'+'slave_fileName_replace.py'
    shutil.copy(python_script_card,target_python_script)
    """

    python_instruction='slave_fileName_replace.py '
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
    
    ######################################################################## 
    temp_index=0    
    for insar_file_name in mini_info:
        
        print '\n\n\n'
        #insar_file_name=',insar_file_name
        file_time    =insar_file_name.split('\t')[3].strip()
        
        master_date=insar_file_name.split('\t')[2].strip()
        
        print 'ks_stack_index=',ks_stack_index
        if (file_time.strip()!=master_date.strip()) and\
            int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):
            temp_index=temp_index+1

            print "Interfergrom_combination_info[temp_index-1]=",Interfergrom_combination_info[temp_index-1]
            [dirname,filename]  =os.path.split(Interfergrom_combination_info[temp_index-1].strip())
            print(dirname,"\n",filename) 
            temp_salve_file_time =os.path.basename(Interfergrom_combination_info[temp_index-1]).split('_')[1].strip()
            original_CInt_fold             =dirname+'/'+temp_salve_file_time+'/'+'CInt_Crop'
           
            
            target_fold  =work_dir+'/'+Master_Date_fold+'/'+file_time
    
            if not os.path.exists(target_fold):
                os.makedirs(target_fold)
            else:
                screen_print=target_fold+' is alread here'
                print screen_print 
                print '\n'

            original_source=original_CInt_fold+'/'+'srd_crop.raw' 
            target_source  =target_fold+'/'+'cint.minrefdem.raw'
        
            #shutil.copy(original_source,target_source)
            if  (os.path.isfile(original_source)):        
                shutil.copy(original_source,target_source)
            else:
                print  target_source+' is already here!!!\n'
            #*******************************************************************************
           
           


    #time.sleep(100000000000000000) 
    geo_path=work_dir+'/'+Master_Date_fold+'/'+Geo_date
    os.chdir(geo_path)
    python_instruction='step_geo'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy 

    master_path=work_dir+'/'+Master_Date_fold
    os.chdir(master_path)
    """
    python_instruction='mt_prep_mini_stack '+Amplitude_dispersion+' '+'4'+' '+'8'+' '+' 50 50'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy
    """ 
    """   
    python_instruction='DS_Nad_Para.py '
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy 
    """

    os.chdir(work_dir)

