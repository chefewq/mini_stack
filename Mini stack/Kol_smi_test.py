#!/usr/bin/env python
import os,sys,time
import math
import numpy as np
from numpy import *
import re
import scipy as Sci
import scipy.linalg
import scipy.io as sio
import linecache
#import matplotlib.pyplot as plt

import gdal,gdalconst
from gdalconst import *
#Version 20190730 QQ:460249274
#################################################
def usage():
    print '\nUsage: python concatenateIfg.py <burstList> [Degree] [do_ESD_flag] [Plotflag] '
    print '  where burstList  is the list of bursts you want to merge        '
    print '        Degree     is the  order of the polynomial interpolation  '
    print '        do_ESD_flag is a boolean to save the result of azimuth shift by ESD '
    print '        Plotflag    is a boolean to plot, if true, plot the ESD  interferogram  '
    print ' For example (default instructon)                                '
    print '    python ESD_Azimuth.py 1:9  1(default)                                '
    print ' or\n    python ESD.py 1:9  2 True True                           '
    print ' \n                                                              '
    print " Freek's Suggestion                                               "
    print ' Python: Wu Wenhao    QQ:460249274              '

###############################################################################
def freadbk(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print 'Could not open'+Path_MFF_HDR
        sys.exit(1)
    #print 'Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
    #      thisBurstData_file.GetDriver().LongName
    #print 'Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
    #      'x',thisBurstData_file.RasterCount
    #print 'Projection is ',thisBurstData_file.GetProjection()
    geotransform = thisBurstData_file.GetGeoTransform()
    #if not geotransform is None:
    #    print 'Origin = (',geotransform[0], ',',geotransform[3],')'
    #    print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

    cint_srd=thisBurstData_file.GetRasterBand(1)
    #print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
#*******************************************************************************************
#******************************************************************************
def Linux_cmd(python_instruction):
    #python_instruction='source '+work_dir+'/'+Configfile
    print 'python_instruction=',python_instruction
    with os.popen(python_instruction,"r") as p:
        cmd_profile_copy=p.read()
    #print 'cmd_profile_copy=',cmd_profile_copy       
    return cmd_profile_copy

#******************************************************************************
################################################################################
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


        if format_flag==4:
            
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
################################################################################

#######################################################
Amplitude_dispersion="0.50"

Configfile='StaMPS_CONFIG.bash'

##N1datadir="/home/wuwenhao/data_disk/kunming/data_processing"
master_date="20190720"
Geo_date="00000000"
swath_num='3'
N1datadir="/home/wuwenhao/data_disk/kunming/cangyuan_sentinel/data_processing/"
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


for filename in open('mini_stack_combine.txt'):
    print 'filename=',filename.split('\t')
    master_date=filename.strip().split('\t')[2].strip()
    filename_list.append(N1datadir+'/'+filename.strip()) 
    large_mini_stack.append('stamps_'+filename.split('\t')[5].strip())
#print 'filename_list=',filename_list
#print 'large_mini_stack=',large_mini_stack

large_mini_stack =list(set(large_mini_stack))
large_mini_stack=sorted(large_mini_stack)
print 'large_mini_stack=',large_mini_stack
#******************************************************************************#

#******************************************************************************
geo_initialization_index=0

for ks_stack_num in large_mini_stack:
    large_mini_stack_fold_name = ks_stack_num
    if  os.path.exists(large_mini_stack_fold_name):
        #python_instruction='rm -fr '+os.path.abspath(large_mini_stack_fold_name)
        #os.system(python_instruction)
        print "This is the target fold"
    else:
        #os.makedirs(large_mini_stack_fold_name) 
        print "can not find the fold"
 
    Master_Date_fold=large_mini_stack_fold_name+'/'+'INSAR_'+master_date
    if not os.path.exists(Master_Date_fold):
       #os.makedirs(Master_Date_fold)
       print "can not find the fold"

    geo_initialization_index=0
    Insar_master_check=0
    ks_stack_index=ks_stack_num.split('_')[-1].strip()
##################################################################
    
    mini_stack_index_list=[]
    for insar_file_name in filename_list:
        
        #print '\n\n\n'
        #insar_file_name=',insar_file_name
        file_time    =insar_file_name.split('\t')[3].strip()
        
        master_date=insar_file_name.split('\t')[2].strip()
        
        #print 'ks_stack_index=',ks_stack_index
        if int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):
                      
            #int(insar_file_name.split('\t')[-1].strip())==int(ks_stack_index):
            if int(insar_file_name.split('\t')[1].strip()) >= 0:
                mini_stack_index_list.append(insar_file_name.split('\t')[1].strip())
    mini_stack_index_list =list(set(mini_stack_index_list))
    mini_stack_index_list=sorted(mini_stack_index_list)
    #print 'mini_stack_index_list=',mini_stack_index_list  

    Relative_master_mini_stack_list=[]
    mini_stack_combine_list=[]
    for mini_stack_index_num in mini_stack_index_list:
        mini_stack_file_name_list=[]
        
        for insar_file_name in filename_list:
              
            file_time    =insar_file_name.split('\t')[3].strip()
        
            master_date=insar_file_name.split('\t')[2].strip()
           
            #print ' mini_stack_index_num=', mini_stack_index_num
            if int(insar_file_name.split('\t')[1].strip())==int(mini_stack_index_num)\
               and int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index):
                print int(insar_file_name.split('\t')[1].strip())
                print insar_file_name.split('\t')[4].strip()
                print insar_file_name
                print '\n'
                mini_stack_file_name_list.append(file_time) 
                Relative_master_mini_stack_list.append(insar_file_name.split('\t')[4].strip()) 
        #print '\n\n\n'
        #print ' mini_stack_index_num=', mini_stack_index_num          
        #print 'mini_stack_file_name_list=',mini_stack_file_name_list
    Relative_master_mini_stack_list=list(set(Relative_master_mini_stack_list))
        
    #print 'Relative_master_mini_stack_list=',Relative_master_mini_stack_list
    abs_path_Master_Date_fold=os.path.abspath(Master_Date_fold)

    
    
    os.chdir(Master_Date_fold)
    master_res=os.path.abspath('master.res')
    if get_parameter('First_line (w.r.t. ovs_image)',master_res,1): #oversampled data
        l0 = int(get_parameter('First_line (w.r.t. ovs_image)',master_res,1))  
        lN = int(get_parameter('Last_line (w.r.t. ovs_image)',master_res,1))  
        p0 = int(get_parameter('First_pixel (w.r.t. ovs_image)',master_res,1))
        pN = int(get_parameter('Last_pixel (w.r.t. ovs_image)',master_res,1))
    

    else: #original data  
        l0 = int(get_parameter('First_line (w.r.t. original_image)',master_res,1))
        lN = int(get_parameter('Last_line (w.r.t. original_image)',master_res,1))
        p0 = int(get_parameter('First_pixel (w.r.t. original_image)',master_res,1))
        pN = int(get_parameter('Last_pixel (w.r.t. original_image)',master_res,1))

    Image_width=pN-p0+1

    
    python_instruction='mt_prep_mini_stack '+Amplitude_dispersion+" "+"3 2 50 50"
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy

   
    #calamp calamp.in $width $WORKDIR/calamp.out

    master_cr4_path=abs_path_Master_Date_fold+'/'+master_date+'_crop.slc.cr4'
    output_calamp_file= abs_path_Master_Date_fold+'/calamp.in'          
    python_file=open(output_calamp_file,'w')
    python_file.writelines(master_cr4_path +'\n')
    for i_index_name in  Relative_master_mini_stack_list:
        if int(i_index_name.strip())!=int(master_date.strip()):
            python_file.writelines(abs_path_Master_Date_fold+'/'+i_index_name.strip()+'/slave_res.slc' +'\n')
    python_file.close()
    

    python_instruction='calamp calamp.in '+str(Image_width)+' calamp.out'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy

    python_file=open('selpsc.in','w')
    python_file.writelines(Amplitude_dispersion+'\n')
    python_file.writelines(str(Image_width) +'\n')
    for i_index_name in  open('calamp.out'):
        python_file.writelines(i_index_name.strip() +'\n')
    python_file.close()


    python_instruction='mt_extract_cands_mini_stack'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy
     
    python_instruction='DS_Nad_Para.py'
    cmd_profile_copy=Linux_cmd(python_instruction)
    print 'cmd_profile_copy=',cmd_profile_copy
    
    os.chdir(work_dir)
    
