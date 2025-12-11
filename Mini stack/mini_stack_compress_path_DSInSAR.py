 #!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize
#from scipy.optimize import Bounds
import scipy
print scipy.version.version
import multiprocessing
from scipy.optimize import minimize
import os,sys,time
from scipy.stats import ks_2samp
import math
import numpy as np
from numpy import *
from numpy import dot
from numpy import cos
import numpy
#import matplotlib.pyplot as plt
import gdal,gdalconst
from gdalconst import *
from struct import *
import re
import scipy.linalg
import scipy.io as sio
import linecache
import struct,shutil
import datetime
import copy
import gc
#global Observation_data_array
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
 
#---------------------------------------------------------------#
###############################################################################
def freadbk(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        #print 'Could not open'+Path_MFF_HDR
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
    print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
###############################################################################
###############################################################################
def SAR_write(output_file,Outdata,outformat):
    cols = Outdata.shape[1]
    rows = Outdata.shape[0]
    #print "cols = ",cols
    #print "rows = ", rows
    fid = open(output_file,'wb')
    if fid < 0:
        print "ERROR :cannot open file for writing\n"
        sys.exit(1)
    else: 
        if outformat=='complex_real4':
            Outdata=Outdata.astype(complex64)
            for i_temp in range(0,rows):
                fid.write(Outdata[i_temp,:])
            fid.close()
        elif outformat=='complex_short':
            data_short=zeros((2*cols,1),dtype='int16')
            for i_temp in arange(0,rows):        
                data_short[0::2]= (Outdata[i_temp,:].real.astype(np.int16)).reshape(-1,1,order='F').copy()
                data_short[1::2]= (Outdata[i_temp,:].imag.astype(np.int16)).reshape(-1,1,order='F').copy() 
                fid.write(data_short)
            fid.close() 
        elif outformat=='real4':
            Outdata=Outdata.real.astype(float32)
            for i_temp in range(0,rows):
                fid.write(Outdata[i_temp,:])
            fid.close()
        else:
             print "Sorry! not yet completed !"
             usage()
             sys.exit(1)
#****************************************************************************************
def SAR_data_read(SLC_path_file,data_format,nofLines,nofPixels,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels):
    
    print SLC_path_file
    [fname,fename]=os.path.splitext(SLC_path_file)
    Path_MFF_HDR   =fname+'.hdr'
    if data_format=='complex_real4':
        Link_DATA      =fname+'.x00'
    elif data_format=='complex_short':
        Link_DATA      =fname+'.j00'
    elif data_format=='real4':
        Link_DATA      =fname+'.r00'
    else:
        print "Sorry! not yet completed !"
        usage()
        sys.exit(1)

    print 'Path_MFF_HDR=',Path_MFF_HDR 
    #if os.path.exists(Path_MFF_HDR)==True:
    #   os.remove(Path_MFF_HDR)
    if os.path.exists(Path_MFF_HDR)==False:
       
        outStream      = open(Path_MFF_HDR,'w')
        outStream.write('IMAGE_FILE_FORMAT = MFF\n')
        outStream.write('FILE_TYPE = IMAGE\n')
        outStream.write('IMAGE_LINES = %d\n' % int(nofLines))
        outStream.write('LINE_SAMPLES = %d\n'% int(nofPixels))
        outStream.write('BYTE_ORDER = LSB\n')
        outStream.write('END\n')
        outStream.close()

    if (os.path.isfile(Link_DATA)):
            #os.remove(Link_DATA)
        print  'Link_DATA=' ,Link_DATA
    else:        
        os.symlink(SLC_path_file ,Link_DATA)
            
    Cint_Data=freadbk(Path_MFF_HDR,Patch_First_Line,Patch_First_Pixel,Patch_nofLines,Patch_nofPixels)
    return Cint_Data

#############################################################################################
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
master_date="20190720"
Geo_date="00000000"
swath_num='3'
N1datadir="/home/wuwenhao/disk_external/Hunan/Path_11/Frame_86_91/data_processing/"
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
#if os.path.exists(SAR_Interfer_Data_path):    
#    file_dir=os.listdir(SAR_Interfer_Data_path)
#else:
#    print "This path is wrong !!!"
#    sys.exit(1)



Target_fold_Stamps_compress=work_dir+'/'+'DS_Stamps_compress_4'
if  os.path.exists(Target_fold_Stamps_compress):
    python_instruction='rm -fr '+os.path.abspath(Target_fold_Stamps_compress)
    os.system(python_instruction)
    print "This is the target file fold"
    os.makedirs(Target_fold_Stamps_compress) 
else:
    os.makedirs(Target_fold_Stamps_compress) 

filename_list=[]
large_mini_stack=[]
mini_stack_index_list=[]
Relative_master_mini_stack_list=[]
for filename in open('mini_stack_combine.txt'):
    print 'filename=',filename.split('\t')
    master_date=filename.strip().split('\t')[2].strip()
    filename_list.append(N1datadir+'/'+filename.strip()) 
    large_mini_stack.append('stamps_'+filename.split('\t')[5].strip())
    if (filename.strip()!=master_date.strip()) and\
       (int(filename.split('\t')[1].strip()) >= 0):
            
        mini_stack_index_list.append(int(filename.split('\t')[1].strip()))
        Relative_master_mini_stack_list.append(int(filename.split('\t')[4].strip()))

mini_stack_index_list =list(set(mini_stack_index_list))
mini_stack_index_list=sorted(mini_stack_index_list)
Relative_master_mini_stack_list=list(set(Relative_master_mini_stack_list))
Relative_master_mini_stack_list =sorted(Relative_master_mini_stack_list)
print 'Relative_master_mini_stack_list=',Relative_master_mini_stack_list
#print 'filename_list=',filename_list
#print 'large_mini_stack=',large_mini_stack

large_mini_stack =list(set(large_mini_stack))
large_mini_stack=sorted(large_mini_stack)
print 'large_mini_stack=',large_mini_stack

 
Target_fold_Stamps_compress_InSAR=Target_fold_Stamps_compress+'/'+'INSAR_'+master_date
if  os.path.exists(Target_fold_Stamps_compress_InSAR):
    python_instruction='rm -fr '+os.path.abspath(Target_fold_Stamps_compress_InSAR)
    os.system(python_instruction)
    print "This is the target file fold"
    os.makedirs(Target_fold_Stamps_compress_InSAR)
else:
    os.makedirs(Target_fold_Stamps_compress_InSAR) 
             
geo_initialization_index=0               
#******************************************************************************#
for master_mini_stack in Relative_master_mini_stack_list:                   
    if geo_initialization_index==0 and (int(master_date)!=int(master_mini_stack)):                         
        if abs(int(master_date)-int(Geo_date)) > abs(int(master_date)-int(master_mini_stack)):               
            Geo_date=int(master_mini_stack)
print 'Geo_date=',Geo_date
print 'master_date=',master_date                
geo_initialization_index=1  
#******************************************************************************

Relative_master_mini_stack_list=[]
for ks_stack_num in large_mini_stack:
    large_mini_stack_fold_name = ks_stack_num
    if  os.path.exists(large_mini_stack_fold_name):
        #python_instruction='rm -fr '+os.path.abspath(large_mini_stack_fold_name)
        #os.system(python_instruction)
        print "This is the target file fold"
   
        ori_Master_Date_fold=large_mini_stack_fold_name+'/'+'INSAR_'+master_date
        if  os.path.exists(ori_Master_Date_fold):
       
            geo_initialization_index=0
            Insar_master_check=0
            ks_stack_index=ks_stack_num.split('_')[-1].strip()
##################################################################
            mini_stack_index_list=[]
            Relative_master_mini_stack_list=[]
            for insar_file_name in filename_list:
        
               #print '\n\n\n'
               #insar_file_name=',insar_file_name
               file_time    =insar_file_name.split('\t')[3].strip()        
               master_date=insar_file_name.split('\t')[2].strip()        
               #print 'ks_stack_index=',ks_stack_index               
               if (file_time.strip()!=master_date.strip()) and\
                  (int(insar_file_name.split('\t')[5].strip())==int(ks_stack_index)) and\
                   (int(insar_file_name.split('\t')[1].strip()) >= 0):
            
                   mini_stack_index_list.append(int(insar_file_name.split('\t')[1].strip()))
                   Relative_master_mini_stack_list.append(int(insar_file_name.split('\t')[4].strip()))
            mini_stack_index_list =list(set(mini_stack_index_list))
            mini_stack_index_list=sorted(mini_stack_index_list)
            print 'mini_stack_index_list=',mini_stack_index_list

            Relative_master_mini_stack_list=list(set(Relative_master_mini_stack_list))
            Relative_master_mini_stack_list =sorted(Relative_master_mini_stack_list)
            print 'Relative_master_mini_stack_list=',Relative_master_mini_stack_list
      
           
            for master_mini_stack in Relative_master_mini_stack_list:
                print 'master_mini_stack=',master_mini_stack
                if master_mini_stack !=int(master_date): 
                    original_fold  =work_dir+'/'+ori_Master_Date_fold+'/'+str(master_mini_stack)
                    print 'original_fold=' ,original_fold

                    if os.path.exists(original_fold):
                        
                        target_fold  =Target_fold_Stamps_compress+'/'+'INSAR_'+master_date+'/'+str(master_mini_stack)
                        if  not os.path.exists(target_fold):
                            os.makedirs(target_fold)
                        else:
                            screen_print=target_fold+' is alread here'
                            print screen_print 
                            print '\n'

                        original_data_file  =original_fold+'/'+'slave_res.slc'
                        original_source_res =work_dir+'/'+ori_Master_Date_fold+'/'+str(master_mini_stack)+'/'+'slave.res'
 
                        target_source_res  =target_fold+'/'+'slave.res'
                        target_data_file   =target_fold+'/'+'slave_res.slc'
                        shutil.copy(original_source_res,target_source_res)
                        shutil.copy(original_data_file,target_data_file)


                        original_slave_data=get_parameter('Data_output_file',original_source_res,2,'*_Start_crop','* End_crop:_NORMAL')
                        #print 'original_slave_data=',
                        
                        original_source=work_dir+'/'+ori_Master_Date_fold+'/'+str(master_mini_stack)+'/'+original_slave_data.strip()
                        #print 'original_source=',original_source
                        target_source  =target_fold+'/'+original_slave_data

                        if os.path.isfile(original_source):
                            if not (os.path.isfile(target_source)) :        
                                os.symlink(original_source,target_source)
                                print "link finished!!!\n\n\n"
                            else:
                                print  target_source+' is already here!!!\n'


                        ###################################################################
                        original_source=original_fold+'/'+'merged_h2ph_srd.raw' 
                        target_source  =target_fold+'/'+'merged_h2ph_srd.raw'
        
                        #shutil.copy(original_source,target_source)
                        if not (os.path.isfile(target_source)):        
                            os.symlink(original_source,target_source)
                        else:
                            print  target_source+' is already here!!!\n'
                        ###################################################################

                if master_mini_stack ==int(master_date):
                    original_fold  =work_dir+'/'+ori_Master_Date_fold
                    print 'original_fold=' ,original_fold
                    if os.path.exists(original_fold):
                        
                        target_fold  =Target_fold_Stamps_compress+'/'+'INSAR_'+master_date+'/'
                        if  not os.path.exists(target_fold):
                            os.makedirs(target_fold)
                        else:
                            screen_print=target_fold+' is alread here'
                            print screen_print 
                            print '\n'

                        original_data_file  =original_fold+'/'+master_date+'_crop.slc.cr4'
                        original_source_res =work_dir+'/'+ori_Master_Date_fold+'/'+'master.res'

                        target_data_file    =Target_fold_Stamps_compress_InSAR+'/'+master_date+'_crop.slc'
                        target_data_file_cr4=Target_fold_Stamps_compress_InSAR+'/'+master_date+'_crop.slc.cr4'
                        
                        shutil.copy(original_data_file,target_data_file)
                        shutil.copy(original_data_file,target_data_file_cr4)
                        
                        print 'original_source_res =',original_source_res 
                        if os.path.isfile(original_source_res):
                            print 'original_source_res =',original_source_res
                            f=open(original_source_res)
                            lines=f.readlines()
                            f.close()

                            master_res=os.path.abspath(original_source_res)
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

                            target_source_res =Target_fold_Stamps_compress_InSAR+'/'+'master.res'
                            formatData1='complex_real4'
                            outStream=open(target_source_res,'w')
                            for line in lines:
                                outStream.write(line)
                                if not (line.find('*_Start_crop:')):
                                    break
                            output_file=target_data_file
                            outStream.write('*******************************************************************\n')
                            outStream.write('Data_output_file:                          %s\n' % output_file )
                            outStream.write('Data_output_format: 			    %s\n' % formatData1) 
       
                            outStream.write('First_line (w.r.t. original_image): 	    %s\n' % str(l0))
                            outStream.write('Last_line (w.r.t. original_image): 	    %s\n' % str(lN))
                            outStream.write('First_pixel (w.r.t. original_image): 	    %s\n' % str(p0))
                            outStream.write('Last_pixel (w.r.t. original_image): 	    %s\n' % str(pN))
                            outStream.write('*******************************************************************\n')
                            outStream.write('* End_crop:_NORMAL\n')
                            outStream.write('*******************************************************************\n')

                            outStream.write('\n')
                            outStream.write('    Current time: %s\n' % time.asctime())
                            outStream.write('\n')
                            outStream.close()
                    
                        


   
python_file=open('timetableN1.txt','w')
for insar_file_name in Relative_master_mini_stack_list:       
    print '\n\n'
    
    file_time    =insar_file_name
    print 'file_time=',file_time
          
    if int(file_time)!=int(master_date.strip()):

        ori_source_res =Target_fold_Stamps_compress_InSAR+'/'+'master.res'  
        # original_fold_res=work_dir+'/'+Master_Date_fold+'/'+'master.res'
        target_fold_res  =Target_fold_Stamps_compress_InSAR+'/'+str(file_time)+'/'+'master.res'
        shutil.copy(ori_source_res ,target_fold_res )
        python_file.writelines(str(file_time)+'\n')
python_file.close()



os.chdir(Target_fold_Stamps_compress_InSAR)

     #print os.path.split(original_master_data)
     #time_tables=os.path.split(original_master_data)[0]+'/slcs.list'
     #print time_tables
time_tables=work_dir+'/timetableN1.txt'
slc_list   = Target_fold_Stamps_compress_InSAR+'/'+'slcs.list'
day_in     = Target_fold_Stamps_compress_InSAR+'/'+'day.1.in'
shutil.copy(time_tables,slc_list)
shutil.copy(time_tables,day_in )

     #
dem_input_card  = work_dir+'/'+'dem.dorisin'
target_dem_card = Target_fold_Stamps_compress_InSAR+'/'+'dem.dorisin'
shutil.copy(dem_input_card,target_dem_card)

python_script_card   = work_dir+'/'+'slave_fileName_replace.py'
target_python_script = Target_fold_Stamps_compress_InSAR+'/'+'slave_fileName_replace.py'
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
#############################################################################
for master_mini_stack in Relative_master_mini_stack_list:
    if master_mini_stack !=int(master_date): 
        #if int(master_date) < int(master_mini_stack):
        #original_fold  =work_dir+'/'+ori_Master_Date_fold+'/'+'Relative_master'+'/'+master_date+'_'+str(master_mini_stack)
        #else:
        original_fold  =work_dir+'/'+ori_Master_Date_fold+'/'+'Relative_master'+'/'+str(master_mini_stack)
        print 'original_fold=' ,original_fold

        if os.path.exists(original_fold):
                        
            target_fold  =Target_fold_Stamps_compress+'/'+'INSAR_'+master_date+'/'+str(master_mini_stack)
            if  not os.path.exists(target_fold):
                os.makedirs(target_fold)
            else:
                screen_print=target_fold+' is alread here'
                print screen_print 
                print '\n'

            original_data_file  =original_fold+'/'+'cint.relative_master.raw'
           

            target_data_file   =target_fold+'/'+'cint.minrefdem.raw'
                        
            shutil.copy(original_data_file,target_data_file)


########################################################################### 
    

     
geo_path=Target_fold_Stamps_compress_InSAR+'/'+str(Geo_date)
os.chdir(geo_path)
python_instruction='step_geo'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy 

os.chdir(Target_fold_Stamps_compress_InSAR)
python_instruction='mt_prep '+Amplitude_dispersion+' 2 2 '+' 50 50'
cmd_profile_copy=Linux_cmd(python_instruction)
print 'cmd_profile_copy=',cmd_profile_copy
       
#python_instruction='DS_Nad_Para.py '
#cmd_profile_copy=Linux_cmd(python_instruction)
#print 'cmd_profile_copy=',cmd_profile_copy 
   

os.chdir(work_dir)

