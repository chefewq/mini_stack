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
from datetime import date
import re
import shutil
import copy
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
def Caltime(date1,date2):
    date1=time.strptime(date1,"%Y%m%d")
    date2=time.strptime(date2,"%Y%m%d")

    date1=datetime.datetime(date1[0],date1[1],date1[2])
    date2=datetime.datetime(date2[0],date2[1],date2[2]) 
    return date2-date1
################################################################
def convertstringtodate(stringtime):
 
    if stringtime[0:2] == "20":
        time_strptime=time.strptime(stringtime,"%Y%m%d")
       
        begintime=datetime.datetime(time_strptime[0],time_strptime[1],time_strptime[2])
        return begintime
    else :
        year="20"+stringtime[0:2]
        month=stringtime[2:4]
        day=stringtime[4:6]
        begintime=date(int(year),int(month),int(day))
    return begintime

def comparetime(nowtime,stringtime):

    if isinstance(nowtime,date):
        pass
    else:
        nowtime=convertstringtodate(nowtime)
    if isinstance(stringtime,date):
        pass
    else:
        stringtime=convertstringtodate(stringtime)
    result=nowtime-stringtime
    return result.days
###############################################################
def is_date(str):
    try:
       time.strptime(str,"%Y%m%d")
       return True
    except:
       return False
################################################################################################
def main():
     
    #cores=multiprocessing.cpu_count()
    #print 'Total cores is ', cores
    #pool=multiprocessing.Pool(processes=(2))
    Amplitude_dispersion="0.49"
#######################################################

    Configfile='StaMPS_CONFIG.bash'

    N1datadir="/home/wuwenhao/data_disk/kunming/data_processing"
    master_date="20190329"
    Geo_date="20190410"
    swath_num='3'
    Mini_stack_size=10
    mini_stack_time_baseline=120
    ks_test_size=10009
    work_dir=os.getcwd();
    print 'work_dir=' ,work_dir


    filename_list=[]
    time_list=[]
    slave_date_list=[]
    Relative_master_list=[]
    Relative_master_list_prepose=[]
    Relative_master_list_post=[]
    mini_stack_index_list=[]
    ks_test_index_list=[]
    for filename in open('Int_Data_merge.list'):
        #print 'filename=',filename.strip()         
        #file_name, file_extension =os.path.splitext(filename.strip())
        #print 'file_name=',file_name
        #print 'file_extension=',file_extension
        file_basename=os.path.basename(filename.strip())
        master_date=file_basename.strip().split('_')[0]
        slave_date =file_basename.strip().split('_')[1]
        #print 'master_date=',master_date
        #print 'slave_date=',slave_date
        time_list.append(slave_date.strip())
        slave_date_list.append(slave_date.strip())
        filename_list.append(N1datadir+'/'+filename.strip()) 
    time_list.append(master_date.strip())
    #print 'filename_list=',filename_list

    sorted_nums = sorted(enumerate(time_list), key=lambda x: x[1])
    idx = [i[0] for i in sorted_nums]
    nums = [i[1] for i in sorted_nums]
    Image_time_list= [time_list[(i[0])] for i in sorted_nums]
    print 'Image_time_list=',Image_time_list
    Image_stamp_time_list=[]
    for temp_image_time in Image_time_list:
        date_stamp=time.strptime(temp_image_time.strip(),"%Y%m%d")
        #print 'date_stamp=',date_stamp
        Image_stamp_time_list.append(date_stamp)
    print Image_stamp_time_list
  
    print '\n\n\n\n'
    Relative_master_list=copy.deepcopy(Image_time_list)
    Relative_master_list_prepose=copy.deepcopy(Image_time_list)
    Relative_master_list_post=copy.deepcopy(Image_time_list)
    mini_stack_index_list=copy.deepcopy(Image_time_list)
    ks_test_index_list=copy.deepcopy(Image_time_list)

    sorted_nums = sorted(enumerate(slave_date_list), key=lambda x: x[1])
    idx = [i[0] for i in sorted_nums]
    nums = [i[1] for i in sorted_nums]
    slave_time_list= [slave_date_list[(i[0])] for i in sorted_nums]
    print 'slave_time_list=',slave_time_list
    print '\n\n\n\n'
    #print 'slave_time_list=',slave_time_list[-1]
    #print 'slave_time_list=',slave_time_list[0]
    Time_begin_stamps = time.strptime(slave_time_list[0],"%Y%m%d")
    Time_begin_time   = datetime.datetime(Time_begin_stamps[0],Time_begin_stamps[1],Time_begin_stamps[2])
  
    if is_date(slave_time_list[0]) and is_date(slave_time_list[-1]):
        print 'slave_time_list=',slave_time_list[-1]
        Slave_time_baseline=Caltime(slave_time_list[0],slave_time_list[-1])
        print 'Slave_time_Baseline=',Slave_time_baseline.days
    #***************************************************************************#
        if Slave_time_baseline.days%mini_stack_time_baseline > 0:
            Total_Numer_Mini_Stack=(Slave_time_baseline.days/mini_stack_time_baseline)+1
        else:
            Total_Numer_Mini_Stack=(Slave_time_baseline.days/mini_stack_time_baseline) 
        Mini_Stack_Real_Index=0
        for Mini_Stack_Index in range(0,Total_Numer_Mini_Stack):

            cint_minrefdem_file=[]
            mini_stack_file=[]

            Mini_stack_start=Mini_Stack_Index*mini_stack_time_baseline
            Mini_stack_end  =(Mini_Stack_Index+1)*mini_stack_time_baseline
            
            #datetime.datetime
            
            Mini_stack_start_stamp=Time_begin_time+datetime.timedelta(days=Mini_stack_start)
            Mini_stack_end_stamp  =Time_begin_time+datetime.timedelta(days=Mini_stack_end)
            print 'Mini_stack_start=',Mini_stack_start
            print 'Mini_stack_end=',Mini_stack_end

            print 'Mini_stack_start_stamp=',Mini_stack_start_stamp
            print 'Mini_stack_end_stamp=',Mini_stack_end_stamp
            print 'Mini_stack_end_stamp=',type(Mini_stack_end_stamp)

            if Mini_stack_end > Slave_time_baseline.days:
                Mini_stack_end = Slave_time_baseline.days
 
            
            for temp_image_time in Image_time_list:
                
                print 'temp_image_time=',convertstringtodate(temp_image_time)

                if convertstringtodate(temp_image_time) < Mini_stack_end_stamp \
                   and convertstringtodate(temp_image_time) >= Mini_stack_start_stamp:
                    mini_stack_file.append(temp_image_time.strip())
            print 'mini_stack_file=',mini_stack_file

            for i_index_mini_stack in mini_stack_file:     

                for j_index_mini_stack in mini_stack_file:
                    if int(j_index_mini_stack) > int(i_index_mini_stack):
                        relative_master_time=i_index_mini_stack.strip()
                        relative_slave_time =j_index_mini_stack.strip() 
                        #print master_time.strip(), slave_time.strip()
                        temp_file=relative_master_time.strip()+'_'+relative_slave_time.strip()
                        cint_minrefdem_file.append(temp_file)
            print 'cint_minrefdem_file=',cint_minrefdem_file
            #mini_stack_file=Image_time_list[Mini_stack_start:Mini_stack_end]
            #print 'mini_stack_file=',mini_stack_file
         
            Relative_MASTER_DATE_INDEX=int((len(mini_stack_file)-1)/2)
            print 'Relative_MASTER_DATE_INDEX=',Relative_MASTER_DATE_INDEX

            Relative_master_mini_stack=mini_stack_file[Relative_MASTER_DATE_INDEX]
            print 'Relative_master_mini_stack=',Relative_master_mini_stack
            print '\n\n\n\n'

        
            
            
            for i_index_mini_stack in mini_stack_file:
                temp_index_num=Image_time_list.index(i_index_mini_stack)
                if len(mini_stack_file) <= 2:
                    mini_stack_index_list[temp_index_num]=-1
                else:
                    mini_stack_index_list[temp_index_num]=Mini_Stack_Real_Index
                    
                if master_date.strip() in mini_stack_file: 
                    Relative_master_list[temp_index_num]        =master_date.strip()
                    Relative_master_list_prepose[temp_index_num]=master_date.strip()
                    Relative_master_list_post[temp_index_num]   =master_date.strip()
                else:
                    Relative_master_list[temp_index_num]           =Relative_master_mini_stack
                    if len(mini_stack_file) > 2:
                        Relative_master_mini_stack_prepose=mini_stack_file[Relative_MASTER_DATE_INDEX+1]
                        print 'Relative_master_mini_stack_prepose=',Relative_master_mini_stack_prepose
                        #print '\n\n\n\n'
                        Relative_master_list_prepose[temp_index_num]   =Relative_master_mini_stack_prepose
                        Relative_master_mini_stack_post=mini_stack_file[Relative_MASTER_DATE_INDEX-1]
                        print 'Relative_master_mini_stack_post=',Relative_master_mini_stack_post
                        #print '\n\n'
                        Relative_master_list_post[temp_index_num]      =Relative_master_mini_stack_post
                    else:
                        Relative_master_list_prepose[temp_index_num]=-1
                        Relative_master_list_post[temp_index_num]   =-1
            
                ks_test_index_list[temp_index_num]=int(temp_index_num/ks_test_size) 
                print 'Mini_Stack_Index=',Mini_Stack_Index
                print 'ks_test_index_list[temp_index_num]=',ks_test_index_list[temp_index_num]
   
            if len(mini_stack_file) > 2:
                Mini_Stack_Real_Index= Mini_Stack_Real_Index+1 

    print 'Relative_master_list=',Relative_master_list
    print '\n\n\n\n'
    #Relative_master=

    temp=0
    fileinfo = open('mini_stack_combine.txt','w')
    for i in range(len(Image_time_list)):
        temp=temp+1
        if (str(Image_time_list[i])!='NULL'):
        
            fileinfo.write('%d	%d	%s	%s	%s	%d	%s	%s\n'\
            % (temp,mini_stack_index_list[i],master_date, (Image_time_list[i]),\
              Relative_master_list[i],ks_test_index_list[i],\
              Relative_master_list_prepose[i], Relative_master_list_post[i]))
    
    fileinfo.close()
    
    
if __name__== '__main__':
    main()
