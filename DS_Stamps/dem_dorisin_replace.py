#!/usr/bin/env python
import os,sys,time

def usage():
    print '\nUsage: python  dem_dorisin_replace.py  Input_Filename  Output_Filename'
    print '  where                                           '
    print '        Input_Filename     is the dem input parameter file ' 
    print '        Output_Filename    is the dem input parameter file  '
    print '                                                  '
    print ' This python applies cpxint16 transform  to  cpxcomplex64'
    print '                                                          '
    print ' for example                                              '
    print ' python   dem_dorisin_replace.py   dem.dorisin input.dem                '
    print '                                                          '
    print ' Python: Wu Wenhao   Wuhan university   QQ:460249274'
try:
    Input_Filename           = sys.argv[1]
    Output_Filename          = sys.argv[2]
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
if len(sys.argv) == 3:
    print 'Expected input'
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1) 


################################################################################


work_dir=os.getcwd()

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

#f=open('day.1.in')
#cint_minrefdem_file=f.readlines()
#for file_name in cint_minrefdem_file:
    
    #slave_resfile=file_name.strip()+'/slave.res'
    #print slave_resfile


    """
    First_line=get_parameter('First_line (w.r.t. original_master)',slave_resfile,2,'*_Start_resample','* End_resample:_NORMAL') 
    Last_line =get_parameter('Last_line (w.r.t. original_master)',slave_resfile,2,'*_Start_resample','* End_resample:_NORMAL') 
    First_pixel=get_parameter('First_pixel (w.r.t. original_master)',slave_resfile,2,'*_Start_resample','* End_resample:_NORMAL') 
    Last_pixel=get_parameter('Last_pixel (w.r.t. original_master)',slave_resfile,2,'*_Start_resample','* End_resample:_NORMAL')
 
    
    formatData1=get_parameter('Data_output_format',slave_resfile,2,'*_Start_resample','* End_resample:_NORMAL')
    """



if os.path.exists(Input_Filename.strip()):
    #f=open(Input_Filename.strip())
    #lines=f.readlines()
    #f.close()
    inputStream  = open(Input_Filename.strip(),'r')
    textStream   = inputStream.read()
    inputStream.close()

    sourceText   = "I_RESFILE       dem.out                         // parameter file"
    replaceText  = "I_RESFILE       interferogram.out               // parameter file"

    outputStream = open(Output_Filename.strip(), "w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()
    """
    outStream=open(Output_Filename,'w')
    
    for line in lines:
        outStream.write(line)
        if not (line.find('Shifted azimuth spectrum:')):
            break
    output_file='slave_res.slc'
    outStream.write('*******************************************************************\n')
    outStream.write('Data_output_file:                          %s\n' % output_file )
    outStream.write('Data_output_format: 			    %s\n' % formatData1) 
    outStream.write('Interpolation kernel:                 		12 point raised cosine kernel\n')    
    outStream.write('First_line (w.r.t. original_master): 	    %s\n' % First_line)
    outStream.write('Last_line (w.r.t. original_master): 	    %s\n' % Last_line)
    outStream.write('First_pixel (w.r.t. original_master): 	    %s\n' % First_pixel)
    outStream.write('Last_pixel (w.r.t. original_master): 	    %s\n' % Last_pixel)
    outStream.write('*******************************************************************\n')
    outStream.write('* End_resample:_NORMAL\n')
    outStream.write('*******************************************************************\n')

    outStream.write('\n')
    outStream.write('    Current time: %s\n' % time.asctime())
    outStream.write('\n')
    outStream.close()
    """
    """
    # replace crop tag in result file
    sourceText   = "Data_output_file:                     		slave_rsmp.raw"
    replaceText  = "Data_output_file:                     		slave_res.slc"
    #replaceText  =replaceText+ work_dir + "/slave_res.slc"
    print replaceText
    inputStream  = open(slave_resfile,'r')
    textStream   = inputStream.read()
    inputStream.close()
    outputStream = open(slave_resfile, "w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()

    sourceText   = "Data_output_file:                     		slave_rsmp_crop.raw"
    replaceText  = "Data_output_file:                     		slave_res.slc"
    print replaceText
    inputStream  = open(slave_resfile,'r')
    textStream   = inputStream.read()
    inputStream.close()
    outputStream = open(slave_resfile, "w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()
    """

