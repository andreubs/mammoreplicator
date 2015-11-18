#!/usr/bin/python
 
#
#                 ** Reproducing 2D mammograms with 3D printed surrogates **
#                                  [Andreu Badal, 2015-07-10, 2015-02-03]
#
#   Process a DICOM mammography to generate a printable triangle mesh
#
#

    # Code adapted for faster execution with pypy just-in-time compilation    !!PYPY!!  [AB, 2015-11-16]
    # Importing "numpypy" instead of "numpypy". The "dicom" package setup file has to be recompiled with pypy

    # Water filled phantom: combine plastic and water to speed up print and reduce artifacts.  !!WaterFill!!  [AB, 2015-11-16]

#NOTE:To disable PYPY execution and use the standard python and numpy, comment next line and uncomment following: 
import numpypy; from numpypy import zeros, array
#import numpy as numpypy; from numpy import zeros, array     # Not using pypy

import time
import math
from math import log, sqrt, exp
import sys

#   Instal DICOM reader for Python: http://www.pydicom.org/ (just run ez_setup.py as root). Compiled with PYPY
import dicom

# Simple function to convert a float into a string with 4 decimal places (problem: 0.0 becomes 0.0000):
#def str4(n):
  #return "{0:.4f}".format(float(n))


# Input data:
print '\n** "mammoreplicator": Replicating 2D mammograms with 3D printed surrogates **'
print '                                               [Andreu Badal, 2015-07-10]\n'

print ' Reading input data...'

image_file_name = "../CancerImaging_TCGA-BRCA_TCGA-AO-A0JI_Raw_0.dcm"
compression_thickness    = 6.0    # cm
average_breast_intensity = 1400
air_intensity            = 16383     # Pixel value outside the breast (air attenuation only)      Matthew Clark 

output_binning = 5               # Rebin the pixels to reduce the print resolution (1=no-binning, 2=average_4pixels...)
subtract_layer = 0.0 # cm
output_dcm = 0     # 1==true, 0==false->do not output dicom file with thickness
output_mm  = 1     # 1==true, 0==false->output triangles in cm
output_base= 1     # 1==true, 0==false->add a flat support at the base of the object to define a solid object, or output just an empty 2D shell

rows       = [0,2294] 
columns    = [0, 800]    # [0,1914]  

pixel_size      = 0.0100  # cm  Real pixel size: 0.00941 cm     # Desired pixel size in the printed phantom
mfp_breast      = 1.284   # cm     # 16.5keV=1.284cm, 20keV=1.933cm;      // Molybdenum target, 28 kVp -> average energy 16.5 keV
mfp_plastic     = 1.429   # cm  = 1/0.7 -> Approx ABS MFP at 16.5keV
  #mfp_plastic = 1.0  # cm -> Approx SLA polymer MFP
  #mfp_plastic = (1.0/21.47) # cm  <-- Aluminum, 15 keV; attenuation coeff. http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html (density=2.699g/cm3)
  #mfp_plastic = 0.863  # cm     # 16.5keV=0.863cm, 20keV=1.315cm
  
source_height_Z = 66.0   # cm


# -- Input variables for the water fill:        !!WaterFill!!  
#    If "min_height>50" water fill is completely disabled.
min_height = 4.5        # cm
height_water_fill = 4.3 # cm
border_cm  = 0.4        # cm
mfp_water  = 0.7675     # cm at 16.5 keV      
  

num_rows   = int((rows[1]-rows[0])/output_binning)*output_binning         # Processed sizwe will be a multiple of binning value 
num_columns= int((columns[1]-columns[0])/output_binning)*output_binning
rows[1]    = rows[0]+num_rows
columns[1] = columns[0]+num_columns

center_pixel    = numpypy.array([0,0])   #([num_columns/2,num_rows])      # Location x-ray field of view center in pixel units [x,y]=[Row,Column]
center_coord    = numpypy.array([0.0, 0.0, 0.0])  # Location x-ray field of view center in cm (printed phantom origin)
source_coord    = numpypy.array([0.0, (rows[0]+(rows[1]-rows[0])*0.5)*pixel_size, source_height_Z])  # Location x-ray source focal spot. Using NUMPY arrays for calculations.

air_threshold   = 0   #!!DeBuG!! Currently not removing air pixels: drawing 0 height columns instead       # Any pixel below this threshold will be considered air and asigned 0 thickness    


print "\n -- Input conversion parameters ("+time.strftime("%c")+"):"
print "      compression_thickness    =", compression_thickness
print "      average_breast_intensity =", average_breast_intensity
print "      air_intensity            =", air_intensity
print "      columns=["+str(columns[0])+","+str(columns[1])+"], num_columns = ", num_columns
print "         rows=["+str(rows[0])+","+str(rows[1])+"], num_rows = ", num_rows
#print "      max_pixel_value =", max_pixel_value
#print "      air_threshold   =", air_threshold
print "      mfp_breast      =", mfp_breast
print "      mfp_plastic     =", mfp_plastic
print "      center_pixel    =", center_pixel
print "      center_coord    =", center_coord
print "      source_coord    =", source_coord
print "      pixel_size      =", pixel_size
print "      subtract_layer  =", subtract_layer
print "      output_binning  =", output_binning
print "      output_mm       =", output_mm
print "      output_base     =", output_base
print "      binned columns  = "+str(num_columns/output_binning)+", binned rows = "+str(num_rows/output_binning)
print " "

# Read image data:
print ' -- Reading input file: ',image_file_name,'...'
mammo = dicom.read_file(image_file_name) 


# Process pixel values:
print '\n Processing pixel values; rows = ',num_rows,', columns = ',num_columns,'...'
sys.stdout.write(' Current row: ');               # Used instead of print bc no new line added after text
pix = numpypy.zeros((rows[1], columns[1]))        # We will work on this copy of the data as floats, init to 0
num_air_pix = 0

conversion_factor = (compression_thickness/log(float(average_breast_intensity/float(air_intensity)))) * (mfp_plastic/mfp_breast) 
# Test1:    conversion_factor = (compression_thickness*log(float(average_breast_intensity))) * (mfp_plastic/mfp_breast)           #   Matthew Clark operation change 7/15/2015
# Original: conversion_factor = (compression_thickness/log(float(average_breast_intensity))) * (mfp_plastic/mfp_breast) 

for j in range(rows[0], rows[1], output_binning):
  sys.stdout.write(str(j)+', '); sys.stdout.flush()
  for i in range(columns[0], columns[1], output_binning):
  
    # If requested, re-bin the output to reduce the print resolution (1=no-binning, 2=average_4pixels...):
    if output_binning>1:
      avg = 0.0
      for k1 in range(0, output_binning):
        for k2 in range(0, output_binning):
            #print 'j='+str(j)+' i='+str(i)+':  j+k1='+str(j+k1)+' i+k2='+str(i+k2)    # !!Verbose!!
          avg = avg + float(mammo.pixel_array[j+k1][i+k2])
      mammo.pixel_array[j][i] = int(round(avg/(output_binning**2)))
  
    # Skip air pixels (leave them as 0 thickness):
    if mammo.pixel_array[j][i]<air_threshold:
      num_air_pix = num_air_pix + 1
      pix[j][i] = 0
    else:  

      if (mammo.pixel_array[j][i]>1):

        # Conversion for raw mammograms "for processing":
        pix[j][i] = conversion_factor * log(float(mammo.pixel_array[j][i])/float(air_intensity))    # Inverted operation: assuming air value = incident intensity (Matt Clark)
        
        # Conversion for post-processed mammograms "for presentation": 
        #pix[j][i] = log(float(mammo.pixel_array[j][i])) * conversion_factor  
        
      else:
        pix[j][i] = 0

    if (output_dcm==1):  
      bin_value = int(round(pix[j][i]*100))        # Convert cm floats into tenths of mm integers to store in dicom image
      for k1 in range(0, output_binning):
        for k2 in range(0, output_binning):      
          mammo.pixel_array[j+k1][i+k2] = bin_value
          
sys.stdout.write('\n'); sys.stdout.flush()
  
if (output_dcm==1):  
  # Save final data as dicom:      help(dicom.write_file)  
  print '\n\n Saving final data as DICOM image in integer units of tenths of mm of plastic...'
  mammo.PixelData = mammo.pixel_array.tostring()   # Save the data array into the actual pixel data string    !!DeBuG!! "tostring" is not implemented yet in numpypy!! Dicom output will only work with regular python
  
  dicom.write_file(image_file_name+"_bin"+str(output_binning)+".dcm", mammo)


# Generate triangle mesh:
print '\n Generating the output triangle mesh focused to the x-ray source focal spot (num_air_pix='+str(num_air_pix)+')...'
num_binned_pixels = (num_columns/output_binning) * (num_rows/output_binning) - num_air_pix
num_vertices  = num_binned_pixels * 8     #  8 vertices per printed column  

if (output_base==1):
  num_triangles = num_binned_pixels * 8 + 2*(num_rows/output_binning) + 2*(num_columns/output_binning)        #!!ReducedTriangles!!
else:    
  num_triangles = num_binned_pixels * 6 + 2*(num_rows/output_binning) + 2*(num_columns/output_binning)

  
image_file_name2 = image_file_name+"_RawData_bin"+str(output_binning)+"_focusedZ"+str(source_coord[2])
if (subtract_layer>0.00001):
  image_file_name2 = image_file_name2+"_-"+str(subtract_layer)+"cm"
if (output_mm==1):
  image_file_name2 = image_file_name2+"_mm"

if (min_height<50.0):
  image_file_name2 = image_file_name2+"_WaterFill"+str(height_water_fill)+"cm"      # !!WaterFill!!
  
image_file_name2 = image_file_name2+".ply"

print '\n -- Writing ',num_triangles,' triangles to output file: '+image_file_name2+'\n'

ply = open(image_file_name2, 'w')  # Create PLY file and write header:

ply.write('ply\n')
ply.write('format ascii 1.0\n')
ply.write('comment ** Geometry created by mammoreplicator.py using the input file: '+image_file_name+'\n')
ply.write('comment      [Andreu Badal, 2015-02-03]\n')
ply.write('comment \n')
if (subtract_layer>0.00001):
  ply.write('comment NOTE: subtracting '+str(subtract_layer)+' cm from each column for easier printing.\n')
  ply.write('comment \n')
if (output_mm==1):
  ply.write('comment NOTE: vertex coordinates in mm.\n')
  ply.write('comment \n')
  
ply.write('comment    -- Input conversion parameters ('+time.strftime("%c")+')\n')
ply.write('comment          compression_thickness    = ' +str(compression_thickness)+'\n')
ply.write('comment          average_breast_intensity = ' +str(average_breast_intensity)+'\n')
ply.write('comment          air_intensity  = ' +str(air_intensity)+'\n')
ply.write('comment          mfp_breast     = ' +str(mfp_breast)+'\n')
ply.write('comment          mfp_plastic    = ' +str(mfp_plastic)+'\n')
ply.write('comment          center_pixel   = ' +str(center_pixel)+'\n')
ply.write('comment          center_coord   = ' +str(center_coord)+'\n')
ply.write('comment          source_coord   = ' +str(source_coord)+'\n')
ply.write('comment          pixel_size     = ' +str(pixel_size)+'\n')
ply.write('comment          subtract_layer = ' +str(subtract_layer)+'\n')
ply.write('comment          output_binning = ' +str(output_binning)+'\n')
ply.write('comment          output_mm      = ' +str(output_mm)+'\n')
ply.write('comment          output_base    = ' +str(output_base)+'\n')
ply.write('comment          binned columns = ' +str(num_columns/output_binning)+", binned rows = "+str(num_rows/output_binning)+'\n')
ply.write('comment \n')

# - Report variables for the water fill:        !!WaterFill!!  
if (min_height<50.0):
  border_pixels = output_binning * int(border_cm/(pixel_size*output_binning) + 0.5)   # convert cm into pixels
  pix_max=0.01; pix_max_water=0.0; i_max=0; j_max=0; pix_counter=0; volume_subtracted=0.0; volume_kept=0.0; overfill_warning=0
  ply.write('comment     !!WaterFill!!\n')
  ply.write('comment         height_water_fill = '+str(height_water_fill)+' cm, min_height = '+str(min_height)+' cm, border = '+str(border_cm)+' cm\n')   #!!WaterFill!!  
  ply.write('comment         mfp_water = ' +str(mfp_water)+' cm\n')
  ply.write('comment \n')

ply.write('element vertex '+str(num_vertices)+'\n')
ply.write('property float x\n')
ply.write('property float y\n')
ply.write('property float z\n')
ply.write('element face '+str(num_triangles)+'\n')
ply.write('property list uchar int vertex_index\n')
ply.write('end_header\n')

sys.stdout.write(' Current row: ')

for j in range(rows[0], rows[1], output_binning):
  sys.stdout.write(str(j)+', '); sys.stdout.flush()
  for i in range(columns[0], columns[1], output_binning):    

    if pix[j][i]>-999990.01:     # Skip air pixels
      # Define vectors:
      r0 = numpypy.zeros(3)
      #r1 = numpypy.zeros(3)
      d  = numpypy.zeros(3)
      vertices = numpypy.zeros((8,3))
      offset = 0.5*pixel_size*output_binning
      offset0 = array([+offset,+offset,0.0])      #  3 ----- 0   (vertex order clockwise)
      offset1 = array([+offset,-offset,0.0])      #    - c -
      offset2 = array([-offset,-offset,0.0])      #  2 ----- 1
      offset3 = array([-offset,+offset,0.0])       
      
      # Set pixel center coordinates in real space:

      r0[0] = center_coord[0] + float(i-center_pixel[0]) * pixel_size    
      r0[1] = center_coord[1] + float(j-center_pixel[1]) * pixel_size
      r0[2] = center_coord[2]

      vertices[0] = r0 + offset0   # 4 vertices at the base
      vertices[1] = r0 + offset1
      vertices[2] = r0 + offset2
      vertices[3] = r0 + offset3
      
      # - Shorten object for printing (and make sure all pixels are 0 or more):
      pix[j][i] = pix[j][i] - subtract_layer    # !!PYPY!! Old: pix[j][i] = max(pix[j][i] - subtract_layer, 0.0)
      if (pix[j][i] < 0):
        pix[j][i] = 0



#!!WaterFill!! 
      # - Create a container for water in the middle of the phantom. Water attenuates more than ABS!       #!!WaterFill!!      
      if (pix[j][i] > min_height):   # Check minimum height
        if (j>=rows[0]+border_pixels) and (i>=columns[0]+border_pixels) and (j<rows[1]-border_pixels) and (i<columns[1]-border_pixels):    # Check pixel interval
          if (pix[j+border_pixels][i+border_pixels]>0.1) and (pix[j-border_pixels][i+border_pixels]>0.1) and (pix[j+border_pixels][i-border_pixels]>0.1) and (pix[j-border_pixels][i-border_pixels]>0.1):   # Check if not in periphery
            
            # - Process internal pixels:
            
            #   Compute the effective water thickness at current pixel in the slanted path to the source:
            d = source_coord-r0
            norm = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])
            height_water_fill_pixel = height_water_fill/source_height_Z*norm
            #print '   height_water_fill='+str(height_water_fill)+', height_water_fill_pixel='+str(height_water_fill_pixel)


            pix_tmp = pix[j][i]
            pix[j][i] = ( height_water_fill_pixel/mfp_water - pix[j][i]/mfp_plastic ) / ( 1.0/mfp_water - 1.0/mfp_plastic )    # Compute new pixel plastic height assuming filled with water up to height_water_fill   !!WaterFill!!
            
            pix_counter+=1
            volume_subtracted+=pix_tmp-pix[j][i]
            volume_kept+=pix[j][i]
            if (pix_tmp>pix_max):
              pix_max=pix_tmp; pix_max_water=pix[j][i]; i_max=i; j_max=j   # Store thickest plastic pixel to report
            if ((pix[j][i]<0.0) and (overfill_warning==0)):
              overfill_warning=1
 #!!WaterFill!! 
 

      d = source_coord-vertices[0]    # 4 vertices at the top focused to the source: calculate mormalized direction vector for each vertex
      norm = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])    # !!PYPY!! Not using function: numpypy.linalg.norm(d)
      vertices[4] = vertices[0] + pix[j][i]*d/norm
      d = source_coord-vertices[1]
      norm = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])
      vertices[5] = vertices[1] + pix[j][i]*d/norm
      d = source_coord-vertices[2]
      norm = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])
      vertices[6] = vertices[2] + pix[j][i]*d/norm
      d = source_coord-vertices[3]
      norm = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])
      vertices[7] = vertices[3] + pix[j][i]*d/norm

      # Write vertices coordinates:
      for n in range(0, 8):
        
        if (output_mm==1):
          ply.write(str(10.0*vertices[n][0])+' '+str(10.0*vertices[n][1])+' '+str(10.0*vertices[n][2])+'\n')  # Outputing mesh in mm
        else:
          ply.write(str(vertices[n][0])+' '+str(vertices[n][1])+' '+str(vertices[n][2])+'\n')                 # Outputing mesh in cm
    
sys.stdout.write('\n'); sys.stdout.flush()


# -- Report info on water fill:           !!WaterFill!!
if(pix_max>0.02):
  pixel_area=(pixel_size*output_binning)*(pixel_size*output_binning)
  print "\n!!WaterFill!!  Input:  height_water_fill="+str(height_water_fill)+", min_height=" + str(min_height) + ", border_pixels=" + str(border_pixels) + ", mfp_water=" + str(mfp_water)  + ", mfp_plastic=" + str(mfp_plastic)
  print "               Number pixels transformed = "+str(pix_counter)+" = "+str(pix_counter*pixel_area)+" cm^2. Plastic volume = "+str(volume_kept*pixel_area)+", Water volume = "+str((height_water_fill*pix_counter-volume_kept)*pixel_area)+", Volume saved = "+str(volume_subtracted*pixel_area)+" cm^3"
  print "               Original plastic pix_max = "+str(pix_max)+" cm, plastic thickness with water = pix["+str(j_max)+"]["+str(i_max)+"] = "+str(pix[j_max][i_max])
  print "               exp(-pix_max/mfp_plastic) = "+str(exp(-pix_max/mfp_plastic))+", exp(-pix/mfp_plastic-(height_water_fill-pix)/mfp_water)="+str(exp(-pix[j_max][i_max]/mfp_plastic-(height_water_fill-pix[j_max][i_max])/mfp_water))  
  print "               (exponentials differ bc pixel water thickness >= height_water_fill: we account for slanted path to source)"  
  if(pix[j_max][i_max]<0) or (overfill_warning!=0):
    print "\n***WARNING*** The water subtraction failed! Input plastic/water thicknesses combination not valid (results in negative water thickness).   !!WaterFill!!\n\n"
#!!WaterFill!! 


# Write triangles in the 6 sides of each printed column giving groups of 3 vertices:
for n in range(0, num_binned_pixels):

  # - Create a triangle mesh combining the walls of the columns of consecutive pixels:
  
  if (output_base==1):
    ply.write('3 '+str(0+n*8)+' '+str(1+n*8)+' '+str(2+n*8)+'\n')
    ply.write('3 '+str(2+n*8)+' '+str(3+n*8)+' '+str(0+n*8)+'\n')     #  3 ----- 0   (vertex order clockwise)    
  ply.write('3 '+str(4+n*8)+' '+str(5+n*8)+' '+str(6+n*8)+'\n')       #    - c -
  ply.write('3 '+str(6+n*8)+' '+str(7+n*8)+' '+str(4+n*8)+'\n')       #  2 ----- 1
  
  if (n<(num_columns/output_binning)):                                #  7 ----- 4
    ply.write('3 '+str(1+n*8)+' '+str(5+n*8)+' '+str(6+n*8)+'\n')     #    -top-
    ply.write('3 '+str(6+n*8)+' '+str(2+n*8)+' '+str(1+n*8)+'\n')     #  6 ----- 5
    
  elif (n>(num_binned_pixels-(num_columns/output_binning)-1)):
    ply.write('3 '+str(1+n*8+3-(num_columns/output_binning)*8)+' '+str(5+n*8)+' '+str(6+n*8)+'\n')
    ply.write('3 '+str(6+n*8)+' '+str(2+n*8+5-(num_columns/output_binning)*8)+' '+str(1+n*8+3-(num_columns/output_binning)*8)+'\n')
    ply.write('3 '+str(3+n*8)+' '+str(7+n*8)+' '+str(4+n*8)+'\n')
    ply.write('3 '+str(4+n*8)+' '+str(0+n*8)+' '+str(3+n*8)+'\n')
  else:
    ply.write('3 '+str(1+n*8+3-(num_columns/output_binning)*8)+' '+str(5+n*8)+' '+str(6+n*8)+'\n')
    ply.write('3 '+str(6+n*8)+' '+str(2+n*8+5-(num_columns/output_binning)*8)+' '+str(1+n*8+3-(num_columns/output_binning)*8)+'\n')
    #ply.write('3 '+str(3+n*8)+' '+str(7+n*8)+' '+str(4+n*8)+'\n')      #  No need to build this wall, already written in past column
    #ply.write('3 '+str(4+n*8)+' '+str(0+n*8)+' '+str(3+n*8)+'\n')  
  
  if ((n+1)%(num_columns/output_binning))==0:   #!!ReducedTriangles!!    
    ply.write('3 '+str(0+n*8)+' '+str(4+n*8)+' '+str(5+n*8)+'\n')
    ply.write('3 '+str(5+n*8)+' '+str(1+n*8)+' '+str(0+n*8)+'\n')
  elif (n%(num_columns/output_binning))==0:       #!!ReducedTriangles!!     
    ply.write('3 '+str(0+n*8+7+8)+' '+str(4+n*8)+' '+str(5+n*8)+'\n')
    ply.write('3 '+str(5+n*8)+' '+str(1+n*8+5+8)+' '+str(0+n*8+7+8)+'\n')
    ply.write('3 '+str(2+n*8)+' '+str(6+n*8)+' '+str(7+n*8)+'\n')
    ply.write('3 '+str(7+n*8)+' '+str(3+n*8)+' '+str(2+n*8)+'\n')
  else:    
    ply.write('3 '+str(0+n*8+7+8)+' '+str(4+n*8)+' '+str(5+n*8)+'\n')   # Connect top right with top left next column (7)   #!!ReducedTriangles!!
    ply.write('3 '+str(5+n*8)+' '+str(1+n*8+5+8)+' '+str(0+n*8+7+8)+'\n')
    #ply.write('3 '+str(2+n*8)+' '+str(6+n*8)+' '+str(7+n*8)+'\n')      #  No need to build this wall, already written in past column
    #ply.write('3 '+str(7+n*8)+' '+str(3+n*8)+' '+str(2+n*8)+'\n')
  
# ply.write(str(n)+'\n')   # !!DeBuG!! Why writing str(n)?? Error?
ply.write('\n')
ply.close()
