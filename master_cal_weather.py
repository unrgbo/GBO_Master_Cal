#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 07:36:28 2018

@author: jfausett
"""

import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import robust as rb
#from length import *
import math
from tqdm import tqdm

# This is a test comment for looking at branch

# This is a test comment to check that file is pushed to gitlab
# =============================================================================
# Make Calibration Frames:
# =============================================================================
def make_cals(bias=False,dark=False,flat=False,
              bins=['1X1','2X2','3X3','4X4'],
              temps=['-25','-30','-35'],
              exposures=['030','060','120','180','240','300','600'],
              filters=['B','Blue',"g'",'Grating','Green','Ha','I',"i'",
                       'Luminance','OIII','R',"r'",'Red','SII','V','Y',"z'"],
              root ='Z:\Calibration Files',
              suffix='.fts'):
              #root ='/Users/jfausett/Dropbox/Calibration Files/',suffix='.fts'):
    """
    Overview:
    ---------
    Create master calibration frames from a large dataset. 
    Output is written to "outdir" in FITS format.

    Inputs:
    -------
    root        : Parent directory to search for files (default 'Z:\Calibration Files')
    suffix      : Type of extension to search for (default '.fts')
    bins        : Type of Binning to use (default ['1X1','2X2','3X3'])
    temps       : CCD temperature setting (default ['-25','-35'])
    exposures   : Dark exposure times (default ['030','060','120','180','240','300','600'])
    
    Keyword arguments:
    ------------------
    bias        : Weather or not to create master bias (default False)
    dark        : Weather or not to create master darks (default False)
    flat        : Weather or not to create master flats (default False)

    Calling sequence:
    -----------------
    make_cals(bias=True,dark=True,flat=True)

    """

    mstdir = 'Z:\Calibration Master Frames\\'
    #mstdir = '/Users/jfausett/Dropbox/Calibration Master Frames/'
    if not os.path.exists(mstdir):
        os.makedirs(mstdir)
        print 'Creating master directory'

    #  Generate a list of all files in root directory with given suffix
    print '\nFinding all calibration files\n'
    list_of_files = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            if name.endswith(suffix):
                list_of_files.append(os.path.join(path, name))
                
    print 'Found %d calibration files \n' % len(list_of_files)
    
    types = []                
    names = []
    temp = []
    binning = []
    JD = []
    exposure = []
    bands = []
    bad_files = []
    no_temp = []
    
    #  Determine weather or not a dataframe has already been created
    if os.path.exists(mstdir+'all_files.pkl'):
        done_files = pd.read_pickle(mstdir+'all_files.pkl')
        
        print 'Found existing dataframe with %d files' % len(done_files)
        done_list = done_files['name'].tolist()
        new_list = [x for x in list_of_files if x not in done_list]
        #  Append existing data frame with info from new files
        if len(new_list) > 0:
            print 'Creating dataframe for new files'
            for fname in tqdm(new_list):
                with fits.open(fname) as hdu:
                    types.append(hdu[0].header['IMAGETYP'])
                    names.append(fname)
                    try:
                        temp.append(str(int(hdu[0].header['SET-TEMP'])))
                    except:
                        print 'No temperature keyword for '+fname+'\nFilling with NaN'
                        temp.append(np.nan)
                        no_temp.append(fname)
                    binning.append(str(hdu[0].header['XBINNING'])+'X'+str(hdu[0].header['YBINNING']))
                    JD.append(int(hdu[0].header['JD']))
                    try:
                        exposure.append(str(int(hdu[0].header['EXPOSURE'])))
                    except:
                        exposure.append(np.nan)
                    try:
                        bands.append(hdu[0].header['FILTER'])
                    except:
                        bands.append(np.nan)

            d = {'type': types, 'name': names, 'temp': temp, 'binning':
                binning, 'JD': JD, 'exp': exposure, 'filter': bands}
            new_files = pd.DataFrame(data=d)
        
            del types,names,temp,binning,JD,exposure,bands

            done_files = pd.concat([done_files, new_files], ignore_index=True)
            
            print '\nWriting out new dataframe with %d files' % len(done_files)
            done_files.to_pickle(mstdir+'all_files.pkl')
            print '\nHere are the files without temperature info in header:\n'
            print no_temp
        
        del done_files
        
    #  Create a master dataframe for all files
    else:
        print 'Creating master dataframe to parse \n'
        for fname in tqdm(list_of_files):
            try:
                with fits.open(fname) as hdu:
                    types.append(hdu[0].header['IMAGETYP'])
                    names.append(fname)
                    try:
                        temp.append(str(int(hdu[0].header['SET-TEMP'])))
                    except:
                        print 'No temperature keyword for '+fname+'\nFilling with NaN'
                        temp.append(np.nan)
                        no_temp.append(fname)
                    binning.append(str(hdu[0].header['XBINNING'])+'X'+str(hdu[0].header['YBINNING']))
                    JD.append(int(hdu[0].header['JD']))
                    try:
                        exposure.append(str(int(hdu[0].header['EXPOSURE'])))
                    except:
                        exposure.append(np.nan)
                    try:
                        bands.append(hdu[0].header['FILTER'])
                    except:
                        bands.append(np.nan)
            except:
                print '\nThere is a problem with the file:\n'+fname
                bad_files.append(fname)
                
        print '\nHere are the files without temperature info in header:\n'
        print no_temp

        print '\nHere is a list of the files that raised errors:\n'
        print bad_files
        
        d = {'type': types, 'name': names, 'temp': temp, 'binning':
            binning, 'JD': JD, 'exp': exposure, 'filter': bands}
        all_files = pd.DataFrame(data=d)
    
        del types,names,temp,binning,JD,exposure,bands,bad_files
            
        print 'Writing out dataframe to : '+mstdir+'all_files.pkl'
        all_files.to_pickle(mstdir+'all_files.pkl')
        
        del all_files

    if bias:
        
        all_files = pd.DataFrame()
        all_files = pd.read_pickle(mstdir+'all_files.pkl')
        
        bias_files = all_files.where(all_files['type'] == 'Bias Frame')
        bias_files.dropna(how='all',inplace=True)   
        
        del all_files
        
        dates = np.unique(bias_files.loc[:,'JD'])
        dates = [int(x) for x in dates]  
        
        for date in dates:
            year, month, day = jd_to_date(date)
            day = int(day)
            year, month, day = str(year), str(month), str(day)
            if len(month) == 1:
                month = '0'+month
            if len(day) == 1:
                day = '0'+day
            tagdate = year+month+day
            print '\nChecking frames for '+tagdate
            for tmp in temps:
                for bns in bins:
                    filter1 = bias_files['temp'] == tmp
                    filter2 = bias_files['binning'] == bns
                    filter3 = bias_files['JD'] == date
                    
                    files = bias_files.where(filter1 & filter2 & filter3)
        
                    files.dropna(how='all', inplace=True)

                    fnames = files['name'].tolist()

                    fct = len(fnames)

                    if fct > 50:
                        print 'There are %d total files\n'%fct
                        print 'Too many files to make master\n'
                        break

                    elif fct > 35 and fct <= 50 and bns == '1X1':
                        sub_frames = [fnames[x:x+15] for x in xrange(0, len(fnames), 15)]
                        print 'Creating sub_master_bias for '+tagdate+' with '+bns+' and a temperature of '+tmp
                        if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                            os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                         if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                             os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                        for i in range(3):
                            master_bias(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_%d'% (i+1))
                            #master_bias(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_%d'% (i+1))
                    elif fct>=20 and fct <= 35 and bns == '1X1':
                        sub_frames = [fnames[x*(fct/2):(x+1)*(fct/2)] for x in range((len(fnames)+ (fct/2)-1)//(fct/2))]
                        print 'Creating sub_master_bias for '+tagdate+' with '+bns+' and a temperature of '+tmp
                        if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                            os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                         if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                             os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                        for i in range(2):
                            master_bias(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_%d'% (i+1))
                            #master_bias(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_%d'% (i+1))
                    elif fct > 1:
                        if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                            os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                         if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                             os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                        print 'Creating master_bias for '+tagdate+' with '+bns+' and a temperature of '+tmp
                        master_bias(files=fnames, outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns)
                        #master_bias(files=fnames, outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns)
        del bias_files
    
    if dark:
        
        all_files = pd.DataFrame()
        all_files = pd.read_pickle(mstdir+'all_files.pkl')
        
        dark_files = all_files.where(all_files['type'] == 'Dark Frame')
        dark_files.dropna(how='all',inplace=True)

        del all_files

        dates = np.unique(dark_files.loc[:,'JD'])
        dates = [int(x) for x in dates]

        for date in dates:
            year, month, day = jd_to_date(date)
            day = int(day)
            year, month, day = str(year), str(month), str(day)
            if len(month) == 1:
                month = '0'+month
            if len(day) == 1:
                day = '0'+day
            tagdate = year+month+day
            print '\nChecking frames for '+tagdate
            for tmp in temps:
                for bns in bins:
                    for exp in exposures:
                        filter1 = dark_files['temp'] == tmp
                        filter2 = dark_files['binning'] == bns
                        filter3 = dark_files['JD'] == date
                        filter4 = dark_files['exp'] == exp
                        
                        files = dark_files.where(filter1 & filter2 & filter3 & filter4)
                        files.dropna(how='all',inplace=True)
                    
                        fnames = files['name'].tolist()

                        fct = len(fnames)

                        if fct > 50:
                            print 'There are %d total files\n'%fct
                            print 'Too many files to make master\n'
                            break

                        elif fct > 35 and fct <= 50 and bns == '1X1':
                            print 'There are %d total files\n'%fct
                            sub_frames = [fnames[x:x+15] for x in xrange(0, len(fnames), 15)]
                            print 'Creating sub_master_dark for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                            for i in range(3):
                                master_dark(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_'+exp+'_%d'% (i+1))
                                #master_dark(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_'+exp+'_%d'% (i+1)
                        elif fct>=20 and fct <= 35 and bns == '1X1':
                            print 'There are %d total files\n'%fct
                            sub_frames = [fnames[x*(fct/2):(x+1)*(fct/2)] for x in range((len(fnames)+ (fct/2)-1)//(fct/2))]
                            print 'Creating sub_master_dark for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                            for i in range(2):
                                master_dark(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_'+exp+'_%d'% (i+1))
                                #master_dark(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_'+exp+'_%d'% (i+1)
                 
                        elif fct > 1:
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                            print 'Creating master_dark for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            master_dark(files=fnames, outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_'+exp)
                            #master_dark(files=fnames, outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_'+exp)
        del dark_files
                            
    if flat:

        all_files = pd.DataFrame()
        all_files = pd.read_pickle(mstdir+'all_files.pkl')
        
        flat_files = all_files.where(all_files['type'] == 'Flat Field')
        flat_files.dropna(how='all',inplace=True)     
        
        del all_files

        dates = np.unique(flat_files.loc[:,'JD'])
        dates = [int(x) for x in dates]

        for date in dates:
            year, month, day = jd_to_date(date)
            day = int(day)
            year, month, day = str(year), str(month), str(day)
            if len(month) == 1:
                month = '0'+month
            if len(day) == 1:
                day = '0'+day
            tagdate = year+month+day
            print '\nChecking frames for '+tagdate
            for tmp in temps:
                for bns in bins:
                    for band in filters:
                        filter1 = flat_files['temp'] == tmp
                        filter2 = flat_files['binning'] == bns
                        filter3 = flat_files['JD'] == date
                        filter4 = flat_files['filter'] == band
                        
                        files = flat_files.where(filter1 & filter2 & filter3 & filter4)
                        files.dropna(how='all',inplace=True)
                    
                        fnames = files['name'].tolist()

                        fct = len(fnames)
                    
                        if fct > 50:
                            print 'There are %d total files\n'%fct
                            print 'Too many files to make master\n'
                            break

                        elif fct > 35 and fct <= 50 and bns == '1X1':
                            print 'There are %d total files\n'%fct
                            sub_frames = [fnames[x:x+15] for x in xrange(0, len(fnames), 15)]
                            print 'Creating sub_master_flat for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                            for i in range(3):
                                master_flat(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_%d'% (i+1), band=band)
                                #master_flat(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_%d'% (i+1), band=band)
                        elif fct>=20 and fct <= 35 and bns == '1X1':
                            print 'There are %d total files\n'%fct
                            sub_frames = [fnames[x*(fct/2):(x+1)*(fct/2)] for x in range((len(fnames)+ (fct/2)-1)//(fct/2))]
                            print 'Creating sub_master_flat for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================
                            for i in range(2):
                                master_flat(files=sub_frames[i], outdir=mstdir+tagdate+'\\'+bns+'\\', tag=tagdate+'_'+bns+'_%d'% (i+1), band=band)
                                #master_flat(files=sub_frames[i], outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns+'_%d'% (i+1), band=band)
                                           
                        if fct > 1:
                            if not os.path.exists(mstdir+tagdate+'\\'+bns+'\\'):
                                os.makedirs(mstdir+tagdate+'\\'+bns+'\\')
# =============================================================================
#                             if not os.path.exists(mstdir+tagdate+'/'+bns+'/'):
#                                 os.makedirs(mstdir+tagdate+'/'+bns+'/')
# =============================================================================                                
                            print '\nCreating master_flat for '+tagdate+' with '+bns+' and a temperature of '+tmp
                            master_flat(files=fnames, outdir=mstdir+tagdate+'/'+bns+'/', tag=tagdate+'_'+bns, band=band)
        del flat_files

#----------------------------------------------------------------------#
# master_bias:
#----------------------------------------------------------------------#
def master_bias(files,write=True,outdir='/',readnoise=False,clobber=False,verbose=True,
                float32=True,tag='',median=True):

    """
   Overview:
    ---------
    Create master bias frame from series of biases (median filter).
    Returns a master_bias frame and writes FITS file to disk in specified
    directory.

    Optionally, the read noise is calculated from the variance of each
    pixel in the bias stack. This is *very* slow. So only use this option
    if you really need to. The readnoise image is also written to disk.

    Inputs:
    -------
    files       : List of flat field files from which a master bias will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    readnoise   : Do readnoise calculation (very slow! default False)
    verbose     : Print out progress (default True)

    Calling sequence:
    -----------------
    master_bias = master_bias(biasfiles,write=True,readnoise=False,
                              outdir='/home/users/bob/stuff/')

    """

# Don't redo master_bias unless clobber keyword set
    name  = outdir+'master_bias_'+tag+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master bias already exists!")
        master_bias = fits.getdata(name,0,header=False)
        return master_bias

# Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))
    temps = []

# Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        temps.append(header["CCD-TEMP"])
        stack[i,:,:] = image

# Calculate read noise directly from bias frames if prompted
    if readnoise:
        rn = np.zeros((ysz,xsz))
        print("Starting readnoise calculation")
        pbar = tqdm(desc = 'Calculating readnoise', total = ysz, unit = 'rows')
        for i in np.arange(ysz):
            for j in np.arange(xsz):
                rn[i,j] = rb.std(stack[:,i,j])
            pbar.update(1)

# Make a nice plot (after all that hard work)
        aspect = np.float(xsz)/np.float(ysz)
        plt.figure(39,figsize=(5*aspect*1.2,5))
        plt.clf()
        sig = rb.std(rn)
        med = np.median(rn)
        mean = np.mean(rn)
        vmin = med - 2*sig
        vmax = med + 2*sig
        plt.imshow(rn,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
        plt.colorbar()
        plt.annotate(r'$\bar{\sigma}$ = %.2f cts' % mean, [0.95,0.87],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.SimpleLineShadow(linewidth=3,foreground="w")])
        plt.annotate(r'med($\sigma$) = %.2f cts' % med, [0.95,0.8],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.annotate(r'$\sigma_\sigma$ = %.2f cts' % sig,
                     [0.95,0.73],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.title("Read Noise")
        plt.xlabel("pixel number")
        plt.ylabel("pixel number")

        if write:
            plt.savefig(outdir+'readnoise'+tag+'.png',dpi=300)

# Calculate master bias frame by median filter
    print('Calculating median of stacked frames...')
    if median:
        master_bias = np.median(stack,axis=0)
    else:
        master_bias = np.mean(stack,axis=0)
        
    # Make a plot
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(38,figsize=(5*aspect*1.2,5))
    plt.clf()
    sig = rb.std(master_bias)
    med = np.median(master_bias)
    vmin = med - 2*sig
    vmax = med + 2*sig
    plt.imshow(master_bias,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.annotate('Bias Level = %.2f cts' % med, [0.95,0.87],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large',color='k')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts' % sig, [0.95,0.8],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.95,0.73],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Bias")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

# Write out bias, readnoise and plot
    if write:
        name  = outdir+'master_bias_'+tag
        plt.savefig(name+'.png',dpi=300)

        hout = fits.Header()
        hout['CCDTEMP'] = (np.median(temps), "Median CCD temperature")
        hout["TEMPSIG"] = (np.std(temps), "CCD temperature RMS")
        hout["BIAS"] = (med, "Median bias level (cts)")
        hout["BIASSIG"] = (sig, "Bias RMS (cts)")
        if len(glob.glob(name+'.fits')) == 1:
            os.system('rm '+name+'.fits')
        if float32:
            fits.writeto(name+'.fits', np.float32(master_bias), hout)
        else:
            fits.writeto(name+'.fits', master_bias, hout)

        if readnoise:
            name  = outdir+'readnoise'+tag
            if len(glob.glob(name+'.fits')) == 1:
                os.system('rm '+name+'.fits')
            if float32:
                fits.writeto(name+'.fits', np.float32(rn), hout)
            else:
                fits.writeto(name+'.fits', rn, hout)

    return master_bias

#----------------------------------------------------------------------#
# master_dark:
#----------------------------------------------------------------------#
def master_dark(files,bias=None,write=True,outdir='/',clobber=False,float32=True,tag='',
                median=True):
    """
    Overview:
    ---------
    Create master dark frame from series of darks (median filter).
    Returns a master dark frame. If write is specified, a FITS file
    will be written to "outdir" (default is pwd).

    Inputs:
    -------
    files       : List of flat field files from which a master dark will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)

    Calling sequence:
    -----------------
    master_dark = master_dark(darkfiles,bias=master_bias,write=True,
                              outdir='/home/users/bob/stuff/')

    """

# Don't redo master_dark unless clobber keyword set
    name  = outdir+'master_dark_'+tag+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master dark already exists!")
        master_dark = fits.getdata(name,0,header=False)
        return master_dark

 # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))
    temps = []
    exps = []

# Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        exp = header["EXPOSURE"]
        exps.append(exp)
        temps.append(header["CCD-TEMP"])   
#        if length(bias) == 1:
#        image = np.float(image)/exp
# =============================================================================
#         else:
#             image = (image-bias)/exp
# =============================================================================
        stack[i,:,:] = image

# Obtain statistics for the master dark image header
    # Temperature
    tmax = np.max(temps)
    tmin = np.min(temps)
    tmean = np.mean(temps)
    tmed = np.median(temps)
    tsig = np.std(temps)
    # Exposure times
    expmax = np.max(exps)
    expmin = np.min(exps)
    print('')
    print("Minimum CCD Temp. %.2f C" % tmin)
    print("Maximum CCD Temp. %.2f C" % tmax)
    print("CCD Temp. rms: %.3f C" % tsig)
    print("CCD Temp. mean: %.2f C" % tmean)
    print("CCD Temp. median: %.2f C" % tmed)

# Create master dark by median filter or mean
    if median:
        master_dark = np.median(stack,axis=0)
    else:
        master_dark = np.mean(stack,axis=0)
        
# Make a plot
    sig = rb.std(master_dark)
    med = np.median(master_dark)
    vmin = med - 2*sig
    vmax = med + 2*sig
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(37,figsize=(5*aspect*1.2,5))
    plt.clf()
    plt.imshow(master_dark,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.annotate('Dark Current = %.2f cts/sec' % med, [0.72,0.8],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts/sec' % sig, [0.72,0.75],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.72,0.7],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Dark")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

# Write out plot and master dark array
    if write:
        name = outdir+'master_dark_'+tag

        plt.savefig(name+'.png',dpi=300)

        hout = fits.Header()
        hout["TEMPMAX"] = (tmax, "Maximum CCD temperature")
        hout["TEMPMIN"] = (tmin, "Minimum CCD temperature")
        hout["TEMPMED"] = (tmed, "Median CCD temperature")
        hout["TEMPMN"] = (tmean, "Mean CCD temperature")
        hout["TEMPSIG"] = (tsig, "CCD temperature RMS")
        hout["EXPMAX"] = (expmax,"Maximum exposure time")
        hout["EXPMIN"] = (expmin, "Minimum exposure time")
        hout["DARKCNT"] = (med, "Median dark current (cts/sec)")
        hout["DARKSIG"] = (sig, "Dark current RMS (cts/sec)")
        if len(glob.glob(name)) == 1:
            os.system('rm '+name+'.fits')
        if float32:
            fits.writeto(name+'.fits', np.float32(master_dark), hout)
        else:
            fits.writeto(name+'.fits', master_dark, hout)

    return master_dark

#----------------------------------------------------------------------#
# master_flat:
#----------------------------------------------------------------------#
def master_flat(files,bias=None,dark=None,write=True,outdir='/',band='V',
                tag='',clobber=False,stretch=3,float32=True,median=True):

    """
    Overview:
    ---------
    Create a master flat using (optionally) a provided bias and dark frame. Output
    is written to "outdir" in FITS format.

    Inputs:
    -------
    files       : List of flat field files from which a master flat will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    dark        : Master dark frame calibrated in ADU/sec (default None)
    band        : Band from which flatis being produced
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    stretch     : Multiple of the noise RMS to stretch image (default 3)


    Calling sequence:
    -----------------
    master_flat = master_flat(flatfiles,bias=master_bias,dark=master_dark,write=True,
                              outdir='/home/users/bob/stuff/')

    """

# Don't redo master_flat unless clobber keyword set
    name = outdir+'master_flat_'+tag+'_'+band+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master flat already exists!")
        master_flat = fits.getdata(name,0, header=False)
        return master_flat

 # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    filter = header["filter"]

    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))

# Load stack array and get CCD temperatures
    meds = []
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        image = np.float32(image)
        if header["filter"] != filter:
            sys.exit("Filters do not match!")
# =============================================================================
#         if length(bias) > 1:
#             image -= bias
#         if length(dark) > 1:
#             exptime = header['EXPTIME']
#             image -= dark*exptime
# =============================================================================
        meds.append(np.median(image))
        stack[i,:,:] = image/np.median(image)

# Obtain statistics for the master dark image header
    med = np.median(meds)
    sig = np.std(meds)

# Create master flat by median filter
    if median:
        master_flat = np.median(stack,axis=0)
    else:
        master_flat = np.mean(stack,axis=0)
        
# Make a plot
    sig = rb.std(master_flat)
    med = np.median(master_flat)
    vmin = med - stretch*sig
    vmax = med + stretch*sig
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(40,figsize=(5*aspect*1.2,5))
    plt.clf()
    plt.imshow(master_flat,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.title("Master Flat")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

# Write out plot and master flat array
    if write:
        plt.savefig(outdir+'master_flat'+tag+'_'+band+'.png',dpi=300)
        hout = fits.Header()
        hout["FILTER"] = (filter, "Filter used when taking image")
        hout["MEDCTS"] = (med, "Median counts in individual flat frames")
        hout["MEDSIG"] = (sig, "Median count RMS in individual flat frames")
# =============================================================================
#         if length(bias) > 1:
#             hout.add_comment("Bias subtracted")
#         if length(dark) > 1:
#             hout.add_comment("Dark subtracted")
# =============================================================================

        if len(glob.glob(outdir+'master_flat'+tag+'.fits')) == 1:
            os.system('rm '+outdir+'master_flat'+tag+'.fits')
        if float32:
            fits.writeto(outdir+'master_flat_'+tag+'_'+band+'.fits', np.float32(master_flat), hout)
        else:
            fits.writeto(outdir+'master_flat_'+tag+'_'+band+'.fits', master_flat, hout)

    return master_flat

# =============================================================================
# Date Convert:
# =============================================================================
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day
