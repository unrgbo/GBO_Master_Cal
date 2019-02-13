#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 07:36:28 2018

@author: jfausett
"""

import sys, os, glob, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import robust as rb
# from length import *
import math
from tqdm import tqdm
from jdcal import gcal2jd, jd2gcal
import itertools

# This is a test comment for looking at branch

# This is a test comment to check that file is pushed to gitlab
# =============================================================================
# Make Calibration Frames:
# =============================================================================
def make_cals(bias=False, dark=False, flat=False,
              root='Z:\Calibration Files\\',
              suffix='.fts',
              clobber=False):
    # root ='/Users/jfausett/Dropbox/Calibration Files/',
    # suffix='.fts'):
    # root = '/Users/jakrin/Calibration Files/',
    # suffix='.fts'):
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
    bias        : Whether or not to create master bias (default False)
    dark        : Whether or not to create master darks (default False)
    flat        : Whether or not to create master flats (default False)

    Calling sequence:
    -----------------
    make_cals(bias=True,dark=True,flat=True)

    """

    sjd = 2400000.5

    mstdir = 'Z:\Calibration Master Frames\\'
    # mstdir = '/Users/jfausett/Dropbox/Calibration Master Frames/'
    # mstdir = '/Users/jakrin/Calibration Master Frames/'
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
    tagdate = []
    exposure = []
    bands = []
    bad_files = []

    #  Determine weather or not a dataframe has already been created
    if os.path.exists(mstdir + 'all_files.pkl'):
        done_files = pd.read_pickle(mstdir + 'all_files.pkl')

        print 'Found existing dataframe with {} files'.format(len(done_files))
        done_list = done_files['name'].tolist()
        new_list = [x for x in list_of_files if x not in done_list]
        #  Append existing data frame with info from new files
        if len(new_list) > 0:
            print 'Creating dataframe for new files'
            for fname in tqdm(new_list):
                names.append(fname)
                try:
                    types.append(fits.getval(fname, 'IMAGETYP'))
                except:
                    types.append(np.nan)
                    bad_files.append(fname)
                try:
                    binning.append(str(fits.getval(fname, 'XBINNING')) + 'X' + str(fits.getval(fname, 'YBINNING')))
                except:
                    binning.append(np.nan)
                    bad_files.append(fname)
                try:
                    date = int(fits.getval(fname, 'JD'))
                    JD.append(date)
                    year, month, day, sec = jd2gcal(sjd, (date - sjd))
                    caldate = str(year) + str(month).zfill(2) + str(day).zfill(2)
                    tagdate.append(caldate)
                except:
                    JD.append(np.nan)
                    tagdate.append(np.nan)
                    bad_files.append(fname)
                try:
                    temp.append(str(int(fits.getval(fname, 'SET-TEMP'))))
                except:
                    temp.append(np.nan)
                    bad_files.append(fname)
                try:
                    exposure.append(str(int(fits.getval(fname, 'EXPOSURE'))).zfill(3))
                except:
                    exposure.append(np.nan)
                    bad_files.append(fname)
                try:
                    bands.append(fits.getval(fname, 'FILTER'))
                except:
                    bands.append(np.nan)
                    bad_files.append(fname)

            print len(types), len(names), len(temp), len(binning), len(JD), len(exposure), len(bands)

            if (len(names) != len(types) or len(names) != len(temp) or len(names) != len(binning)
                    or len(names) != len(JD) or len(names) != len(exposure) or len(names) != len(bands)):
                print 'Not all lists are the same length'
                return

            d = {'type': types, 'name': names, 'temp': temp, 'binning':
                binning, 'JD': JD, 'exp': exposure, 'filter': bands}

            new_files = pd.DataFrame(data=d)
            new_files.dropna(how='any', inplace=True)

            for fname in bad_files:
                baddir = '{}bad_files\\'.format(mstdir)
                if not os.path.exists(os.path.dirname(baddir)):
                    os.mkdir(os.path.dirname(baddir))
                shutil.move(fname, fname.replace(root, baddir))

            del types, names, temp, binning, JD, exposure, bands, bad_files

            done_files = pd.concat([done_files, new_files], ignore_index=True)

            print '\nWriting out new dataframe with {} files'.format(len(done_files))
            done_files.to_pickle(mstdir + 'all_files.pkl')

        del done_files, new_files

    #  Create a master dataframe for all files
    else:
        print 'Creating master dataframe to parse \n'
        for fname in tqdm(list_of_files):
            names.append(fname)
            try:
                types.append(fits.getval(fname, 'IMAGETYP'))
            except:
                types.append(np.nan)
                bad_files.append(fname)
            try:
                binning.append(str(fits.getval(fname, 'XBINNING')) + 'X' + str(fits.getval(fname, 'YBINNING')))
            except:
                binning.append(np.nan)
                bad_files.append(fname)
            try:
                date = int(fits.getval(fname, 'JD'))
                JD.append(date)
                year, month, day, sec = jd2gcal(sjd, (date - sjd))
                caldate = str(year) + str(month).zfill(2) + str(day).zfill(2)
                tagdate.append(caldate)
            except:
                JD.append(np.nan)
                tagdate.append(np.nan)
                bad_files.append(fname)
            try:
                temp.append(str(int(fits.getval(fname, 'SET-TEMP'))))
            except:
                temp.append(np.nan)
                bad_files.append(fname)
            try:
                exposure.append(str(int(fits.getval(fname, 'EXPOSURE'))).zfill(3))
            except:
                exposure.append(np.nan)
                bad_files.append(fname)
            try:
                bands.append(fits.getval(fname, 'FILTER'))
            except:
                bands.append(np.nan)
                bad_files.append(fname)

        print len(types), len(names), len(temp), len(binning), len(JD), len(exposure), len(bands)

        if (len(names) != len(types) or len(names) != len(temp) or len(names) != len(binning)
            or len(names) != len(JD) or len(names) != len(exposure) or len(names) != len(bands)):
            print 'Not all lists are the same length'
            return

        d = {'type': types, 'name': names, 'temp': temp, 'binning':
        binning, 'JD': JD, 'exp': exposure, 'filter': bands}

        all_files = pd.DataFrame(data=d)
        all_files.dropna(how='any', inplace=True)

        for fname in bad_files:
            baddir = '{}bad_files\\'.format(mstdir)
            if not os.path.exists(os.path.dirname(baddir)):
                os.mkdir(os.path.dirname(baddir))
            shutil.move(fname, fname.replace(root, baddir))

        del types, names, temp, binning, JD, exposure, bands, bad_files

        print 'Writing out dataframe to : {}all_files.pkl'.format(mstdir)
        all_files.to_pickle('{}all_files.pkl'.format(mstdir))

        del all_files

    del list_of_files

    if bias:

        all_files = pd.DataFrame()
        all_files = pd.read_pickle('{}all_files.pkl'.format(mstdir))

        bias_files = all_files.where(all_files['type'] == 'Bias Frame')
        bias_files.dropna(how='all', inplace=True)

        del all_files

        total_files = 0
        bins = bias_files['binning'].unique()
        for binning in bins:
            bin = bias_files.where(bias_files['binning'] == binning)
            bin.dropna(how='all', implace=True)
            dates = bin['tagdate'].unique()
            for date in tqdm(dates):
                biasdir = '{}{}\\Bias\\{}\\'.format(mstdir, binning, date)


                files = bin.where(bin['tagdate'] == date)
                files.dropna(how='all', inplace=True)

                fnames = files['name'].tolist()
                fct = len(fnames)

                if fct > 50:
                    print 'There are %d total files\n' % fct
                    print 'Too many files to make master\n'
                    pass
                elif fct > 35 and fct <= 50 and binning == '1X1':
                    sub_frames = [fnames[x:x + 15] for x in xrange(0, len(fnames), 15)]
                    print 'Creating sub_master_bias: binning is {}, date is {}'.format(binning, tagdate)
                    if not os.path.exists(biasdir):
                        os.makedirs(biasdir)
                    for i in range(3):
                        master_bias(files=sub_frames[i], outdir=biasdir, tag=tagdate + '_' + binning + '_%d' % (i + 1))
                elif fct >= 20 and fct <= 35 and binning == '1X1':
                    sub_frames = [fnames[x * (fct / 2):(x + 1) * (fct / 2)] for x in
                                  range((len(fnames) + (fct / 2) - 1) // (fct / 2))]
                    print 'Creating sub_master_bias: binning is {}, date is {}'.format(binning, tagdate)
                    if not os.path.exists(biasdir):
                        os.makedirs(biasdir)
                    for i in range(2):
                        master_bias(files=sub_frames[i], outdir=biasdir, tag=tagdate + '_' + binning + '_%d' % (i + 1))
                elif fct > 1:
                    if not os.path.exists(biasdir):
                        os.makedirs(biasdir)
                    print 'Creating master_bias: binning is {}, date is {}'.format(binning, tagdate)
                    master_bias(files=fnames, outdir=biasdir, tag=tagdate + '_' + binning)
                else:
                    print 'No Bias Frames found'

                total_files += fct

        print 'There are {} files of type: "Bias Frame"'.format(len(bias_files))
        print 'Processed {} bias files'.format(total_files)
        del bias_files

    if dark:

        all_files = pd.DataFrame()
        all_files = pd.read_pickle(mstdir + 'all_files.pkl')

        dark_files = all_files.where(all_files['type'] == 'Dark Frame')
        dark_files.dropna(how='all', inplace=True)

        exposures = dark_files['exp'].tolist()
        exposures = np.unique(exposures)
        exposures = [str(int(x)).zfill(3) for x in exposures]

        print exposures

        del all_files

        dates = dark_files['JD'].tolist()
        dates = [int(x) for x in dates]
        dates = np.unique(dates)

        print '\nParsing existing data frames to determine new Dark data\n'
        donedates = []
        for date in tqdm(dates):
            year, month, day, sec = jd2gcal(sjd, (date - sjd))
            year, month, day = str(year), str(month).zfill(2), str(day).zfill(2)
            tagdate = year + month + day
            for bns in bins:
                for exp in exposures:
                    darkpath = '{}{}\\Dark\\{}'.format(mstdir, bns, tagdate)
                    if os.path.exists(darkpath):
                        donedates.append(date)

        newdates = [x for x in dates if x not in donedates]
        if clobber:
            newdates = dates

        if len(newdates) == 0:
            print '\nDid not find any new dates with Dark Frames\n'
        else:
            print '\nBeginning master dark creation for new dates\n'
        for date in tqdm(newdates):
            year, month, day, sec = jd2gcal(sjd, (date - sjd))
            year, month, day = str(year), str(month).zfill(2), str(day).zfill(2)
            tagdate = year + month + day
            for bns in bins:
                for exp in exposures:

                    darkdir = mstdir + bns + '\\Dark\\' + tagdate + '\\'
                    # darkdir = mstdir+bns+'/Dark/'+tagdate+'/'

                    filter1 = dark_files['binning'] == bns
                    filter2 = dark_files['JD'] == date
                    filter3 = dark_files['exp'] == exp

                    files = dark_files.where(filter1 & filter2 & filter3)
                    files.dropna(how='any', inplace=True)

                    fnames = files['name'].tolist()
                    fct = len(fnames)

                    if fct > 50:
                        print 'There are %d total files\n' % fct
                        print 'Too many files to make master\n'
                        pass
                    elif fct > 35 and fct <= 50 and bns == '1X1':
                        print 'There are %d total files\n' % fct
                        sub_frames = [fnames[x:x + 15] for x in xrange(0, len(fnames), 15)]
                        print 'Creating sub_master_dark: binning is {}, date is {}'.format(bns, tagdate)
                        if not os.path.exists(darkdir):
                            os.makedirs(darkdir)
                        for i in range(3):
                            master_dark(files=sub_frames[i], outdir=darkdir,
                                        tag=tagdate + '_' + bns + '_' + exp + '_%d' % (i + 1))
                    elif fct >= 20 and fct <= 35 and bns == '1X1':
                        print 'There are %d total files\n' % fct
                        sub_frames = [fnames[x * (fct / 2):(x + 1) * (fct / 2)] for x in
                                      range((len(fnames) + (fct / 2) - 1) // (fct / 2))]
                        print 'Creating sub_master_dark: binning is {}, date is {}'.format(bns, tagdate)
                        if not os.path.exists(darkdir):
                            os.makedirs(darkdir)
                        for i in range(2):
                            master_dark(files=sub_frames[i], outdir=darkdir,
                                        tag=tagdate + '_' + bns + '_' + exp + '_%d' % (i + 1))
                    elif fct > 1:
                        if not os.path.exists(darkdir):
                            os.makedirs(darkdir)
                        print 'Creating master_dark: binning is {}, date is {}'.format(bns, tagdate)
                        master_dark(files=fnames, outdir=darkdir, tag=tagdate + '_' + bns + '_' + exp)
                    else:
                        print 'No Dark Frames found'
        del dark_files

    # if flat:
    #
    #     all_files = pd.DataFrame()
    #     all_files = pd.read_pickle(mstdir + 'all_files.pkl')
    #
    #     flat_files = all_files.where(all_files['type'] == 'Flat Field')
    #     flat_files.dropna(how='all', inplace=True)
    #
    #     del all_files
    #
    #     dates = flat_files['JD'].tolist()
    #     dates = [int(x) for x in dates]
    #     dates = np.unique(dates)
    #
    #     print '\nParsing existing data frames to determine new Flat data\n'
    #     donedates = []
    #     for date in tqdm(dates):
    #         year, month, day, sec = jd2gcal(sjd, (date - sjd))
    #         year, month, day = str(year), str(month).zfill(2), str(day).zfill(2)
    #         tagdate = year + month + day
    #         for bns in bins:
    #             flatpath = '{}{}\\Flat\\{}'.format(mstdir, bns, tagdate)
    #             if os.path.exists(flatpath):
    #                 donedates.append(date)
    #
    #     newdates = [x for x in dates if x not in donedates]
    #
    #     if clobber:
    #         newdates = dates
    #
    #     if len(newdates) == 0:
    #         print '\nDid not find any new dates with Flat Frames\n'
    #     else:
    #         print '\nBeginning master flat creation for new dates\n'
    #     for date in tqdm(newdates):
    #         year, month, day, sec = jd2gcal(sjd, (date - sjd))
    #         year, month, day = str(year), str(month).zfill(2), str(day).zfill(2)
    #         tagdate = year + month + day
    #         for bns in bins:
    #             for band in filters:
    #                 flatdir = mstdir + bns + '\\Flat\\' + tagdate + '\\'
    #                 # flatdir = mstdir+bns+'/Flat/'+tagdate+'/'
    #
    #                 filter1 = flat_files['binning'] == bns
    #                 filter2 = flat_files['JD'] == date
    #                 filter3 = flat_files['filter'] == band
    #
    #                 files = flat_files.where(filter1 & filter2 & filter3)
    #                 files.dropna(how='all', inplace=True)
    #
    #                 fnames = files['name'].tolist()
    #                 fct = len(fnames)
    #
    #                 if fct > 50:
    #                     print 'There are %d total files\n' % fct
    #                     print 'Too many files to make master\n'
    #                     pass
    #                 elif fct > 35 and fct <= 50 and bns == '1X1':
    #                     print 'There are %d total files\n' % fct
    #                     sub_frames = [fnames[x:x + 15] for x in xrange(0, len(fnames), 15)]
    #                     print 'Creating sub_master_flat: binning is {}, date is {}'.format(bns, tagdate)
    #                     if not os.path.exists(flatdir):
    #                         os.makedirs(flatdir)
    #                     for i in range(3):
    #                         master_flat(files=sub_frames[i], outdir=flatdir, tag=tagdate + '_' + bns + '_%d' % (i + 1),
    #                                     band=band)
    #                 elif fct >= 20 and fct <= 35 and bns == '1X1':
    #                     print 'There are %d total files\n' % fct
    #                     sub_frames = [fnames[x * (fct / 2):(x + 1) * (fct / 2)] for x in
    #                                   range((len(fnames) + (fct / 2) - 1) // (fct / 2))]
    #                     print 'Creating sub_master_flat: binning is {}, date is {}'.format(bns, tagdate)
    #                     if not os.path.exists(flatdir):
    #                         os.makedirs(flatdir)
    #                     for i in range(2):
    #                         master_flat(files=sub_frames[i], outdir=flatdir, tag=tagdate + '_' + bns + '_%d' % (i + 1),
    #                                     band=band)
    #                 elif fct > 1:
    #                     if not os.path.exists(flatdir):
    #                         os.makedirs(flatdir)
    #                     print '\nCreating master_flat: binning is {}, date is {}'.format(bns, tagdate)
    #                     master_flat(files=fnames, outdir=flatdir, tag=tagdate + '_' + bns, band=band)
    #                 else:
    #                     print 'No Flat Frames found'
    #     del flat_files


# ----------------------------------------------------------------------#
# master_bias:
# ----------------------------------------------------------------------#
def master_bias(files, write=True, outdir='/', readnoise=False, clobber=False, verbose=True,
                float32=True, tag='', median=True):
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
    name = outdir + 'master_bias_' + tag + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master bias already exists!")
        master_bias = fits.getdata(name, 0, header=False)
        return master_bias

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))
    temps = []

    # Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                          str(i + 1), str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        temps.append(header["CCD-TEMP"])
        stack[i, :, :] = image

    # Calculate read noise directly from bias frames if prompted
    if readnoise:
        rn = np.zeros((ysz, xsz))
        print("Starting readnoise calculation")
        pbar = tqdm(desc='Calculating readnoise', total=ysz, unit='rows')
        for i in np.arange(ysz):
            for j in np.arange(xsz):
                rn[i, j] = rb.std(stack[:, i, j])
            pbar.update(1)

        # Make a nice plot (after all that hard work)
        aspect = np.float(xsz) / np.float(ysz)
        plt.figure(39, figsize=(5 * aspect * 1.2, 5))
        plt.clf()
        sig = rb.std(rn)
        med = np.median(rn)
        mean = np.mean(rn)
        vmin = med - 2 * sig
        vmax = med + 2 * sig
        plt.imshow(rn, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.annotate(r'$\bar{\sigma}$ = %.2f cts' % mean, [0.95, 0.87], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.SimpleLineShadow(linewidth=3,foreground="w")])
        plt.annotate(r'med($\sigma$) = %.2f cts' % med, [0.95, 0.8], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.annotate(r'$\sigma_\sigma$ = %.2f cts' % sig,
                     [0.95, 0.73], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.title("Read Noise")
        plt.xlabel("pixel number")
        plt.ylabel("pixel number")

        if write:
            plt.savefig(outdir + 'readnoise' + tag + '.png', dpi=300)

    # Calculate master bias frame by median filter
    print('Calculating median of stacked frames...')
    if median:
        master_bias = np.median(stack, axis=0)
    else:
        master_bias = np.mean(stack, axis=0)

    # Make a plot
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(38, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    sig = rb.std(master_bias)
    med = np.median(master_bias)
    vmin = med - 2 * sig
    vmax = med + 2 * sig
    plt.imshow(master_bias, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.annotate('Bias Level = %.2f cts' % med, [0.95, 0.87], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large', color='k')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts' % sig, [0.95, 0.8], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.95, 0.73], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Bias")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out bias, readnoise and plot
    if write:
        name = outdir + 'master_bias_' + tag
        plt.savefig(name + '.png', dpi=300)

        hout = fits.Header()
        hout['CCDTEMP'] = (np.median(temps), "Median CCD temperature")
        hout["TEMPSIG"] = (np.std(temps), "CCD temperature RMS")
        hout["BIAS"] = (med, "Median bias level (cts)")
        hout["BIASSIG"] = (sig, "Bias RMS (cts)")
        if len(glob.glob(name + '.fits')) == 1:
            os.system('rm ' + name + '.fits')
        if float32:
            fits.writeto(name + '.fits', np.float32(master_bias), hout)
        else:
            fits.writeto(name + '.fits', master_bias, hout)

        if readnoise:
            name = outdir + 'readnoise' + tag
            if len(glob.glob(name + '.fits')) == 1:
                os.system('rm ' + name + '.fits')
            if float32:
                fits.writeto(name + '.fits', np.float32(rn), hout)
            else:
                fits.writeto(name + '.fits', rn, hout)

    return master_bias


# ----------------------------------------------------------------------#
# master_dark:
# ----------------------------------------------------------------------#
def master_dark(files, bias=None, write=True, outdir='/', clobber=False, float32=True, tag='',
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
    name = outdir + 'master_dark_' + tag + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master dark already exists!")
        master_dark = fits.getdata(name, 0, header=False)
        return master_dark

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))
    temps = []
    exps = []

    # Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                          str(i + 1), str(fct))
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
        stack[i, :, :] = image

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
        master_dark = np.median(stack, axis=0)
    else:
        master_dark = np.mean(stack, axis=0)

    # Make a plot
    sig = rb.std(master_dark)
    med = np.median(master_dark)
    vmin = med - 2 * sig
    vmax = med + 2 * sig
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(37, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    plt.imshow(master_dark, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.annotate('Dark Current = %.2f cts/sec' % med, [0.72, 0.8], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts/sec' % sig, [0.72, 0.75], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.72, 0.7], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Dark")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out plot and master dark array
    if write:
        name = outdir + 'master_dark_' + tag

        plt.savefig(name + '.png', dpi=300)

        hout = fits.Header()
        hout["TEMPMAX"] = (tmax, "Maximum CCD temperature")
        hout["TEMPMIN"] = (tmin, "Minimum CCD temperature")
        hout["TEMPMED"] = (tmed, "Median CCD temperature")
        hout["TEMPMN"] = (tmean, "Mean CCD temperature")
        hout["TEMPSIG"] = (tsig, "CCD temperature RMS")
        hout["EXPMAX"] = (expmax, "Maximum exposure time")
        hout["EXPMIN"] = (expmin, "Minimum exposure time")
        hout["DARKCNT"] = (med, "Median dark current (cts/sec)")
        hout["DARKSIG"] = (sig, "Dark current RMS (cts/sec)")
        if len(glob.glob(name)) == 1:
            os.system('rm ' + name + '.fits')
        if float32:
            fits.writeto(name + '.fits', np.float32(master_dark), hout)
        else:
            fits.writeto(name + '.fits', master_dark, hout)

    return master_dark


# ----------------------------------------------------------------------#
# master_flat:
# ----------------------------------------------------------------------#
def master_flat(files, bias=None, dark=None, write=True, outdir='/', band='V',
                tag='', clobber=False, stretch=3, float32=True, median=True):
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
    name = outdir + 'master_flat_' + tag + '_' + band + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master flat already exists!")
        master_flat = fits.getdata(name, 0, header=False)
        return master_flat

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    filter = header["filter"]

    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))

    # Load stack array and get CCD temperatures
    meds = []
    for i in np.arange(fct):
        output = '\nReading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                          str(i + 1), str(fct))
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
        stack[i, :, :] = image / np.median(image)

    # Obtain statistics for the master dark image header
    med = np.median(meds)
    sig = np.std(meds)

    # Create master flat by median filter
    if median:
        master_flat = np.median(stack, axis=0)
    else:
        master_flat = np.mean(stack, axis=0)

    # Make a plot
    sig = rb.std(master_flat)
    med = np.median(master_flat)
    vmin = med - stretch * sig
    vmax = med + stretch * sig
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(40, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    plt.imshow(master_flat, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.title("Master Flat")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out plot and master flat array
    if write:
        plt.savefig(outdir + 'master_flat' + tag + '_' + band + '.png', dpi=300)
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

        if len(glob.glob(outdir + 'master_flat' + tag + '.fits')) == 1:
            os.system('rm ' + outdir + 'master_flat' + tag + '.fits')
        if float32:
            fits.writeto(outdir + 'master_flat_' + tag + '_' + band + '.fits', np.float32(master_flat), hout)
        else:
            fits.writeto(outdir + 'master_flat_' + tag + '_' + band + '.fits', master_flat, hout)

    return master_flat
