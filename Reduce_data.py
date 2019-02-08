import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import robust as rb
#from length import *
import math
from tqdm import tqdm
from jdcal import gcal2jd, jd2gcal

datapath = raw_input('Type the exact path to the directory containing your data:\n' +
                     '\nexample, C:\Dropbox\User\Object\n')

prefix = raw_input(('Type the prefix for your image file names\n' +
                    '\nexample, NGC1976\n'))

suffix = raw_input(('Type the suffix for your image file names\n' +
                    '\nexample, .fts\n'))

sjd = 2400000.5

mstdir = 'Z:\Calibration Master Frames\\'

caltypes = ['Bias', 'Dark', 'Flat']

datafiles = glob.glob(datapath + '\\{}*{}'.format(prefix, suffix))

names = []
binning = []
JD = []
exposure = []
bands = []
bad_files = []
for fname in tqdm(datafiles):
    try:
        with fits.open(fname) as hdu:
            names.append(fname)
            binning.append(str(hdu[0].header['XBINNING']) + 'X' + str(hdu[0].header['YBINNING']))
            JD.append(int(hdu[0].header['JD']))
            exposure.append(str(int(hdu[0].header['EXPOSURE'])))
            bands.append(hdu[0].header['FILTER'])
    except:
        print '\nThe file: {}\nIs missing key header information'.format(fname)
        bad_files.append(fname)

print '\nHere is a list of the files that raised errors:\n'
print bad_files

d = {'name': names, 'binning': binning, 'JD': JD, 'exp': exposure, 'filter': bands}
dataframe = pd.DataFrame(data=d)

exposures = np.unique(dataframe.loc[:,'exp'])
bands = np.unique(dataframe.loc[:,'filter'])
bins = np.unique(dataframe.loc[:,'binning'])

print 'Found data with {} binning, {} filters, and {} exposure times'.format(bins, bands, exposures)

print 'Beginning reduction'

cal_files = pd.read_pickle(mstdir + 'all_files.pkl')

for fname in datafiles:
    with fits.open(fname) as hdu:
        data, hdr = hdu[0].data, hdu[0].header

        date = int(hdr['JD'])
        exp = str(int(hdr['EXPOSURE']))
        band = hdr['FILTER']
        binning = str(hdr['XBINNING']) + 'X' + str(hdr['YBINNING'])

        cal_files = cal_files.where(cal_files['binning'] == binning)

        print 'Parsing Bias Data'
        bias_files = cal_files.where(cal_files['type'] == 'Bias Frame')
        bias_files.dropna(how='all', inplace=True)
        biasdates = bias_files['JD'].tolist()
        biasdates = np.unique(biasdates)
        biasdate = min(biasdates, key=lambda x: abs(int(x) - date))
        year, month, day, sec = jd2gcal(sjd, (biasdate - sjd))
        year, month, day = str(year), str(month), str(day)
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        biastag = year + month + day
        biaspath = '{}{}\\Bias\\{}\\master_bias_{}_{}.fits'.format(mstdir, binning, biastag, biastag, binning)
        bias = fits.getdata(biaspath)
        print 'Using {}'.format(biaspath)

        filter1 = cal_files['type'] == 'Dark Frame'
        filter2 = cal_files['exp'] == exp
        dark_files = cal_files.where(filter1 & filter2)
        dark_files.dropna(how='all', inplace=True)
        darkdates = dark_files['JD'].tolist()
        darkdates = np.unique(darkdates)
        darkdate = min(darkdates, key=lambda x: abs(int(x) - date))
        year, month, day, sec = jd2gcal(sjd, (darkdate - sjd))
        year, month, day = str(year), str(month), str(day)
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        darktag = year + month + day
        darkpath = '{}{}\\Dark\\{}\\master_dark_{}_{}_{}.fits'.format(mstdir, binning, darktag, darktag, binning, exp)
        dark = fits.getdata(darkpath)*int(exp)
        print 'Using {} '.format(darkpath)

        filter1 = cal_files['type'] == 'Flat Field'
        filter2 = cal_files['filter'] == band
        flat_files = cal_files.where(filter1 & filter2)
        flat_files.dropna(how='all', inplace=True)
        flatdates = flat_files['JD'].tolist()
        flatdates = np.unique(flatdates)
        flatdate = min(flatdates, key=lambda x: abs(int(x) - date))
        year, month, day, sec = jd2gcal(sjd, (flatdate - sjd))
        year, month, day = str(year), str(month), str(day)
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        flattag = year + month + day
        flatpath = '{}{}\\Flat\\{}\\master_flat_{}_{}_{}.fits'.format(mstdir, binning, flattag, flattag, binning, band)
        flat = fits.getdata(flatpath)
        print 'Using {} '.format(flatpath)



