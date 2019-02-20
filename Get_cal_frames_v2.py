import shutil, os
import pandas as pd
import numpy as np
from jdcal import jd2gcal, gcal2jd
import calgui as cg

types, outpath, binning, temp, dates, exposures, filters = cg.gui()

userJD = []
for date in dates:
    jd1, jd2 = gcal2jd(int(date[0:4]), int(date[4:6]), int(date[6:]))
    userJD.append(jd1 + jd2)


sjd = 2400000.5
mstdir = 'Z:\Calibration Master Frames\\'

if not os.path.exists(outpath):
    os.mkdir(outpath)

cal_files = pd.DataFrame()
cal_files = pd.read_pickle(mstdir+'all_files.pkl')

binfiles = cal_files.where(cal_files['binning'] == binning)
binfiles.dropna(how='all', inplace=True)

tempfiles = binfiles.where(binfiles['temp'] == temp)
tempfiles.dropna(how='all', inplace=True)


if 'Bias' in types:
    print 'Copying master bias files'
    bias_files = tempfiles.where(tempfiles['type'] == 'Bias Frame')
    bias_files.dropna(how='all', inplace=True)
    biasJD = bias_files['JD']
    for date in userJD:
        year, month, day, sec = jd2gcal(sjd, (date - sjd))
        usertag = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
        biasdate = min(biasJD, key=lambda x: abs(x - date))
        year, month, day, sec = jd2gcal(sjd, (biasdate - sjd))
        tagdate = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
        fname = '{}{}\\Bias\\{}\\master_bias_{}_{}.fits'.format(mstdir, binning, tagdate, tagdate, binning)
        outname = fname.replace(mstdir + binning + '\\Bias\\' + tagdate,\
                                outpath + usertag + '\\')
        if not os.path.exists(outpath + usertag + '\\'):
            os.mkdir(outpath + usertag + '\\')
        try:
            shutil.copy(fname, outname)
        except:
            shutil.copy(fname.replace('.fits', '_1.fits'), outname)


if 'Dark' in types:
    print 'Copying master dark files'
    dark_files = tempfiles.where(tempfiles['type'] == 'Dark Frame')
    dark_files.dropna(how='all', inplace=True)
    for exp in exposures:
        darkexp = dark_files.where(dark_files['exp'] == exp)
        darkexp.dropna(how='all', inplace=True)
        darkJD = darkexp['JD']
        for date in userJD:
            year, month, day, sec = jd2gcal(sjd, (date - sjd))
            usertag = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
            darkdate = min(darkJD, key=lambda x: abs(x - date))
            year, month, day, sec = jd2gcal(sjd, (darkdate - sjd))
            tagdate = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
            fname = '{}{}\\Dark\\{}\\master_dark_{}_{}_{}.fits'.format(mstdir, binning, tagdate, tagdate, binning, exp)
            outname = fname.replace(mstdir + binning + '\\Dark\\' + tagdate,\
                                    outpath + usertag + '\\')
            if not os.path.exists(outpath + usertag + '\\'):
                os.mkdir(outpath + usertag + '\\')
            try:
                shutil.copy(fname, outname)
            except:
                shutil.copy(fname.replace('.fits', '_1.fits'), outname)


if 'Flat' in types:
    print 'Copying master flat files'
    flat_files = tempfiles.where(tempfiles['type'] == 'Flat Field')
    flat_files.dropna(how='all', inplace=True)
    for band in filters:
        flatband = flat_files.where(flat_files['filter'] == band)
        flatband.dropna(how='all', inplace=True)
        flatJD = flatband['JD']
        for date in userJD:
            year, month, day, sec = jd2gcal(sjd, (date - sjd))
            usertag = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
            flatdate = min(flatJD, key=lambda x: abs(x - date))
            year, month, day, sec = jd2gcal(sjd, (flatdate - sjd))
            tagdate = '{}{}{}'.format(str(year), str(month).zfill(2), str(day).zfill(2))
            fname = '{}{}\\Flat\\{}\\master_flat_{}_{}_{}.fits'.format(mstdir, binning, tagdate, tagdate, binning, band)
            outname = fname.replace(mstdir + binning + '\\Flat\\' + tagdate,\
                                    outpath + usertag + '\\')
            if not os.path.exists(outpath + usertag + '\\'):
                os.mkdir(outpath + usertag + '\\')
            try:
                shutil.copy(fname, outname)
            except:
                shutil.copy(fname.replace('.fits', '_1.fits'), outname)