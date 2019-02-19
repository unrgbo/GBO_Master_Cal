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
    bias_files = tempfiles.where(tempfiles['type'] == 'Bias Frame')
    bias_files.dropna(how='all', inplace=True)
    biasJD = bias_files['JD']
    biasJD = pd.series([int(x) for x in biasJD])
    biasJD = biasJD.unique()
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
        shutil.copy(fname, outname)