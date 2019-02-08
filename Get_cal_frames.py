import shutil, os
import pandas as pd
import numpy as np
from jdcal import jd2gcal, gcal2jd

dates = input('What dates do you need calibration frames for\n'
              '*Must enter as a list of strings*\n'
              'Example: ["20180901","20190105",...]\n')

bins = input('What binning do you need\n'
             '*Must enter as a list of strings*\n'
             'Example: ["1X1"]\n')

filters = input('What filters do you need\n'
                '*Must enter as a list of strings*\n'
                'Example: ["V","Red","Luminance",...]\n')

exposures = input('What exposure times do you need\n'
                  '*Must enter as a list of strings*\n'
                  'Example: ["030", "300"]\n')

print dates, bins, filters, exposures
print type(dates), type(bins), type(filters), type(exposures)

sjd = 2400000.5
mstdir = 'Z:\Calibration Master Frames\\'
outpath = 'C:\Users\user\Desktop\Test_out'

if not os.path.exists(outpath):
    os.mkdir(outpath)

cal_files = pd.DataFrame()
cal_files = pd.read_pickle(mstdir+'all_files.pkl')


for bns in bins:
    for date in dates:
        filter1 = cal_files['type'] == 'Bias Frame'
        filter2 = cal_files['binning'] == bns
        bias_files = cal_files.where(filter1 & filter2)
        bias_files.dropna(how='all', inplace=True)
        year, month, day = int(date[0:4]), int(date[4:6]), int(date[6:])
        jd1, jd2 = gcal2jd(year, month, day)
        jddate = jd1 + jd2
        bdates = bias_files['JD'].tolist()
        bdates = np.unique(bdates)
        biasdate = min(bdates, key=lambda x: abs(int(x) - jddate))
        year, month, day, sec = jd2gcal(sjd, (biasdate - sjd))
        year, month, day = str(year), str(month), str(day)
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        tagdate = year + month + day
        try:
            sourcepath = '{}{}\Bias\{}\master_bias_{}_{}.fits'.format(mstdir, bns, tagdate, tagdate, bns)
            shutil.copy(sourcepath, outpath)
        except:
            sourcepath = '{}{}\Bias\{}\master_bias_{}_{}_1.fits'.format(mstdir, bns, tagdate, tagdate, bns)
            shutil.copy(sourcepath, outpath)

for bns in bins:
    for exp in exposures:
        for date in dates:
            filter1 = cal_files['type'] == 'Dark Frame'
            filter2 = cal_files['binning'] == bns
            filter3 = cal_files['exp'] == exp
            dark_files = cal_files.where(filter1 & filter2 & filter3)
            dark_files.dropna(how='all', inplace=True)
            year, month, day = int(date[0:4]), int(date[4:6]), int(date[6:])
            jd1, jd2 = gcal2jd(year, month, day)
            jddate = jd1 + jd2
            ddates = dark_files['JD'].tolist()
            ddates = np.unique(ddates)
            darkdate = min(ddates, key=lambda x: abs(int(x) - jddate))
            year, month, day, sec = jd2gcal(sjd, (darkdate - sjd))
            year, month, day = str(year), str(month), str(day)
            if len(month) == 1:
                month = '0' + month
            if len(day) == 1:
                day = '0' + day
            tagdate = year + month + day
            try:
                sourcepath = '{}{}\Dark\{}\master_dark_{}_{}_{}.fits'.format(mstdir, bns, tagdate, tagdate, bns, exp)
                shutil.copy(sourcepath, outpath)
            except:
                sourcepath = '{}{}\Dark\{}\master_dark_{}_{}_{}_1.fits'.format(mstdir, bns, tagdate, tagdate, bns, exp)
                shutil.copy(sourcepath, outpath)

for bns in bins:
    for band in filters:
        for date in dates:
            filter1 = cal_files['type'] == 'Flat Field'
            filter2 = cal_files['binning'] == bns
            filter3 = cal_files['filter'] == band
            flat_files = cal_files.where(filter1 & filter2 & filter3)
            flat_files.dropna(how='all', inplace=True)
            year, month, day = int(date[0:4]), int(date[4:6]), int(date[6:])
            jd1, jd2 = gcal2jd(year, month, day)
            jddate = jd1 + jd2
            fdates = flat_files['JD'].tolist()
            fdates = np.unique(fdates)
            flatdate = min(fdates, key=lambda x: abs(int(x) - jddate))
            year, month, day, sec = jd2gcal(sjd, (flatdate - sjd))
            year, month, day = str(year), str(month), str(day)
            if len(month) == 1:
                month = '0' + month
            if len(day) == 1:
                day = '0' + day
            tagdate = year + month + day
            try:
                sourcepath = '{}{}\Flat\{}\master_flat_{}_{}_{}.fits'.format(mstdir, bns, tagdate, tagdate, bns, band)
                shutil.copy(sourcepath, outpath)
            except:
                sourcepath = '{}{}\Flat\{}\master_flat_{}_{}_{}_1.fits'.format(mstdir, bns, tagdate, tagdate, bns, band)
                shutil.copy(sourcepath, outpath)