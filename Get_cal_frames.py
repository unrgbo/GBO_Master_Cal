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


bias_files = cal_files.where(cal_files['type'] == 'Bias Frame')
bias_files.dropna(how='all', inplace=True)


biasdates = []
biaspaths = []
for bns in bins:
    files = bias_files.where(bias_files['binning'] == bns)
    files.dropna(how='all', inplace=True)
    for date in dates:
        year, month, day = int(date[0:4]), int(date[4:6]), int(date[6:])
        jd1, jd2 = gcal2jd(year, month, day)
        jddate = jd1 + jd2
        bdates = files['JD'].tolist()
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