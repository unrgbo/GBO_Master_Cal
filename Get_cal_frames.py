import shutil, os
import pandas as pd
import numpy as np
from jdcal import jd2gcal

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
dark_files = cal_files.where(cal_files['type'] == 'Dark Frame')
dark_files.dropna(how='all', inplace=True)
flat_files = cal_files.where(cal_files['type'] == 'Flat Field')
flat_files.dropna(how='all', inplace=True)

bdates = bias_files['JD'].tolist()
bdates = np.unique(bdates)
biasdates = []
for date in dates:
    biasdate = min(bdates, key=lambda x: abs(int(x) - date))
    biasdates.append(biasdate)
biaspaths = []
for date in biasdates:
    for bns in bins:
        year, month, day, sec = jd2gcal(sjd, (date - sjd))
        year, month, day = str(year), str(month), str(day)
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        tagdate = year + month + day
        sourcepath = '{}{}\\Bias\\{}\\master_bias_{}_{}.fits'.format(mstdir, bns, tagdate, tagdate, bns)
        shutil.copy(sourcepath, outpath)

darkdates = []
for exp in exposures:

    filter1 = dark_files['exp'] == exp

    files = dark_files.where(filter1)
    files.dropna(how='all', inplace=True)

    ddates = files['JD'].tolist()
    ddates = np.unique(ddates)


    fnames = files['name'].tolist()
    fct = len(fnames)