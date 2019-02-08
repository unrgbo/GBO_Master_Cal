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
                     '\nexample, C:\Dropbox\User\Object')

prefix = raw_input(('Type the prefix for your image file names\n' +
                    '\nexample, NGC1976'))

suffix = raw_input(('Type the suffix for your image file names\n' +
                    '\nexample, .fts'))

msdir = 'Z:\Calibration Master Frames\\'

caltypes = ['Bias', 'Dark', 'Flat']

done_files = pd.read_pickle(mstdir + 'all_files.pkl')

datafiles = glob.glob(datapath + '\\{}*{}'.format(prefix, suffix))



dates = []
i = 0
for type in caltypes:
    type = done_files.where(all_files['type'] == '{} Frame'.format(type))
    type.dropna(how='all', inplace=True)
    dates[i] = np.unique(type.loc[:, 'JD'])
    dates[i] = [int(x) for x in dates]
    i += 1

biasdates = dates[0]
darkdates = dates[1]
flatdates = dates[2]


try:
    caldates = []
    print '\nChecking given path\n'
    data_files = []
    for path, subdirs, files in os.walk(datapath):
        for name in files:
            if name.endswith(suffix) and name.startswith(prefix):
                data_files.append(os.path.join(path, name))
    print '\nFound {} files in {}'.format(len(list_of_files), datapath)


    cal_list

    for root, dirs, files in os.walk(datapath):
        for file in files:
            if file.endswith(".fits"):
                caldates.append(root.split('\\')[-1])


    types = []
    names = []
    temp = []
    binning = []
    JD = []
    exposure = []
    bands = []
    bad_files = []
    no_temp = []



except:
    print '\nCould not find {}'.format(datapath)


