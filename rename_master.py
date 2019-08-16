import os


bins = ['1X1', '2X2', '3X3', '4X4']

for bn in bins:
    root = 'Z:\Calibration Master Frames\\{}\\'.format(bn)
    list_of_files = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            if name.endswith('V.fits'):
                list_of_files.append(os.path.join(path, name))
    for fname in list_of_files:
        os.rename(fname, fname.replace('_V.fits', '.fits'))