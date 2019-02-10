import calgui as gui
#import Get_cal_frames as gc
user_input = gui.gui()

outpath = user_input[0]
dates = user_input[1:9]
bins = user_input[9:13]
exposures = user_input[13:21]
filters = user_input[21:]

print outpath
print dates
print bins
print exposures
print filters