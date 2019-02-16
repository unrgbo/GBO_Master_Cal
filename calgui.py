import PySimpleGUI27 as sg
import pandas as pd
import numpy as np

mstdir = '/Users/jakrin/Dropbox/JFausett/'
#mstdir = 'Z:\Calibration Master Frames\\'
outdir = '/Users/jfausett/Dropbox/JFausett/test_out/'
done_files = pd.read_pickle(mstdir + 'master_files.pkl')

bins = np.sort(done_files['binning'].unique())
bins = np.sort(bins)


dark_files = done_files.where(done_files['type'] == 'Dark Frame')

exposures = []
for bn in bins:
    binfiles = dark_files.where(done_files['binning'] == bn)
    binfiles.dropna(how='all', inplace=True)
    expose = binfiles['exp'].unique()
    for x in expose:
        exposures.append(x)

exposures = pd.Series(exposures)
exposures = np.sort(exposures.unique())


print exposures

def gui():


    #sg.SetOptions(text_justification='right')


    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Enter output path for your files', font=('Helvetica', 12))],
              [sg.InputText('C:\Users\user\Dropbox\Master_Finder_Output', size=(40, 1)), sg.Text('Leave as is for default')],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()
    outpath = values[0]

    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
              [sg.Text(outpath, font=('Helvetica', 14))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
              [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
               sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
              [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
               sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()
    dates = values

    print dates

    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
              [sg.Text(outpath, font=('Helvetica', 14))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
              [sg.Text(str(dates), font=('Helvetica', 12))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Binning', font=('Helvetica', 15), justification='left')],
              [sg.Checkbox('_1X1_', size=(12, 1)), sg.Checkbox('2X2', size=(12, 1)),
               sg.Checkbox('3X3', size=(12, 1)), sg.Checkbox('4X4', size=(12, 1))],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()

    bins = values
    print bins

    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
              [sg.Text(outpath, font=('Helvetica', 14))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
              [sg.Text(str(dates), font=('Helvetica', 12))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Binning', font=('Helvetica', 14), justification='left')],
              [sg.Text(bins, font=('Helvetica', 12))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Checkbox('005', size=(12, 1)), sg.Checkbox('010', size=(12, 1)), sg.Checkbox('015', size=(12, 1)),
               sg.Checkbox('030', size=(12, 1)), sg.Checkbox('040', size=(12, 1))],
              [sg.Checkbox('045', size=(12, 1)), sg.Checkbox('050', size=(12, 1)), sg.Checkbox('060', size=(12, 1)),
               sg.Checkbox('090', size=(12, 1)), sg.Checkbox('120', size=(12, 1))],
              [sg.Checkbox('160', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('200', size=(12, 1)),
               sg.Checkbox('220', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
              [sg.Checkbox('280', size=(12, 1)), sg.Checkbox('300', size=(12, 1)), sg.Checkbox('360', size=(12, 1)),
               sg.Checkbox('600', size=(12, 1))],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()


    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
              [sg.Text(outpath, font=('Helvetica', 14))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
              [sg.Text(str(dates), font=('Helvetica', 12))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Binning', font=('Helvetica', 14), justification='left')],
              [sg.Text(bins, font=('Helvetica', 12))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Checkbox('005', size=(12, 1)), sg.Checkbox('010', size=(12, 1)), sg.Checkbox('015', size=(12, 1)),\
              sg.Checkbox('030', size=(12, 1)), sg.Checkbox('040', size=(12, 1))],
              [sg.Checkbox('045', size=(12, 1)), sg.Checkbox('050', size=(12, 1)), sg.Checkbox('060', size=(12, 1)),\
              sg.Checkbox('090', size=(12, 1)), sg.Checkbox('120', size=(12, 1))],
              [sg.Checkbox('160', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('200', size=(12, 1)),\
               sg.Checkbox('220', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
              [sg.Checkbox('280', size=(12, 1)), sg.Checkbox('300', size=(12, 1)), sg.Checkbox('360', size=(12, 1)),\
               sg.Checkbox('600', size=(12, 1))],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()


    # layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
    #           [sg.Text('Enter output path for your files', font=('Helvetica', 12))],
    #           [sg.InputText('C:\Users\user\Dropbox\Master_Finder_Output', size=(40, 1)), sg.Text('Leave as is for default')],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 12))],
    #           [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
    #            sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
    #           [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
    #            sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Binning', font=('Helvetica', 15), justification='left')],
    #           [sg.Checkbox(sg.T('', key='_1X1_'), size=(12, 1)), sg.Checkbox('2X2', size=(12, 1)),
    #            sg.Checkbox('3X3', size=(12, 1)), sg.Checkbox('4X4', size=(12, 1))],
    #           [sg.Text('_'  * 100, size=(65, 1))],
    #           [sg.Text('Exposures', font=('Helvetica', 15), justification='left')],
    #           [sg.Checkbox('030', size=(12, 1)), sg.Checkbox('060', size=(12, 1)), sg.Checkbox('090', size=(12, 1))],
    #           [sg.Checkbox('120', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
    #           [sg.Checkbox('300', size=(12, 1)), sg.Checkbox('600', size=(12, 1))],
    #           [sg.Text('_'  * 100, size=(65, 1))],
    #           [sg.Text('Filters', font=('Helvetica', 15), justification='left')],
    #           [sg.Checkbox('Luminance', size=(12, 1)), sg.Checkbox('Red', size=(12, 1)), sg.Checkbox('Green', size=(12, 1)),
    #            sg.Checkbox('Blue', size=(12, 1)), sg.Checkbox('Grating', size=(12, 1))],
    #           [sg.Checkbox('B', size=(12, 1)), sg.Checkbox('V', size=(12, 1)), sg.Checkbox('R', size=(12, 1)),
    #            sg.Checkbox('I', size=(12, 1)), sg.Checkbox('Y', size=(12, 1))],
    #           [sg.Checkbox("g'", size=(12, 1)), sg.Checkbox("r'", size=(12, 1)), sg.Checkbox("i'", size=(12, 1)),
    #            sg.Checkbox("z'", size=(12, 1)), sg.Checkbox('OIII', size=(12, 1))],
    #           [sg.Checkbox('Ha', size=(12, 1)), sg.Checkbox('SII', size=(12, 1))],
    #           [sg.Submit(), sg.Button('Exit')]]
    #
    # window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)
    #


    return outpath, dates, bins
