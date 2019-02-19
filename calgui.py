import PySimpleGUI27 as sg
import pandas as pd
import numpy as np

#  Paths to use
mstdir = '/Users/jfausett/Dropbox/JFausett/'
#mstdir = 'Z:\Calibration Master Frames\\'
outdir = 'C:\Users\user\Dropbox\Master_Finder_Output\\'

#  Import data frame to parse
done_files = pd.read_pickle(mstdir + 'master_files.pkl')

bins = done_files['binning'].unique()
temps = done_files['temp'].unique()

dark_files = done_files.where(done_files['type'] == 'Dark Frame')

flat_files = done_files.where(done_files['type'] == 'Flat Field')


def gui():

    #  Gui to get Path, Binning, and Temp from user
    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter output path for your files', font=('Helvetica', 12))],
              [sg.InputText(outdir, size=(40, 1), key='path'), sg.Text('Leave as is for default')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Select the type of master frames to find', font=('Helvetica', 14))],
              [sg.Checkbox('Bias', default=True, size=(12, 1), key='Bias'),
               sg.Checkbox('Dark', default=True, size=(12, 1), key='Dark'),
               sg.Checkbox('Flat', default=True, size=(12, 1), key='Flat')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Select the Binning', font=('Helvetica', 14))],
              [sg.Radio('1X1', "RADIO1", default=True, key='1X1'), sg.Radio('2X2', "RADIO1", key='2X2'),
               sg.Radio('3X3', "RADIO1", key='3X3'), sg.Radio('4X4', "RADIO1", key='4X4')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Select the Temp', font=('Helvetica', 14))],
              [sg.Radio('-25', "RADIO2", default=True, key='-25'), sg.Radio('-30', "RADIO2", key='-30'),
               sg.Radio('-35', "RADIO2", key='-35')],
              [sg.Submit(), sg.Button('Exit')]]

    window = sg.Window('GBO Master Calibration Finder', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()

    window.Close()

    types = ['Bias', 'Dark', 'Flat']
    types = [x for x in types if values[x] == True]
    binning = [x for x in bins if values[x] == True]
    binning = binning[0]
    temp = [x for x in temps if values[x] == True]
    temp = temp[0]

    outpath = values['path']
    if not outpath.endswith('\\'):
        outpath = outpath + '\\'

    #  Gui to get dates from the user
    dates = []
    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images', font=('Helvetica', 14))],
              [sg.T('Format: 20180715'), sg.T('Added Date', size=(19,1), justification='right')],
              [sg.InputText('', size=(12, 1), key='datein'),
              sg.T('', size=(20, 1), justification='right', key='dateout')],
              [sg.Submit('Add Date')],
              [sg.Button('Done'), sg.Button('Exit')]]

    window = sg.Window('Window Title').Layout(layout)

    while True:  # Event Loop
        event, values = window.Read()
        print event, values
        if event is None or event == 'Exit':
            break
        if event == 'Done':
            break
        if event == 'Add Date':
            dates.append(values['datein'])
            window.Element('dateout').Update(values['datein'])

    window.Close()

    #  Gui to get dark exposure times from the user
    if 'Dark' in types:
        darkbin = dark_files.where(done_files['binning'] == binning)
        darkbin.dropna(how='all', inplace=True)
        darktemp = darkbin.where(darkbin['temp'] == temp)
        darktemp.dropna(how='all', inplace=True)

        exposures = np.sort(darktemp['exp'].unique())

        layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
                  [sg.Text('Your Binning is:'), sg.Text(binning)],
                  [sg.Text('Your Temp is:', size=(12,1)), sg.Text(temp)],
                  [sg.Text('',size=(15,1)), sg.Text('Available Exposures')],
                  [sg.Text('', size=(16,1)),
                   sg.Listbox(exposures, size=(10,12), select_mode='multiple', key='exposures')],
                  [sg.Submit(), sg.Button('Exit')]]

        window = sg.Window('GBO Master Calibration Finder', font=("Helvetica", 12)).Layout(layout)

        event, values = window.Read()

        window.Close()

        exposures = values['exposures']
    else:
        exposures = []

    #  Gui to get filters from the user
    if 'Flat' in types:
        flatbin = flat_files.where(done_files['binning'] == binning)
        flatbin.dropna(how='all', inplace=True)
        flattemp = flatbin.where(flatbin['temp'] == temp)
        flattemp.dropna(how='all', inplace=True)

        filters = np.sort(flattemp['filter'].unique())

        layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
                  [sg.Text('Your Binning is:'), sg.Text(binning)],
                  [sg.Text('Your Temp is:', size=(12,1)), sg.Text(temp)],
                  [sg.Text('',size=(15,1)), sg.Text('Available Filters')],
                  [sg.Text('', size=(16,1)),
                   sg.Listbox(filters, size=(10,12), select_mode='multiple', key='filters')],
                  [sg.Submit(), sg.Button('Exit')]]

        window = sg.Window('GBO Master Calibration Finder', font=("Helvetica", 12)).Layout(layout)

        event, values = window.Read()

        window.Close()

        filters = values['filters']
    else:
        filters = []

    return types, outpath, binning, temp, dates, exposures, filters



