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
temps = done_files['temp'].unique()

exposures = pd.Series(exposures)
exposures = np.sort(exposures.unique())


print exposures

dates = []

def gui():


    #sg.SetOptions(text_justification='right')


    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Enter output path for your files', font=('Helvetica', 12))],
              [sg.InputText('C:\Users\user\Dropbox\Master_Finder_Output', size=(40, 1), key='path'),
               sg.Text('Leave as is for default')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Which master calibration frames are you looking for?', font=('Helvetica', 14))],
              [sg.Checkbox('Bias', default=True, size=(12, 1), key='Bias'),
               sg.Checkbox('Dark', default=True, size=(12, 1), key='Dark'),
               sg.Checkbox('Flat', default=True, size=(12, 1), key='Flat')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14)),
               sg.T('Dates')],
              [sg.InputText(None, size=(12, 1), key='datein'), sg.T('', key='dateout')],
              [sg.Button('Add')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Select the Binning', font=('Helvetica', 14))],
              [sg.Radio('1X1', "RADIO1", default=True, key='1X1'), sg.Radio('2X2', "RADIO1", key='2X2'),
               sg.Radio('3X3', "RADIO1", key='3X3'), sg.Radio('4X4', "RADIO1", key='4X4')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Select the Temp', font=('Helvetica', 14))],
              [sg.Radio('-25', "RADIO2", default=True, key='-25'), sg.Radio('-30', "RADIO2", key='-30'),
               sg.Radio('-35', "RADIO2", key='-35')],
              [sg.Submit(), sg.Button('Exit')]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    while True:  # Event Loop
        event, values = window.Read()
        print event, values
        if event is None or event == 'Exit':
            break
        if event == 'Add':
            dates.append(values['datein'])
            window.Element('dateout').Update(values['datein'])


    print values

    binning = [x for x in bins if values[x] == True]
    temp = [x for x in temps if values[x] == True]

    print binning, temp

    # outpath = values[0]
    # types = values[1:4]
    # dates = values[4:12]
    # binning = values[12:16]
    # temp = values[16:]
    #
    # bins = [x for x in b]
    #
    # print event, outpath, types, dates, binning, temp

    # layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
    #           [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
    #           [sg.Text(outpath, font=('Helvetica', 14))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
    #           [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
    #            sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
    #           [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
    #            sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Submit(), sg.Cancel()]]
    #
    # window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)
    #
    # event, values = window.Read()
    # dates = values
    #
    # print dates
    #
    # layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
    #           [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
    #           [sg.Text(outpath, font=('Helvetica', 14))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
    #           [sg.Text(str(dates), font=('Helvetica', 12))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Binning', font=('Helvetica', 15), justification='left')],
    #           [sg.Checkbox('_1X1_', size=(12, 1)), sg.Checkbox('2X2', size=(12, 1)),
    #            sg.Checkbox('3X3', size=(12, 1)), sg.Checkbox('4X4', size=(12, 1))],
    #           [sg.Submit(), sg.Cancel()]]
    #
    # window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)
    #
    # event, values = window.Read()
    #
    # bins = values
    # print bins
    #
    # layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
    #           [sg.Text('Master Calibration frames will be copied to:', font=('Helvetica', 14))],
    #           [sg.Text(outpath, font=('Helvetica', 14))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 14))],
    #           [sg.Text(str(dates), font=('Helvetica', 12))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Text('Binning', font=('Helvetica', 14), justification='left')],
    #           [sg.Text(bins, font=('Helvetica', 12))],
    #           [sg.Text('_' * 100, size=(65, 1))],
    #           [sg.Checkbox('005', size=(12, 1)), sg.Checkbox('010', size=(12, 1)), sg.Checkbox('015', size=(12, 1)),
    #            sg.Checkbox('030', size=(12, 1)), sg.Checkbox('040', size=(12, 1))],
    #           [sg.Checkbox('045', size=(12, 1)), sg.Checkbox('050', size=(12, 1)), sg.Checkbox('060', size=(12, 1)),
    #            sg.Checkbox('090', size=(12, 1)), sg.Checkbox('120', size=(12, 1))],
    #           [sg.Checkbox('160', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('200', size=(12, 1)),
    #            sg.Checkbox('220', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
    #           [sg.Checkbox('280', size=(12, 1)), sg.Checkbox('300', size=(12, 1)), sg.Checkbox('360', size=(12, 1)),
    #            sg.Checkbox('600', size=(12, 1))],
    #           [sg.Submit(), sg.Cancel()]]
    #
    # window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)
    #
    # event, values = window.Read()
    #
    #
    return dates
