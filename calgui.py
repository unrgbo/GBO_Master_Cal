import PySimpleGUI27 as sg

# Green & tan color scheme
sg.ChangeLookAndFeel('GreenTan')

#sg.SetOptions(text_justification='right')

layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
          [sg.Text('Enter the path to copy files to', font=('Helvetica', 18))],
          [sg.InputText('Default: C:\Users\user\Desktop\Master_Finder_Output', size=(50, 1))],
          [sg.Text('_' * 100, size=(65, 1))],
          [sg.Text('Binning', font=('Helvetica', 15), justification='left')],
          [sg.Checkbox('1X1', size=(12, 1)), sg.Checkbox('3X3', size=(12, 1))],
          [sg.Checkbox('2X2', size=(12, 1)), sg.Checkbox('4X4', size=(12, 1))],
          [sg.Text('_'  * 100, size=(65, 1))],
          [sg.Text('Exposures', font=('Helvetica', 15), justification='left')],
          [sg.Checkbox('030', size=(12, 1)), sg.Checkbox('060', size=(12, 1)), sg.Checkbox('090', size=(12, 1))],
          [sg.Checkbox('120', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
          [sg.Checkbox('300', size=(12, 1)), sg.Checkbox('600', size=(12, 1))],
          [sg.Text('_'  * 100, size=(65, 1))],
          [sg.Text('Filters', font=('Helvetica', 15), justification='left')],
          [sg.Checkbox('Normalize', size=(12, 1), default=True), sg.Checkbox('Verbose', size=(20, 1))],
          [sg.Checkbox('Cluster', size=(12, 1)), sg.Checkbox('Flush Output', size=(20, 1), default=True)],
          [sg.Checkbox('Write Results', size=(12, 1)), sg.Checkbox('Keep Intermediate Data', size=(20, 1))],
          [sg.Submit(), sg.Cancel()]]

window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

event, values = window.Read()