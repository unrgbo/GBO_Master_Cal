import PySimpleGUI27 as sg


def gui():
    # Green & tan color scheme
    sg.ChangeLookAndFeel('GreenTan')

    #sg.SetOptions(text_justification='right')

    layout = [[sg.Text('GBO Master Calibration Finder', font=('Helvetica', 25), justification='center')],
              [sg.Text('Enter output path for your files', font=('Helvetica', 12))],
              [sg.InputText('C:\Users\user\Dropbox\Master_Finder_Output', size=(40, 1)), sg.Text('Leave as is for default')],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Enter the dates of science images: Format: 20180715', font=('Helvetica', 12))],
              [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
               sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
              [sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1)),
               sg.InputText(None, size=(12, 1)), sg.InputText(None, size=(12, 1))],
              [sg.Text('_' * 100, size=(65, 1))],
              [sg.Text('Binning', font=('Helvetica', 15), justification='left')],
              [sg.Checkbox('1X1', size=(12, 1)), sg.Checkbox('2X2', size=(12, 1))],
              [sg.Checkbox('3X3', size=(12, 1)), sg.Checkbox('4X4', size=(12, 1))],
              [sg.Text('_'  * 100, size=(65, 1))],
              [sg.Text('Exposures', font=('Helvetica', 15), justification='left')],
              [sg.Checkbox('030', size=(12, 1)), sg.Checkbox('060', size=(12, 1)), sg.Checkbox('090', size=(12, 1))],
              [sg.Checkbox('120', size=(12, 1)), sg.Checkbox('180', size=(12, 1)), sg.Checkbox('240', size=(12, 1))],
              [sg.Checkbox('300', size=(12, 1)), sg.Checkbox('600', size=(12, 1))],
              [sg.Text('_'  * 100, size=(65, 1))],
              [sg.Text('Filters', font=('Helvetica', 15), justification='left')],
              [sg.Checkbox('Luminance', size=(12, 1)), sg.Checkbox('Red', size=(12, 1)), sg.Checkbox('Green', size=(12, 1)),
               sg.Checkbox('Blue', size=(12, 1)), sg.Checkbox('Grating', size=(12, 1))],
              [sg.Checkbox('B', size=(12, 1)), sg.Checkbox('V', size=(12, 1)), sg.Checkbox('R', size=(12, 1)),
               sg.Checkbox('I', size=(12, 1)), sg.Checkbox('Y', size=(12, 1))],
              [sg.Checkbox("g'", size=(12, 1)), sg.Checkbox("r'", size=(12, 1)), sg.Checkbox("i'", size=(12, 1)),
               sg.Checkbox("z'", size=(12, 1)), sg.Checkbox('OIII', size=(12, 1))],
              [sg.Checkbox('Ha', size=(12, 1)), sg.Checkbox('SII', size=(12, 1))],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('GBO finder gui', font=("Helvetica", 12)).Layout(layout)

    event, values = window.Read()

    return values
