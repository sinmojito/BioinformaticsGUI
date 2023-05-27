import PySimpleGUI as sg
import importlib

sg.theme('DarkTeal9')  # Change the theme here

# Define the functions and layouts for each script
scripts = {
    'Adapter Trimmer': {'function': 'adapter_trimmer.main', 'layout_module': 'adapter_trimmer'},
    'Quality Filter': {'function': 'quality_filter.main', 'layout_module': 'quality_filter'},
    'Quality Trimmer': {'function': 'quality_trimmer.main', 'layout_module': 'quality_trimmer'},
    'DEA - edgeR via Rpy2 implementation': {'function': 'dea_analysis.main', 'layout_module': 'dea_analysis'},
    'DEA - PyDESeq2 Implementation': {'function': 'pydeseq2_gui.main', 'layout_module': 'pydeseq2_gui'}
}



# Define a function to create a script window
def create_script_window(name):
    # Get the function and module corresponding to the button clicked
    function_name = scripts[name]['function']
    layout_module_name = scripts[name]['layout_module']
    function_module = importlib.import_module(function_name.rsplit('.', 1)[0])
    function = getattr(function_module, function_name.rsplit('.', 1)[1])

    # Get the layout of the selected script
    layout_module = importlib.import_module(layout_module_name)
    layout = layout_module.create_layout()  # Call the create_layout function

    # Create the PySimpleGUI window for the script
    window = sg.Window(name, layout)
    window.finalize()

    # Event loop for the script window
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event is not None:
            # Call the function corresponding to the button clicked
            function()

    # Close the script window
    window.close()

# Define the PySimpleGUI layout for the main window
main_layout = [
    [sg.Button(name, size=(50, 2)), sg.Text({
                                                  'Adapter Trimmer': 'Opens the Adapter Trimmer application, which helps in removing adapter sequences from high-throughput sequencing data.',
                                                  'Quality Filter': 'Opens the Quality Filter application, which helps in filtering out low quality reads from your sequencing data to improve downstream analysis.',
                                                  'Quality Trimmer': 'Opens the Quality Trimmer application, which trims low quality bases from the ends of sequences. It helps in maintaining the high quality of the sequencing data.',
                                                  'DEA - edgeR via Rpy2 implementation': 'Opens the DEA - edgeR via Rpy2 implementation application. This is used to identify genes that are differentially expressed between different experimental conditions.',
                                                  'DEA - PyDESeq2 Implementation': 'Opens the DEA - PyDESeq2 Implementation application. This is another method used to identify differentially expressed genes.',
                                              }[name], size=(50, 2))] for name in scripts.keys()
] + [
    [sg.Button('Exit')],
    [sg.Text('Consolidated GUI', key='Application name', size=(None, 1), justification='left', font=("Alike", 11, "bold"))],  # Application name
    [sg.Text('Thesis Project. Created by Mohit Panwar. Supervised by Julia Åkesson.', key='credits', size=(None, 1), justification='left', font=("Alike", 9))]  # Credits
]



# Create the PySimpleGUI window for the main window
window = sg.Window('Script Selector', main_layout)

# Event loop for the main window
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event is not None:
        # Open the script window for the selected script
        create_script_window(event)

# Close the main window
window.close()
