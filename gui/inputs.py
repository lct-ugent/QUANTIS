'''
Created on 12 Oct 2016

@author: shsymoen
'''


import configparser
import os
import sys
from tkinter import (
    BOTTOM,
    Button,
    Checkbutton,
    Entry,
    Frame,
    IntVar,
    LEFT,
    RIGHT,
    TOP,
    Tk,
)
sys.path.append('.')
sys.path.append('..')

from runner.experiment_visualiser import run_visualise


def create_config(config_file_name):
    """
    Create a config file
    
    Parameters
    ----------
    config_file_name : string
        The .ini file
    
    Returns
    -------
    None
    """
    import configparser
    import os
    
    config = configparser.ConfigParser()
    config.add_section("Files")
    config.add_section("Working_Directory")
    config.set(
        "Files", 
        "conditions", 
        "conditions_propane.csv gen_configs_propane.csv"
    )
    config.set("Working_Directory", "cd", os.getcwd())
 
    with open(config_file_name, "w") as config_file:
        config.write(config_file)
      
  
def get_config(config_file_name):
    """
    Returns the config object
    
    Parameters
    ----------
    config_file_name : string 
        The .ini file
    
    Returns
    -------
    a configparser object
    """
    import configparser
    import os
    
    if not os.path.exists(config_file_name):
        create_config(config_file_name)
 
    config = configparser.ConfigParser()
    config.read(config_file_name)
    return config


def get_setting(config_file_name, section, setting):
    """
    Print out a setting of the configuration file
    
    Parameters
    ----------
    config_file_name : string
        The .ini file
    section : string
        The title of the section of interest
    setting : string
        The setting you would like to return
    
    Returns
    -------
    value : string
        The setting from the section
    """
    config = get_config(config_file_name)
    value = config.get(section, setting)
    return value


def update_setting(dir, config_file_name, section, setting, value):
    """
    Update a setting and changes directory 
    
    Parameters
    ----------
    dir : string
        Directory of the .ini file
    config_file_name : string
        The .ini file
    section : string
        The title of the section of interest
    setting : string
        The setting you would like to update
    value : string
        New updated value of the setting
    
    Returns
    -------
    None
    """
    import os
    os.chdir(dir)
    config = get_config(config_file_name)
    config.set(section, setting, value)
    with open(config_file_name, "w") as config_file:
        config.write(config_file)
    
    return None
    
    
def inputWindow():
    """ Graphical user interface to insert the program settings """
    
    current_directory = os.getcwd()
    
    # Create initial window
    window = Tk()

    # Create checklabels
    checkLabels = [IntVar() for _ in range(6)]
    texts = [
        "Process results",
        "Visualise experiment",
        "Simulate CHEMKIN",
        "Simulate COILSIM1D",
        "Compare experiment with simulation",
        "Simulate COILSIM1D v3.9"
    ]

    def run_BSSC_visualizer():
        import os

        results_file_name = fileName_entry.get()
        
        # Change working directory
        os.chdir(location_entry.get())  
        
        #Update the ini file and change directory
        update_setting(current_directory, 
                       "settings.ini", 
                       "Files", 
                       "conditions", 
                       results_file_name
        )
        
        #Update the ini file and change directory     
        update_setting(current_directory, 
                       "settings.ini", 
                       "Working_Directory", 
                       "cd", 
                       location_entry.get()
        )
        
        # Change working directory again
        os.chdir(location_entry.get())  

        chemkin_names = chemkinNames_entry.get().split(' ')

        process_exp_results = checkLabels[0].get()
        write_chemkin = checkLabels[2].get()
        write_coilsim = checkLabels[3].get()
        visualise = checkLabels[1].get()
        compare = checkLabels[4].get()
        write_coilsim_v39 = checkLabels[5].get()
        mol = wt_mol_checklabel.get()
        fid_sum = abs(fid_sum_checklabel.get() - 1)
        window.destroy()
        
        # Run actual framework
        run_visualise(
            process_exp_results,
            write_chemkin,
            write_coilsim,
            visualise,
            compare,
            results_file_name,
            chemkin_names,
            mol,
            fid_sum,
            write_coilsim_v39
        )

    def close_window():
        window.destroy()

    button_calls = [run_BSSC_visualizer, close_window]
    button_texts = ['Run', 'Cancel']
    sides = [LEFT, RIGHT]

    topFrame = Frame(window)
    topFrame.pack()
    bottomFrame = Frame(window)
    bottomFrame.pack(side=BOTTOM)

    fileName_entry = Entry(topFrame, width=50)
    fileName_entry.insert(
        0,
        get_setting(
            "settings.ini", 
            "Files", 
            "conditions"
        )
    )
    fileName_entry.pack(side=TOP)
    location_entry = Entry(topFrame, width=100)
    location_entry.insert(
        0, 
        get_setting(
            "settings.ini", 
            "Working_Directory", 
            "cd"
        )
    )
    location_entry.pack(side=TOP)

    chemkinNames_entry = Entry(topFrame, width=100)
    chemkinNames_entry.insert(0, "aramcomech2")
    chemkinNames_entry.pack(side=BOTTOM)

    # Initialise all checkbuttons
    for checkLabel, text in zip(checkLabels, texts):
        c = Checkbutton(
            topFrame,
            text=text,
            variable=checkLabel,
            onvalue=1,
            offvalue=0,
            padx=10,
            pady=10
        )
        c.pack(side=LEFT)

    # Initialise all buttons
    for button_call, text, side in zip(button_calls, button_texts, sides):
        b = Button(
            bottomFrame,
            text=text,
            command=button_call,
            padx=10,
            pady=10
        )
        b.pack(side=side)

    wt_mol_checklabel = IntVar()
    c = Checkbutton(
        bottomFrame,
        text='Mol.% wet (default wt.%)',
        variable=wt_mol_checklabel,
        onvalue=1,
        offvalue=0,
        padx=10,
        pady=10
    )
    c.pack(side=BOTTOM)

    fid_sum_checklabel = IntVar()
    d = Checkbutton(
        bottomFrame,
        text='TCD > FID (RGA)',
        variable=fid_sum_checklabel,
        onvalue=1,
        offvalue=0,
        padx=10,
        pady=10
    )
    d.pack(side=BOTTOM)

    window.mainloop()


if __name__ == '__main__':
    """ main program """
    
    inputWindow()
