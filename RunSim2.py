import glob
import rocketCEA as cea
import tkinter as tk
from tkinter import ttk
import numpy as np
import math
import AnimRocket
global sim_vars, exception_vars

# Define simulation variables in a dictionary
sim_vars = {
    "Ox_tank_length": 1.7526, "Ox_tank_diameter": 0.10072832596729638, "Ox_tank_vol": 0.013966124892624783,
    "Starting_Tank_Pressure": 5.516e6, "Starting_Ox_Mass": 18,
    "Starting_Chamber_Pressure": 101325, "CC_vol": 0.019661,
    "Grain_ID": 0.1, "Grain_OD": 0.125, "Grain_Length": 1.5,
    "blowing_number": 15, "a": 0.000155, "n": 0.45, "m": 0,
    "Injector_Hole_Diameter": 0.0015, "Number_of_Injector_Holes": 60,
    "Injector_Discharge_Coefficient": 0.55, "Nozzle_Throat_Diameter": 0.0954,
    "Nozzle_Expansion_Ratio": 1.2, "Nozzle_Efficiency": 0.95,
    "Nozzle_Discharge_Ratio": 0.9, "c_eff": 0.9,
    "dry_mass": 40, "viscosity": 3.70e-5, "For_flight": 0,
    "Aluminum_weight_percent": 0, "Carbon_black_weight_percent": 10
}

exception_vars = ["Aluminum_weight_percent", "Carbon_black_weight_percent", "Ox_tank_vol"]

# Function to update values when the user modifies an entry
def update_value(var_name, entry, toggle):
    global sim_vars, selected_option_fuel
    if (toggle ==0):
        try:
            sim_vars[var_name] = float(entry.get())
            print(f"{var_name} updated: {sim_vars[var_name]}")
        except ValueError:
            print("pass")
    elif (toggle == 1):
        try:
            if selected_option_fuel.get() == "Aluminum":
                sim_vars["Aluminum_weight_percent"] = float(entry.get())
                print(f"Aluminum Weight Percent updated: {sim_vars['Aluminum_weight_percent']}")
            elif selected_option_fuel.get() == "Carbon Black":
                sim_vars["Carbon_black_weight_percent"] = float(entry.get())
                print(f"Carbon Black Weight Percent updated: {sim_vars['Carbon_black_weight_percent']}")
        except ValueError:
            print("pass")


def clear(entry, key):
    global sim_vars
    global exception_vars
    global selected_option_fuel
    
    # Handle special case for fuel
    if key == "fuel":
        entry.delete(0, tk.END)
        if selected_option_fuel.get() == "Aluminum":
            sim_vars["Aluminum_weight_percent"] = 0
            entry.insert(0, str(sim_vars["Aluminum_weight_percent"]))
            print(f"Cleared Aluminum_weight_percent")
        elif selected_option_fuel.get() == "Carbon Black":
            sim_vars["Carbon_black_weight_percent"] = 0
            entry.insert(0, str(sim_vars["Carbon_black_weight_percent"]))
            print(f"Cleared Carbon_black_weight_percent")
    else:
        entry.delete(0, tk.END)
        sim_vars[key] = 0  # Update the global variable
        entry.insert(0, str(sim_vars[key]))
        print(f"Cleared {str(key)}")
    
    return

def on_fuel_selected(*args):
    global sim_vars, fuel_entry, selected_option_fuel
    print(f"Fuel selected: {selected_option_fuel.get()}")
    if selected_option_fuel.get() == "Aluminum":
        fuel_entry.delete(0, tk.END)
        try:
            fuel_entry.insert(0, str(sim_vars["Aluminum_weight_percent"]))
        except:
            pass
    elif selected_option_fuel.get() == "Carbon Black":
        fuel_entry.delete(0, tk.END)
        try:
            fuel_entry.insert(0, str(sim_vars["Carbon_black_weight_percent"]))
        except:
            pass
    print(f"Aluminum Weight Percent updated: {sim_vars['Aluminum_weight_percent']}")
    print(f"Carbon Black Weight Percent updated: {sim_vars['Carbon_black_weight_percent']}")

# Function to create labeled entries
def create_entry(parent, row, label_text, var_name, toggle):
    global selected_option_fuel, fuel_entry
    if (toggle == 0):
        label = ttk.Label(parent, text=label_text)
        label.grid(row=row, column=0, sticky="w")
    
        var = tk.StringVar(value=str(sim_vars[var_name]))
        entry = ttk.Entry(parent, width=20, textvariable=var)
        entry.grid(row=row, column=1, padx=10)

        # Bind event instead of using trace_add
        entry.bind("<KeyRelease>", lambda event: update_value(var_name, var, 0))
    
        clear_btn = ttk.Button(parent, text="Clear", command=lambda: clear(entry, var_name))
        clear_btn.grid(row=row, column=2, padx=10)
        return entry
    elif (toggle == 1):
        # Fuel Weight Percent Characterization
        weight_percent = tk.StringVar(value = str(sim_vars[var_name]))
        fuel_entry = ttk.Entry(parent, width=20, textvariable=weight_percent)
        fuel_entry.grid(row=row, column=1, padx=10)
        
        fuel_entry.bind("<KeyRelease>", lambda event: update_value(var_name, weight_percent, 1))
        
        clear_btn = ttk.Button(frame, text="Clear", command=lambda: clear(fuel_entry, "fuel"))
        clear_btn.grid(row=row, column=2, padx=10)
        # Creates options for the user to select
        options = ["Aluminum", "Carbon Black"]
        # Creates a variable to hold the selected option
        selected_option_fuel = tk.StringVar()
        # Default option
        selected_option_fuel.set(options[0])
        # Adds a trace to the variable to call on_fuel_selected whenever the value changes
        selected_option_fuel.trace_add("write", on_fuel_selected)
        # Creates a dropdown menu
        fuel_dropdown = ttk.OptionMenu(frame, selected_option_fuel, options[0], *options)
        fuel_dropdown.grid(row=row, column=0, sticky='w')
        
        return fuel_entry
        
    
def sim():
    global sim_vars
    sim_vars["Ox_tank_vol"] = sim_vars["Ox_tank_length"] * math.pi * (sim_vars["Ox_tank_diameter"]/2)**2
    sim_vars_with_var = {f"{key}_var": value for key, value in sim_vars.items()}
    cea.set_global_variables(**sim_vars_with_var)
    new_syst = cea.on_button_click()
    #AnimRocket.anim(new_syst)


    
# Titles
titles = [
    ("Hybrid Rocket Simulation","Impact", 18, 0),
    ("Ox Tank Dimensions & Properties","Comic Sans MS", 14, 1),
    ("Combustion Chamber Properties","Comic Sans MS", 14, 6),
    ("Fuel Weight Percent Characterization", "Comic Sans MS", 14, 9),
    ("Grain Dimensions and Characterization","Comic Sans MS", 14, 11),
    ("Injector Dimensions","Comic Sans MS", 14, 19),
    ("Nozzle Dimensions and Characterization","Comic Sans MS", 14, 23)
]
    
# Create the main window
root = tk.Tk()

root.title("Advanced UI Example")

# Create a Style
style = ttk.Style()
style.configure("Custom.TFrame", background="lightblue")  # Set the background color

#Create canvas
canvas = tk.Canvas(root, height=720, width=1280, bg = "lightblue")
canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

#Create a scrollbar
scrollbar = ttk.Scrollbar(root, orient=tk.VERTICAL, command=canvas.yview)
scrollbar.grid(row=0, column=1, sticky="ns")

#Configure the canvas
canvas.configure(yscrollcommand=scrollbar.set)
canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1*(event.delta/120)), "units"))


# Create a frame
frame = ttk.Frame(canvas, style= "Custom.TFrame")
frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
x0 = frame.winfo_screenwidth()/2
y0 = frame.winfo_screenheight()/2
canvas.create_window((x0,y0), window=frame, anchor = "center")
#canvas.create_window((0,0), window=frame, anchor="nw")

titles_vars = ["Ox Tank Length [m]", "Ox Tank Diameter [m]", "Starting Tank Pressure [Pa]", "Starting Ox Mass [kg]", 
          "Starting Chamber Pressure [Pa]", "CC Volume [m^3]", "Grain ID [m]", "Grain OD [m]", "Grain Length [m]",
          "Blowing Number", "a", "n", "m", "Injector Hole Diameter [m]", "Number of Injector Holes",
          "Injector Discharge Coefficient", "Nozzle Throat Diameter [m]", "Nozzle Expansion Ratio", "Nozzle Efficiency",
          "Nozzle Discharge Ratio", "c_eff", "Dry Mass [kg]", "Viscosity [Pa*s]", "For Flight"]

# Create entries dynamically
entries = {}
row_counter = 0
title_counter = 0
for key in sim_vars.keys():
    if key in exception_vars:
        continue
    else:
        while row_counter in {0, 1, 6, 9, 10, 11, 19, 23}:
            row_counter += 1
        entries[key] = create_entry(frame, row_counter, titles_vars[title_counter], key, 0)
        title_counter += 1
        row_counter += 1



for title, font, size, row in titles:
    ttk.Label(frame, text=title, font=(font, size)).grid(row=row, column=0, columnspan=3, pady=10)
    

fuel_entry = create_entry(frame, 10, " ", "Aluminum_weight_percent", 1)

# Sim Button
ttk.Button(frame, text="Run Simulation", command=sim).grid(row=row_counter+1, column=1, pady=20)

# Configure the canvas
frame.update_idletasks()
canvas.config(scrollregion=canvas.bbox("all"))

# Configure the grid
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
frame.columnconfigure(0, weight=1)
frame.columnconfigure(1, weight=1)

# Run the application
root.mainloop()