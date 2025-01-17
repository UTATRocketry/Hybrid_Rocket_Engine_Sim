from multiprocessing import Value
import rocketCEA as cea
import tkinter as tk
from tkinter import ttk
import numpy as np
import math

global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent

#Initializes the variables
Ox_tank_vol = 0.01396612489262478177383457064491
Ox_tank_length = 1.7526
Ox_tank_diameter = ((4*Ox_tank_vol)/(math.pi*Ox_tank_length))**0.5
Aluminum_weight_percent = 0
Carbon_black_weight_percent = 10

def on_value_changed_Ox_tank(*args):
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter
    try:
        Ox_tank_vol = float(tank_vol_entry.get())
    except ValueError:
        pass
    try:
        Ox_tank_length = float(tank_len_entry.get())
    except ValueError:
        pass
    try:
        Ox_tank_diameter = float(tank_dia_entry.get())
    except ValueError:
        pass
        
    
    print(f"Tank Volume updated: {Ox_tank_vol}")
    print(f"Tank Length updated: {Ox_tank_length}")
    print(f"Tank Diameter updated: {Ox_tank_diameter}")

def compute_the_third():
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter
    if (Ox_tank_vol == 0):
        Ox_tank_vol = (math.pi*Ox_tank_length*(Ox_tank_diameter**2))/4
        print(f"Ox_tank_vol: {Ox_tank_vol}")
        tank_vol_entry.delete(0, tk.END)
        tank_vol_entry.insert(0, str(Ox_tank_vol))
    elif (Ox_tank_length == 0):
        Ox_tank_length = (4*Ox_tank_vol)/(math.pi*(Ox_tank_diameter**2))
        print(f"Ox_tank_length: {Ox_tank_length}")
        tank_len_entry.delete(0, tk.END)
        tank_len_entry.insert(0, str(Ox_tank_length))
    elif (Ox_tank_diameter == 0):
        Ox_tank_diameter = (4*Ox_tank_vol/(math.pi*Ox_tank_length))**0.5
        print(f"Ox_tank_diameter: {Ox_tank_diameter}")
        tank_dia_entry.delete(0, tk.END)
        tank_dia_entry.insert(0, str(Ox_tank_diameter))
    else:
        print("Error: All values are given")
    return

def on_fuel_selected(*args):
    global Aluminum_weight_percent, Carbon_black_weight_percent, fuel_entry
    print(f"Fuel selected: {selected_option_fuel.get()}")
    if selected_option_fuel.get() == "Aluminum":
        fuel_entry.delete(0, tk.END)
        try:
            fuel_entry.insert(0, str(Aluminum_weight_percent))
        except:
            pass
    elif selected_option_fuel.get() == "Carbon Black":
        fuel_entry.delete(0, tk.END)
        try:
            fuel_entry.insert(0, str(Carbon_black_weight_percent))
        except:
            pass
    print(f"Aluminum Weight Percent updated: {Aluminum_weight_percent}")
    print(f"Carbon Black Weight Percent updated: {Carbon_black_weight_percent}")

def on_value_changed_fuel(*args):
    global Aluminum_weight_percent, Carbon_black_weight_percent, fuel_entry, selected_option_fuel
    if selected_option_fuel.get() == "Aluminum": 
        try:
            Aluminum_weight_percent = float(fuel_entry.get())
        except:
            pass
    elif selected_option_fuel.get() == "Carbon Black":
        try:
            Carbon_black_weight_percent = float(fuel_entry.get())
        except:
            pass
    print(f"Aluminum Weight Percent updated: {Aluminum_weight_percent}")
    print(f"Carbon Black Weight Percent updated: {Carbon_black_weight_percent}")

def clear(var):
    global tank_vol_entry, tank_len_entry, tank_dia_entry, fuel_entry, selected_option_fuel
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent
    if var == "vol":
        tank_vol_entry.delete(0, tk.END)
        Ox_tank_vol = 0
        tank_vol_entry.insert(0, str(Ox_tank_vol))
        print("Cleared Tank Vol")
    elif var == "len":
        tank_len_entry.delete(0, tk.END)
        Ox_tank_length = 0
        tank_len_entry.insert(0, str(Ox_tank_length))
        print("Cleared Tank Len")
    elif var == "dia":
        tank_dia_entry.delete(0, tk.END)
        Ox_tank_diameter = 0
        tank_dia_entry.insert(0, str(Ox_tank_diameter))
        print("Cleared Tank Dia")
    elif var == "fuel":
        fuel_entry.delete(0, tk.END)
        if selected_option_fuel.get() == "Aluminum":
            Aluminum_weight_percent = 0
            fuel_entry.insert(0, str(Aluminum_weight_percent))
        elif selected_option_fuel.get() == "Carbon Black":
            Carbon_black_weight_percent = 0
            fuel_entry.insert(0, str(Carbon_black_weight_percent))
        
        print("Cleared Fuel")
    return

def sim():
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent
    #Ox_tank_vol = 0.01396612489262478177383457064491
    #Ox_tank_length = 1.7526
    #Ox_tank_diameter = ((4*Ox_tank_vol)/(math.pi*Ox_tank_length))**0.5
    #Aluminum_weight_percent = 0
    #Carbon_black_weight_percent = 10

    dry_mass = 40
    viscosity = 3.70e-5
    blowing_number = 15
    a = 0.000155
    n = 0.45
    m = 0

    CC_vol = 0.019661

    Nozzle_Throat_Diameter = 0.0954278
    Nozzle_Expansion_Ratio = 1.2
    Nozzle_Efficiency = 0.95
    Nozzle_Discharge_Ratio = 0.9

    # Assuming showerhead injector for now
    Injector_Hole_Diamter = 0.0015
    Number_of_Injector_Holes = 60
    Injector_Discharge_Coefficient = 0.55


    c_eff = 0.9
    Grain_ID = 0.1
    Grain_OD = 0.125
    Grain_Length = 1.5

    Starting_Tank_Pressure = 5.516e6
    Starting_Chamber_Pressure = 101325
    Starting_Ox_Mass = 18
    For_flight = 0
    cea.set_global_variables(Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent, CC_vol, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, Injector_Hole_Diamter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, c_eff, Grain_ID, Grain_OD, Grain_Length, Starting_Tank_Pressure, Starting_Chamber_Pressure, Starting_Ox_Mass, For_flight, dry_mass, viscosity, blowing_number, a, n, m)
    cea.on_button_click()

# Create the main window
root = tk.Tk()
root.title("Advanced UI Example")

# Create a frame
frame = ttk.Frame(root, padding="10")
frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create a label
title1 = ttk.Label(frame, text="This is a Sim Trust", font=("Impact", 18))
title1.grid(row=0, column=0, columnspan=2, pady=10)

title2 = ttk.Label(frame, text="Tank Dimensions", font = ("Comic Sans MS", 14))
title2.grid(row=1, column=0, columnspan=2, pady=10)

# Create a tank volume label
tank_vol_label = ttk.Label(frame, text="Tank Volume")
tank_vol_label.grid(row=2, column=0, pady=10)
# Create a StringVar to track the value of the entry widget
tank_vol_var = tk.StringVar()
# Set up a trace on the StringVar to call on_value_change whenever the value changes
tank_vol_var.trace_add("write", on_value_changed_Ox_tank)
# Create an entry widget
tank_vol_entry = ttk.Entry(frame, width=20, textvariable = tank_vol_var)
tank_vol_entry.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
# Create a clear button
clear_button1 = ttk.Button(frame, text="Clear", command=lambda: clear("vol"))
clear_button1.grid(row=2, column=1, padx=5, pady=5)

# For Tank length
tank_len_label = ttk.Label(frame, text="Tank Length")
tank_len_label.grid(row=3, column=0, pady=10)
tank_len_var = tk.StringVar()
tank_len_var.trace_add("write", on_value_changed_Ox_tank)
tank_len_entry = ttk.Entry(frame, width=20, textvariable = tank_len_var)
tank_len_entry.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
clear_button2 = ttk.Button(frame, text="Clear", command=lambda: clear("len"))
clear_button2.grid(row=3, column=1, padx=5, pady=5)

# For Tank Diameter
tank_dia_label = ttk.Label(frame, text="Tank Diameter")
tank_dia_label.grid(row=4, column=0, pady=10)
tank_dia_var = tk.StringVar()
tank_dia_var.trace_add("write", on_value_changed_Ox_tank)
tank_dia_entry = ttk.Entry(frame, width=20, textvariable = tank_dia_var)
tank_dia_entry.grid(row=4, column=0, columnspan=2, padx=5, pady=5)
clear_button3 = ttk.Button(frame, text="Clear", command=lambda: clear("dia"))
clear_button3.grid(row=4, column=1, padx=5, pady=5)

#Tank Calculation Button
compute_button = ttk.Button(frame, text="Find 3rd Var", command=compute_the_third)
compute_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

title3 = ttk.Label(frame, text="Fuel Characterization", font = ("Comic Sans MS", 14))
title3.grid(row=6, column=0, columnspan=2, pady=10)

# Fuel Characterization
weight_percent = tk.StringVar()
weight_percent.trace_add("write", on_value_changed_fuel)
fuel_entry = ttk.Entry(frame, width=20, textvariable=weight_percent)
fuel_entry.grid(row=7, column=0, columnspan=2, padx=5, pady=5)
clear_button4 = ttk.Button(frame, text="Clear", command=lambda: clear("fuel"))
clear_button4.grid(row=7, column=1, padx=5, pady=5)

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
fuel_dropdown.grid(row=7, column=0, pady=10)



# Sim Button
compute_button = ttk.Button(frame, text="Click Me", command=sim)
compute_button.grid(row=10, column=0, columnspan=2, padx=5, pady=5)

# Configure the grid
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
frame.columnconfigure(0, weight=1)
frame.columnconfigure(1, weight=1)

#Insert Default Values
tank_vol_entry.insert(0, str(Ox_tank_vol))
tank_len_entry.insert(0, str(Ox_tank_length))
tank_dia_entry.insert(0, str(Ox_tank_diameter))

# Run the application
root.mainloop()
