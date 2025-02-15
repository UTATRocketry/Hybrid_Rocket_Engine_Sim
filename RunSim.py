import rocketCEA as cea
import tkinter as tk
from tkinter import ttk
import numpy as np
import math
import AnimRocket
global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent, Starting_Tank_Pressure, Starting_Ox_Mass, Starting_Chamber_Pressure, CC_vol, Grain_ID, Grain_OD, Grain_Length, blowing_number, a, n, m, Injector_Hole_Diameter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, c_eff, dry_mass, viscosity, For_flight

#Initializes the variables
Ox_tank_vol = 0.01396612489262478177383457064491
Ox_tank_length = 1.7526
Ox_tank_diameter = ((4*Ox_tank_vol)/(math.pi*Ox_tank_length))**0.5
Aluminum_weight_percent = 0
Carbon_black_weight_percent = 10
Starting_Tank_Pressure = 5.516e6
Starting_Ox_Mass = 18
Starting_Chamber_Pressure = 101325
CC_vol = 0.019661
Grain_ID = 0.1
Grain_OD = 0.125
Grain_Length = 1.5
blowing_number = 15
a = 0.000155
n = 0.45
m = 0
Injector_Hole_Diameter = 0.0015
Number_of_Injector_Holes = 60
Injector_Discharge_Coefficient = 0.55
Nozzle_Throat_Diameter = 0.0954278
Nozzle_Expansion_Ratio = 1.2
Nozzle_Efficiency = 0.95
Nozzle_Discharge_Ratio = 0.9
c_eff = 0.9
dry_mass = 40
viscosity = 3.70e-5
For_flight = 0

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
    
def on_value_changed_tank_pressure(*args):
    global Starting_Tank_Pressure
    try:
        Starting_Tank_Pressure = float(tank_pressure_entry.get())
    except ValueError:
        pass
    print(f"Tank Pressure updated: {Starting_Tank_Pressure}")

def on_value_changed_starting_ox_mass(*args):
    global Starting_Ox_Mass
    try:
        Starting_Ox_Mass = float(starting_ox_mass_entry.get())
    except ValueError:
        pass
    print(f"Starting Ox Mass updated: {Starting_Ox_Mass}")

def on_value_changed_starting_chamber_pressure(*args):
    global Starting_Chamber_Pressure
    try:
        Starting_Chamber_Pressure = float(starting_chamber_pressure_entry.get())
    except ValueError:
        pass
    print(f"Starting Chamber Pressure updated: {Starting_Chamber_Pressure}")
    
def on_value_changed_cc_vol(*args):
    global CC_vol
    try:
        CC_vol = float(cc_vol_entry.get())
    except ValueError:
        pass
    print(f"CC Volume updated: {CC_vol}")

def on_value_changed_grain_id(*args):
    global Grain_ID
    try:
        Grain_ID = float(grain_id_entry.get())
    except ValueError:
        pass
    print(f"Grain ID updated: {Grain_ID}")

def on_value_changed_grain_od(*args):
    global Grain_OD
    try:
        Grain_OD = float(grain_od_entry.get())
    except ValueError:
        pass
    print(f"Grain OD updated: {Grain_OD}")

def on_value_changed_grain_length(*args):
    global Grain_Length
    try:
        Grain_Length = float(grain_length_entry.get())
    except ValueError:
        pass
    print(f"Grain Length updated: {Grain_Length}")

def on_value_changed_blowing_number(*args):
    global blowing_number
    try:
        blowing_number = float(blowing_number_entry.get())
    except ValueError:
        pass
    print(f"Blowing Number updated: {blowing_number}")
    
def on_value_changed_a(*args):
    global a
    try:
        a = float(a_entry.get())
    except ValueError:
        pass
    print(f"a updated: {a}")
    
def on_value_changed_n(*args):
    global n
    try:
        n = float(n_entry.get())
    except ValueError:
        pass
    print(f"n updated: {n}")

def on_value_changed_m(*args):
    global m
    try:
        m = float(m_entry.get())
    except ValueError:
        pass
    print(f"m updated: {m}")    

def on_value_changed_injector_hole_diameter(*args):
    global Injector_Hole_Diameter
    try:
        Injector_Hole_Diameter = float(injector_hole_diameter_entry.get())
    except ValueError:
        pass
    print(f"Injector Hole Diameter updated: {Injector_Hole_Diameter}")

def on_value_changed_number_of_injector_holes(*args):
    global Number_of_Injector_Holes
    try:
        Number_of_Injector_Holes = float(number_of_injector_holes_entry.get())
    except ValueError:
        pass
    print(f"Number of Injector Holes updated: {Number_of_Injector_Holes}")
    
def on_value_changed_injector_discharge_coefficient(*args):
    global Injector_Discharge_Coefficient
    try:
        Injector_Discharge_Coefficient = float(injector_discharge_coefficient_entry.get())
    except ValueError:
        pass
    print(f"Injector Discharge Coefficient updated: {Injector_Discharge_Coefficient}")
    
def on_value_changed_nozzle_throat_diameter(*args):
    global Nozzle_Throat_Diameter
    try:
        Nozzle_Throat_Diameter = float(nozzle_throat_diameter_entry.get())
    except ValueError:
        pass
    print(f"Nozzle Throat Diameter updated: {Nozzle_Throat_Diameter}")
    
def on_value_changed_nozzle_expansion_ratio(*args):
    global Nozzle_Expansion_Ratio
    try:
        Nozzle_Expansion_Ratio = float(nozzle_expansion_ratio_entry.get())
    except ValueError:
        pass
    print(f"Nozzle Expansion Ratio updated: {Nozzle_Expansion_Ratio}")
    
def on_value_changed_nozzle_efficiency(*args):
    global Nozzle_Efficiency
    try:
        Nozzle_Efficiency = float(nozzle_efficiency_entry.get())
    except ValueError:
        pass
    print(f"Nozzle Efficiency updated: {Nozzle_Efficiency}")
    
def on_value_changed_nozzle_discharge_ratio(*args):
    global Nozzle_Discharge_Ratio
    try:
        Nozzle_Discharge_Ratio = float(nozzle_discharge_ratio_entry.get())
    except ValueError:
        pass
    print(f"Nozzle Discharge Ratio updated: {Nozzle_Discharge_Ratio}")
    
def on_value_changed_c_eff(*args):
    global c_eff
    try:
        c_eff = float(c_eff_entry.get())
    except ValueError:
        pass
    print(f"c_eff updated: {c_eff}")
    
def on_value_changed_dry_mass(*args):
    global dry_mass
    try:
        dry_mass = float(dry_mass_entry.get())
    except ValueError:
        pass
    print(f"Dry Mass updated: {dry_mass}")
    
def on_value_changed_viscosity(*args):
    global viscosity
    try:
        viscosity = float(viscosity_entry.get())
    except ValueError:
        pass
    print(f"Viscosity updated: {viscosity}")
    
def on_value_changed_for_flight(*args):
    global For_flight
    try:
        For_flight = float(for_flight_entry.get())
    except ValueError:
        pass
    print(f"For Flight updated: {For_flight}")

def clear(var):
    global tank_vol_entry, tank_len_entry, tank_dia_entry, fuel_entry, selected_option_fuel, tank_pressure_entry, starting_ox_mass_entry, starting_chamber_pressure_entry, cc_vol_entry, grain_id_entry, grain_od_entry, grain_length_entry, blowing_number_entry, a_entry, n_entry, m_entry, injector_hole_diameter_entry, number_of_injector_holes_entry, injector_discharge_coefficient_entry, nozzle_throat_diameter_entry, nozzle_expansion_ratio_entry, nozzle_efficiency_entry, nozzle_discharge_ratio_entry, c_eff_entry, dry_mass_entry, viscosity_entry, for_flight_entry
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent, Starting_Tank_Pressure, Starting_Ox_Mass, Starting_Chamber_Pressure, CC_vol, Grain_ID, Grain_OD, Grain_Length, blowing_number, a, n, m, Injector_Hole_Diameter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, c_eff, dry_mass, viscosity, For_flight

    # Mapping of variables to corresponding entries, values, and print messages
    clear_map = {
        "vol": (tank_vol_entry, "Ox_tank_vol", "Cleared Tank Vol"),
        "len": (tank_len_entry, "Ox_tank_length", "Cleared Tank Len"),
        "dia": (tank_dia_entry, "Ox_tank_diameter", "Cleared Tank Dia"),
        "fuel": (fuel_entry, None, "Cleared Fuel"),
        "tank_pressure": (tank_pressure_entry, "Starting_Tank_Pressure", "Cleared Tank Pressure"),
        "starting_ox_mass": (starting_ox_mass_entry, "Starting_Ox_Mass", "Cleared Starting Ox Mass"),
        "starting_chamber_pressure": (starting_chamber_pressure_entry, "Starting_Chamber_Pressure", "Cleared Starting Chamber Pressure"),
        "cc_volume": (cc_vol_entry, "CC_vol", "Cleared CC Volume"),
        "grain_id": (grain_id_entry, "Grain_ID", "Cleared Grain ID"),
        "grain_od": (grain_od_entry, "Grain_OD", "Cleared Grain OD"),
        "grain_length": (grain_length_entry, "Grain_Length", "Cleared Grain Length"),
        "blowing_number": (blowing_number_entry, "blowing_number", "Cleared Blowing Number"),
        "a": (a_entry, "a", "Cleared a"),
        "n": (n_entry, "n", "Cleared n"),
        "m": (m_entry, "m", "Cleared m"),
        "injector_hole_diameter": (injector_hole_diameter_entry, "Injector_Hole_Diameter", "Cleared Injector Hole Diameter"),
        "number_of_injector_holes": (number_of_injector_holes_entry, "Number_of_Injector_Holes", "Cleared Number of Injector Holes"),
        "injector_discharge_coefficient": (injector_discharge_coefficient_entry, "Injector_Discharge_Coefficient", "Cleared Injector Discharge Coefficient"),
        "nozzle_throat_diameter": (nozzle_throat_diameter_entry, "Nozzle_Throat_Diameter", "Cleared Nozzle Throat Diameter"),
        "nozzle_expansion_ratio": (nozzle_expansion_ratio_entry, "Nozzle_Expansion_Ratio", "Cleared Nozzle Expansion Ratio"),
        "nozzle_efficiency": (nozzle_efficiency_entry, "Nozzle_Efficiency", "Cleared Nozzle Efficiency"),
        "nozzle_discharge_ratio": (nozzle_discharge_ratio_entry, "Nozzle_Discharge_Ratio", "Cleared Nozzle Discharge Ratio"),
        "c_eff": (c_eff_entry, "c_eff", "Cleared c_eff"),
        "dry_mass": (dry_mass_entry, "dry_mass", "Cleared Dry Mass"),
        "viscosity": (viscosity_entry, "viscosity", "Cleared Viscosity"),
        "for_flight": (for_flight_entry, "For_flight", "Cleared For Flight")
    }

    if var in clear_map:
        entry, variable_name, message = clear_map[var]

        # Handle special case for fuel
        if var == "fuel":
            entry.delete(0, tk.END)
            if selected_option_fuel.get() == "Aluminum":
                Aluminum_weight_percent = 0
                entry.insert(0, str(Aluminum_weight_percent))
            elif selected_option_fuel.get() == "Carbon Black":
                Carbon_black_weight_percent = 0
                entry.insert(0, str(Carbon_black_weight_percent))
        else:
            entry.delete(0, tk.END)
            globals()[variable_name] = 0  # Update the global variable
            entry.insert(0, str(globals()[variable_name]))
        
        print(message)

    return

def sim():
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent, CC_vol, Starting_Tank_Pressure, Starting_Ox_Mass, Starting_Chamber_Pressure, Grain_ID, Grain_OD, Grain_Length, blowing_number, a, n, m, Injector_Hole_Diameter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, c_eff, dry_mass, viscosity, For_flight
    #Ox_tank_vol = 0.01396612489262478177383457064491
    #Ox_tank_length = 1.7526
    #Ox_tank_diameter = ((4*Ox_tank_vol)/(math.pi*Ox_tank_length))**0.5
    #Starting_Tank_Pressure = 5.516e6
    #Starting_Ox_Mass = 18
    
    #Starting_Chamber_Pressure = 101325
    #CC_vol = 0.019661
    

    #Aluminum_weight_percent = 0
    #Carbon_black_weight_percent = 10

    
    # Grain_ID = 0.1
    # Grain_OD = 0.125
    # Grain_Length = 1.5
    # blowing_number = 15
    # a = 0.000155
    # n = 0.45
    # m = 0

    # Assuming showerhead injector for now
    #Injector_Hole_Diameter = 0.0015
    #Number_of_Injector_Holes = 60
    #Injector_Discharge_Coefficient = 0.55
    
    # Nozzle_Throat_Diameter = 0.0954278
    # Nozzle_Expansion_Ratio = 1.2
    # Nozzle_Efficiency = 0.95
    # Nozzle_Discharge_Ratio = 0.9
    # c_eff = 0.9
    
    # dry_mass = 40
    # viscosity = 3.70e-5
    # For_flight = 0
    
    cea.set_global_variables(Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_percent, CC_vol, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, Injector_Hole_Diameter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, c_eff, Grain_ID, Grain_OD, Grain_Length, Starting_Tank_Pressure, Starting_Chamber_Pressure, Starting_Ox_Mass, For_flight, dry_mass, viscosity, blowing_number, a, n, m)
    new_syst = cea.on_button_click()
    AnimRocket.animation(new_syst)

    
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

# Main Title
title1 = ttk.Label(frame, text="This is a Sim Trust", font=("Impact", 18))
title1.grid(row=0, column=0, columnspan=3, pady=10)

#Ox Tank Dimensions & Properties Title
title2 = ttk.Label(frame, text="Ox Tank Dimensions & Properties", font = ("Comic Sans MS", 14))
title2.grid(row=1, column=0, columnspan=3, pady=10)

# Create a tank volume label
tank_vol_label = ttk.Label(frame, text="Tank Volume")
tank_vol_label.grid(row=2, column=0, sticky='w')
# Create a StringVar to track the value of the entry widget
tank_vol_var = tk.StringVar()
# Set up a trace on the StringVar to call on_value_change whenever the value changes
tank_vol_var.trace_add("write", on_value_changed_Ox_tank)
# Create an entry widget
tank_vol_entry = ttk.Entry(frame, width=20, textvariable = tank_vol_var)
tank_vol_entry.grid(row=2, column=1)
# Create a clear button
clear_button1 = ttk.Button(frame, text="Clear", command=lambda: clear("vol"))
clear_button1.grid(row=2, column=2)

# For Tank length
tank_len_label = ttk.Label(frame, text="Tank Length")
tank_len_label.grid(row=3, column=0, sticky='w')
tank_len_var = tk.StringVar()
tank_len_var.trace_add("write", on_value_changed_Ox_tank)
tank_len_entry = ttk.Entry(frame, width=20, textvariable = tank_len_var)
tank_len_entry.grid(row=3, column=1)
clear_button2 = ttk.Button(frame, text="Clear", command=lambda: clear("len"))
clear_button2.grid(row=3, column=2)

# For Tank Diameter
tank_dia_label = ttk.Label(frame, text="Tank Diameter")
tank_dia_label.grid(row=4, column=0, sticky='w')
tank_dia_var = tk.StringVar()
tank_dia_var.trace_add("write", on_value_changed_Ox_tank)
tank_dia_entry = ttk.Entry(frame, width=20, textvariable = tank_dia_var)
tank_dia_entry.grid(row=4, column=1)
clear_button3 = ttk.Button(frame, text="Clear", command=lambda: clear("dia"))
clear_button3.grid(row=4, column=2)

#Tank Calculation Button
compute_button = ttk.Button(frame, text="Find 3rd Var", command=compute_the_third)
compute_button.grid(row=5, column=1, pady=20)

#Tank Pressure
tank_pressure_label = ttk.Label(frame, text="Tank Pressure")
tank_pressure_label.grid(row=6, column=0, sticky='w')
tank_pressure_var = tk.StringVar()
tank_pressure_var.trace_add("write", on_value_changed_tank_pressure)
tank_pressure_entry = ttk.Entry(frame, width=20, textvariable = tank_pressure_var)
tank_pressure_entry.grid(row=6, column=1)
clear_button4 = ttk.Button(frame, text="Clear", command=lambda: clear("tank_pressure"))
clear_button4.grid(row=6, column=2)

#Starting Ox Mass
starting_ox_mass_label = ttk.Label(frame, text="Starting Ox Mass")
starting_ox_mass_label.grid(row=7, column=0, sticky='w')
starting_ox_mass_var = tk.StringVar()
starting_ox_mass_var.trace_add("write", on_value_changed_starting_ox_mass)
starting_ox_mass_entry = ttk.Entry(frame, width=20, textvariable = starting_ox_mass_var)
starting_ox_mass_entry.grid(row=7, column=1)
clear_button5 = ttk.Button(frame, text="Clear", command=lambda: clear("starting_ox_mass"))
clear_button5.grid(row=7, column=2)

#Combustion Chamber Properties Title
title3 = ttk.Label(frame, text="Combustion Chamber Properties", font = ("Comic Sans MS", 14))
title3.grid(row=8, column=0, columnspan=2, pady=(30,10))

#Starting Chamber Pressure
starting_chamber_pressure_label = ttk.Label(frame, text="Starting Chamber Pressure")
starting_chamber_pressure_label.grid(row=9, column=0, sticky='w')
starting_chamber_pressure_var = tk.StringVar()
starting_chamber_pressure_var.trace_add("write", on_value_changed_starting_chamber_pressure)
starting_chamber_pressure_entry = ttk.Entry(frame, width=20, textvariable = starting_chamber_pressure_var)
starting_chamber_pressure_entry.grid(row=9, column=1)
clear_button6 = ttk.Button(frame, text="Clear", command=lambda: clear("starting_chamber_pressure"))
clear_button6.grid(row=9, column=2)

#CC Volume
cc_vol_label = ttk.Label(frame, text="Combustion Chamber Volume")
cc_vol_label.grid(row=10, column=0, sticky='w')
cc_vol_var = tk.StringVar()
cc_vol_var.trace_add("write", on_value_changed_cc_vol)
cc_vol_entry = ttk.Entry(frame, width=20, textvariable = cc_vol_var)
cc_vol_entry.grid(row=10, column=1)
clear_button7 = ttk.Button(frame, text="Clear", command=lambda: clear("cc_volume"))
clear_button7.grid(row=10, column=2)

#Fuel Weight Percent Characterization Title
title3 = ttk.Label(frame, text="Fuel Weight Percent Characterization", font = ("Comic Sans MS", 14))
title3.grid(row=11, column=0, columnspan=2, pady=(30,10))

# Fuel Weight Percent Characterization
weight_percent = tk.StringVar()
weight_percent.trace_add("write", on_value_changed_fuel)
fuel_entry = ttk.Entry(frame, width=20, textvariable=weight_percent)
fuel_entry.grid(row=12, column=1)
clear_button4 = ttk.Button(frame, text="Clear", command=lambda: clear("fuel"))
clear_button4.grid(row=12, column=2)
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
fuel_dropdown.grid(row=12, column=0, sticky='w')

#Grain Dimensions and Characterization Title
title4 = ttk.Label(frame, text="Grain Dimensions and Characterization", font = ("Comic Sans MS", 14))
title4.grid(row=13, column=0, columnspan=2, pady=(30,10))

#Grain ID
grain_id_label = ttk.Label(frame, text="Grain ID")
grain_id_label.grid(row=14, column=0, sticky='w')
grain_id_var = tk.StringVar()
grain_id_var.trace_add("write", on_value_changed_grain_id)
grain_id_entry = ttk.Entry(frame, width=20, textvariable = grain_id_var)
grain_id_entry.grid(row=14, column=1)
clear_button8 = ttk.Button(frame, text="Clear", command=lambda: clear("grain_id"))
clear_button8.grid(row=14, column=2)

#Grain OD
grain_od_label = ttk.Label(frame, text="Grain OD")
grain_od_label.grid(row=15, column=0, sticky='w')
grain_od_var = tk.StringVar()
grain_od_var.trace_add("write", on_value_changed_grain_od)
grain_od_entry = ttk.Entry(frame, width=20, textvariable = grain_od_var)
grain_od_entry.grid(row=15, column=1)
clear_button9 = ttk.Button(frame, text="Clear", command=lambda: clear("grain_od"))
clear_button9.grid(row=15, column=2)

#Grain Length
grain_length_label = ttk.Label(frame, text="Grain Length")
grain_length_label.grid(row=16, column=0, sticky='w')
grain_length_var = tk.StringVar()
grain_length_var.trace_add("write", on_value_changed_grain_length)
grain_length_entry = ttk.Entry(frame, width=20, textvariable = grain_length_var)
grain_length_entry.grid(row=16, column=1)
clear_button10 = ttk.Button(frame, text="Clear", command=lambda: clear("grain_length"))
clear_button10.grid(row=16, column=2)

#Blowing Number
blowing_number_label = ttk.Label(frame, text="Blowing Number")
blowing_number_label.grid(row=17, column=0, sticky='w')
blowing_number_var = tk.StringVar()
blowing_number_var.trace_add("write", on_value_changed_blowing_number)
blowing_number_entry = ttk.Entry(frame, width=20, textvariable = blowing_number_var)
blowing_number_entry.grid(row=17, column=1)
clear_button11 = ttk.Button(frame, text="Clear", command=lambda: clear("blowing_number"))
clear_button11.grid(row=17, column=2)

#a
a_label = ttk.Label(frame, text="a")
a_label.grid(row=18, column=0, sticky='w')
a_var = tk.StringVar()
a_var.trace_add("write", on_value_changed_a)
a_entry = ttk.Entry(frame, width=20, textvariable = a_var)
a_entry.grid(row=18, column=1)
clear_button12 = ttk.Button(frame, text="Clear", command=lambda: clear("a"))
clear_button12.grid(row=18, column=2)

#n
n_label = ttk.Label(frame, text="n")
n_label.grid(row=19, column=0, sticky='w')
n_var = tk.StringVar()
n_var.trace_add("write", on_value_changed_n)
n_entry = ttk.Entry(frame, width=20, textvariable = n_var)
n_entry.grid(row=19, column=1)
clear_button13 = ttk.Button(frame, text="Clear", command=lambda: clear("n"))
clear_button13.grid(row=19, column=2)

#m
m_label = ttk.Label(frame, text="m")
m_label.grid(row=20, column=0, sticky='w')
m_var = tk.StringVar()
m_var.trace_add("write", on_value_changed_m)
m_entry = ttk.Entry(frame, width=20, textvariable = m_var)
m_entry.grid(row=20, column=1)
clear_button14 = ttk.Button(frame, text="Clear", command=lambda: clear("m"))
clear_button14.grid(row=20, column=2)

# Injector Dimensions Title
title5 = ttk.Label(frame, text="Injector Dimensions", font = ("Comic Sans MS", 14))
title5.grid(row=21, column=0, columnspan=2, pady=(30,10))

# Injector Hole Diameter
injector_hole_diameter_label = ttk.Label(frame, text="Injector Hole Diameter")
injector_hole_diameter_label.grid(row=22, column=0, sticky='w')
injector_hole_diameter_var = tk.StringVar()
injector_hole_diameter_var.trace_add("write", on_value_changed_injector_hole_diameter)
injector_hole_diameter_entry = ttk.Entry(frame, width=20, textvariable = injector_hole_diameter_var)
injector_hole_diameter_entry.grid(row=22, column=1)
clear_button15 = ttk.Button(frame, text="Clear", command=lambda: clear("injector_hole_diameter"))
clear_button15.grid(row=22, column=2)

# Number of Injector Holes
number_of_injector_holes_label = ttk.Label(frame, text="Number of Injector Holes")
number_of_injector_holes_label.grid(row=23, column=0, sticky='w')
number_of_injector_holes_var = tk.StringVar()
number_of_injector_holes_var.trace_add("write", on_value_changed_number_of_injector_holes)
number_of_injector_holes_entry = ttk.Entry(frame, width=20, textvariable = number_of_injector_holes_var)
number_of_injector_holes_entry.grid(row=23, column=1)
clear_button16 = ttk.Button(frame, text="Clear", command=lambda: clear("number_of_injector_holes"))
clear_button16.grid(row=23, column=2)

# Injector Discharge Coefficient
injector_discharge_coefficient_label = ttk.Label(frame, text="Injector Discharge Coefficient")
injector_discharge_coefficient_label.grid(row=24, column=0, sticky='w')
injector_discharge_coefficient_var = tk.StringVar()
injector_discharge_coefficient_var.trace_add("write", on_value_changed_injector_discharge_coefficient)
injector_discharge_coefficient_entry = ttk.Entry(frame, width=20, textvariable = injector_discharge_coefficient_var)
injector_discharge_coefficient_entry.grid(row=24, column=1)
clear_button17 = ttk.Button(frame, text="Clear", command=lambda: clear("injector_discharge_coefficient"))
clear_button17.grid(row=24, column=2)

# Nozzle Dimensions and Characterization Title
title6 = ttk.Label(frame, text="Nozzle Dimensions and Characterization", font = ("Comic Sans MS", 14))
title6.grid(row=25, column=0, columnspan=2, pady=(30,10))

# Nozzle Throat Diameter
nozzle_throat_diameter_label = ttk.Label(frame, text="Nozzle Throat Diameter")
nozzle_throat_diameter_label.grid(row=26, column=0, sticky='w')
nozzle_throat_diameter_var = tk.StringVar()
nozzle_throat_diameter_var.trace_add("write", on_value_changed_nozzle_throat_diameter)
nozzle_throat_diameter_entry = ttk.Entry(frame, width=20, textvariable = nozzle_throat_diameter_var)
nozzle_throat_diameter_entry.grid(row=26, column=1)
clear_button18 = ttk.Button(frame, text="Clear", command=lambda: clear("nozzle_throat_diameter"))
clear_button18.grid(row=26, column=2)

# Nozzle Expansion Ratio
nozzle_expansion_ratio_label = ttk.Label(frame, text="Nozzle Expansion Ratio")
nozzle_expansion_ratio_label.grid(row=27, column=0, sticky='w')
nozzle_expansion_ratio_var = tk.StringVar()
nozzle_expansion_ratio_var.trace_add("write", on_value_changed_nozzle_expansion_ratio)
nozzle_expansion_ratio_entry = ttk.Entry(frame, width=20, textvariable = nozzle_expansion_ratio_var)
nozzle_expansion_ratio_entry.grid(row=27, column=1)
clear_button19 = ttk.Button(frame, text="Clear", command=lambda: clear("nozzle_expansion_ratio"))
clear_button19.grid(row=27, column=2)

# Nozzle Efficiency
nozzle_efficiency_label = ttk.Label(frame, text="Nozzle Efficiency")
nozzle_efficiency_label.grid(row=28, column=0, sticky='w')
nozzle_efficiency_var = tk.StringVar()
nozzle_efficiency_var.trace_add("write", on_value_changed_nozzle_efficiency)
nozzle_efficiency_entry = ttk.Entry(frame, width=20, textvariable = nozzle_efficiency_var)
nozzle_efficiency_entry.grid(row=28, column=1)
clear_button20 = ttk.Button(frame, text="Clear", command=lambda: clear("nozzle_efficiency"))
clear_button20.grid(row=28, column=2)

# Nozzle Discharge Ratio
nozzle_discharge_ratio_label = ttk.Label(frame, text="Nozzle Discharge Ratio")
nozzle_discharge_ratio_label.grid(row=29, column=0, sticky='w')
nozzle_discharge_ratio_var = tk.StringVar()
nozzle_discharge_ratio_var.trace_add("write", on_value_changed_nozzle_discharge_ratio)
nozzle_discharge_ratio_entry = ttk.Entry(frame, width=20, textvariable = nozzle_discharge_ratio_var)
nozzle_discharge_ratio_entry.grid(row=29, column=1)
clear_button21 = ttk.Button(frame, text="Clear", command=lambda: clear("nozzle_discharge_ratio"))
clear_button21.grid(row=29, column=2)

# c_eff
c_eff_label = ttk.Label(frame, text="c_eff")
c_eff_label.grid(row=30, column=0, sticky='w')
c_eff_var = tk.StringVar()
c_eff_var.trace_add("write", on_value_changed_c_eff)
c_eff_entry = ttk.Entry(frame, width=20, textvariable = c_eff_var)
c_eff_entry.grid(row=30, column=1)
clear_button22 = ttk.Button(frame, text="Clear", command=lambda: clear("c_eff"))
clear_button22.grid(row=30, column=2)

# Other Properties Title
title7 = ttk.Label(frame, text="Other Properties", font = ("Comic Sans MS", 14))
title7.grid(row=31, column=0, columnspan=2, pady=(30,10))

# Dry Mass
dry_mass_label = ttk.Label(frame, text="Dry Mass")
dry_mass_label.grid(row=32, column=0, sticky='w')
dry_mass_var = tk.StringVar()
dry_mass_var.trace_add("write", on_value_changed_dry_mass)
dry_mass_entry = ttk.Entry(frame, width=20, textvariable = dry_mass_var)
dry_mass_entry.grid(row=32, column=1)
clear_button23 = ttk.Button(frame, text="Clear", command=lambda: clear("dry_mass"))
clear_button23.grid(row=32, column=2)

# Viscosity
viscosity_label = ttk.Label(frame, text="Viscosity")
viscosity_label.grid(row=33, column=0, sticky='w')
viscosity_var = tk.StringVar()
viscosity_var.trace_add("write", on_value_changed_viscosity)
viscosity_entry = ttk.Entry(frame, width=20, textvariable = viscosity_var)
viscosity_entry.grid(row=33, column=1)
clear_button24 = ttk.Button(frame, text="Clear", command=lambda: clear("viscosity"))
clear_button24.grid(row=33, column=2)

# For Flight
for_flight_label = ttk.Label(frame, text="For Flight")
for_flight_label.grid(row=34, column=0, sticky='w')
for_flight_var = tk.StringVar()
for_flight_var.trace_add("write", on_value_changed_for_flight)
for_flight_entry = ttk.Entry(frame, width=20, textvariable = for_flight_var)
for_flight_entry.grid(row=34, column=1)
clear_button25 = ttk.Button(frame, text="Clear", command=lambda: clear("for_flight"))
clear_button25.grid(row=34, column=2)


# Sim Button
compute_button = ttk.Button(frame, text="Click Me", command=sim)
compute_button.grid(row=35, column=1, pady=20)

#Insert Default Value
tank_vol_entry.insert(0, str(Ox_tank_vol))
tank_len_entry.insert(0, str(Ox_tank_length))
tank_dia_entry.insert(0, str(Ox_tank_diameter))
#fuel_entry.insert(0, str(Aluminum_weight_percent))
tank_pressure_entry.insert(0, str(Starting_Tank_Pressure))
starting_ox_mass_entry.insert(0, str(Starting_Ox_Mass))
starting_chamber_pressure_entry.insert(0, str(Starting_Chamber_Pressure))
cc_vol_entry.insert(0, str(CC_vol))
grain_id_entry.insert(0, str(Grain_ID))
grain_od_entry.insert(0, str(Grain_OD))
grain_length_entry.insert(0, str(Grain_Length))
blowing_number_entry.insert(0, str(blowing_number))
a_entry.insert(0, str(a))
n_entry.insert(0, str(n))
m_entry.insert(0, str(m))
injector_hole_diameter_entry.insert(0, str(Injector_Hole_Diameter))
number_of_injector_holes_entry.insert(0, str(Number_of_Injector_Holes))
injector_discharge_coefficient_entry.insert(0, str(Injector_Discharge_Coefficient))
nozzle_throat_diameter_entry.insert(0, str(Nozzle_Throat_Diameter))
nozzle_expansion_ratio_entry.insert(0, str(Nozzle_Expansion_Ratio))
nozzle_efficiency_entry.insert(0, str(Nozzle_Efficiency))
nozzle_discharge_ratio_entry.insert(0, str(Nozzle_Discharge_Ratio))
c_eff_entry.insert(0, str(c_eff))
dry_mass_entry.insert(0, str(dry_mass))
viscosity_entry.insert(0, str(viscosity))
for_flight_entry.insert(0, str(For_flight))

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


'''
def clear(entry, variable, specialcase):
    
    if specialcase == 1:
        entry.delete(0, tk.END)
        if selected_option_fuel.get() == "Aluminum":
            Aluminum_weight_percent = 0
            entry.insert(0, str(Aluminum_weight_percent))
        elif selected_option_fuel.get() == "Carbon Black":
            Carbon_black_weight_percent = 0
            entry.insert(0, str(Carbon_black_weight_percent))
    else:
        entry.delete(0, tk.END)
        variable = 0  # Update the variable
        entry.insert(0, str(variable))
        print(f"Cleared: {str(variable)}")

    return
    
#Helper function to create labeled entries
def create_labeled_entry(parent, row, label_text, variable, var_trace=None):
    label = ttk.Label(parent, text=label_text)
    label.grid(row=row, column=0, pady=10)
    
    var = tk.StringVar()
    if var_trace:
        var.trace_add("write", var_trace)
    
    entry = ttk.Entry(parent, width=20, textvariable=var)
    entry.grid(row=row, column=1, columnspan=2, padx=5, pady=5)
    
    clear_button = ttk.Button(parent, text="Clear", command=lambda: clear(entry, variable, 0))
    clear_button.grid(row=row, column=2, padx=5, pady=5)
    
    return var, entry
    
# Titles
titles = [
    ("This is a Sim Trust","Impact", 18, 0),
    ("Ox Tank Dimensions & Properties","Comic Sans MS", 14, 1),
    ("Ox Tank Dimensions & Properties","Comic Sans MS", 14, 8),
    ("Fuel Weight Percent Characterization","Comic Sans MS", 14, 11),
    ("Grain Dimensions and Characterization","Comic Sans MS", 14, 13),
    ("Grain Dimensions and Characterization","Comic Sans MS", 14, 19)
]

for title, font, size, row in titles:
    ttk.Label(frame, text=title, font=(font, size)).grid(row=row, column=0, columnspan=2, pady=10)

#Entries
entries = []
entries.append(create_labeled_entry(frame, 2, "Tank Volume", Ox_tank_vol, on_value_changed_Ox_tank))
entries.append(create_labeled_entry(frame, 3, "Tank Length", Ox_tank_length, on_value_changed_Ox_tank))
entries.append(create_labeled_entry(frame, 4, "Tank Diameter", Ox_tank_diameter, on_value_changed_Ox_tank))

entries.append(create_labeled_entry(frame, 6, "Tank Pressure", on_value_changed_tank_pressure))
entries.append(create_labeled_entry(frame, 7, "Starting Ox Mass", on_value_changed_starting_ox_mass))
entries.append(create_labeled_entry(frame, 9, "Starting Chamber Pressure", on_value_changed_starting_chamber_pressure))
entries.append(create_labeled_entry(frame, 10, "Combustion Chamber Volume", on_value_changed_cc_vol))
entries.append(create_labeled_entry(frame, 13, "Grain ID", on_value_changed_grain_id))
entries.append(create_labeled_entry(frame, 14, "Grain OD", on_value_changed_grain_od))
entries.append(create_labeled_entry(frame, 15, "Grain Length", on_value_changed_grain_length))
entries.append(create_labeled_entry(frame, 16, "Blowing Number", on_value_changed_blowing_number))
entries.append(create_labeled_entry(frame, 17, "a", on_value_changed_a))
entries.append(create_labeled_entry(frame, 18, "n", on_value_changed_n))
entries.append(create_labeled_entry(frame, 19, "m", on_value_changed_m))

# Fuel Weight Percent Characterization
weight_percent = tk.StringVar()
weight_percent.trace_add("write", on_value_changed_fuel)
fuel_entry = ttk.Entry(frame, width=20, textvariable=weight_percent)
fuel_entry.grid(row=12, column=0, columnspan=2, padx=5, pady=5)
clear_button4 = ttk.Button(frame, text="Clear", command=lambda: clear(fuel_entry,1))
clear_button4.grid(row=12, column=1, padx=5, pady=5)

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
fuel_dropdown.grid(row=12, column=0, pady=10)

#Tank Calculation Button
compute_button = ttk.Button(frame, text="Find 3rd Var", command=compute_the_third)
compute_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

# Sim Button
compute_button = ttk.Button(frame, text="Click Me", command=sim)
compute_button.grid(row=20, column=0, columnspan=2, padx=5)


values = [Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Starting_Tank_Pressure, Starting_Ox_Mass, Starting_Chamber_Pressure, CC_vol, Grain_ID, Grain_OD, Grain_Length, blowing_number, a, n, m]
for i, value in enumerate(values):
    entries[i][1].delete(0, tk.END)
    entries[i][1].insert(0, str(value))

'''
