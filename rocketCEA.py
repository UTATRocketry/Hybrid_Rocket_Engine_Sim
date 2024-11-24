from rocketcea.cea_obj import CEA_Obj
#from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
import math
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np
import scipy
import math
import copy
import matplotlib.pyplot as plt
import csv
import tkinter as tk
from tkinter import ttk

#The Sim
#in SI units

#Im setting all the variables global
global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_precent, CC_vol, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, Injector_Hole_Diamter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, c_eff, Grain_ID, Grain_OD, Grain_Length, Starting_Tank_Pressure, Starting_Chamber_Pressure, Starting_Ox_Mass, For_flight
global time_propert, dynamic_system_propert, constant_system_properties, overall_system

def set_global_variables(Ox_tank_vol_var, Ox_tank_length_var, Ox_tank_diameter_var, Aluminum_weight_percent_var, Carbon_black_weight_precent_var, CC_vol_var, Nozzle_Throat_Diameter_var, Nozzle_Expansion_Ratio_var, Nozzle_Efficiency_var, Nozzle_Discharge_Ratio_var, Injector_Hole_Diamter_var, Number_of_Injector_Holes_var, Injector_Discharge_Coefficient_var, c_eff_var, Grain_ID_var, Grain_OD_var, Grain_Length_var, Starting_Tank_Pressure_var, Starting_Chamber_Pressure_var, Starting_Ox_Mass_var, For_flight_var):
    global Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_precent, CC_vol, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, Injector_Hole_Diamter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, c_eff, Grain_ID, Grain_OD, Grain_Length, Starting_Tank_Pressure, Starting_Chamber_Pressure, Starting_Ox_Mass, For_flight
    Ox_tank_vol, Ox_tank_length, Ox_tank_diameter, Aluminum_weight_percent, Carbon_black_weight_precent, CC_vol, Nozzle_Throat_Diameter, Nozzle_Expansion_Ratio, Nozzle_Efficiency, Nozzle_Discharge_Ratio, Injector_Hole_Diamter, Number_of_Injector_Holes, Injector_Discharge_Coefficient, c_eff, Grain_ID, Grain_OD, Grain_Length, Starting_Tank_Pressure, Starting_Chamber_Pressure, Starting_Ox_Mass, For_flight = Ox_tank_vol_var, Ox_tank_length_var, Ox_tank_diameter_var, Aluminum_weight_percent_var, Carbon_black_weight_precent_var, CC_vol_var, Nozzle_Throat_Diameter_var, Nozzle_Expansion_Ratio_var, Nozzle_Efficiency_var, Nozzle_Discharge_Ratio_var, Injector_Hole_Diamter_var, Number_of_Injector_Holes_var, Injector_Discharge_Coefficient_var, c_eff_var, Grain_ID_var, Grain_OD_var, Grain_Length_var, Starting_Tank_Pressure_var, Starting_Chamber_Pressure_var, Starting_Ox_Mass_var, For_flight_var


def Oxidizer_Properties(T, fluid):
    '''Returns a dictionary of properties for a given oxidizer at a specified temperature in Kelvin. 
    The oxidizer name should be as referred to by CoolProp (i.e., NO2, O2). 
    Check this website for reference: 
    http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids'''

    properties = {
        'Temperature': T,
        'Pressure': PropsSI('P', 'T', T, 'Q', 0, fluid),
        'Density_liquid': PropsSI('D', 'T', T, 'Q', 0, fluid),
        'Density_vapor': PropsSI('D', 'T', T, 'Q', 1, fluid),
        'Enthalpy_liquid': PropsSI('H', 'T', T, 'Q', 0, fluid),
        'Enthalpy_vapor': PropsSI('H', 'T', T, 'Q', 1, fluid),
        'Compressibility': PropsSI('Z', 'T', T, 'Q', 0, fluid),
        'Latent_heat_vaporization': PropsSI('H', 'T', T, 'Q', 1, fluid) - PropsSI('H', 'T', T, 'Q', 0, fluid),
        'Cp':PropsSI('C', 'T', T, 'Q', 0, fluid)
    }
    
    return properties

def find_temp_for_vapor_pressure(T, P_target, fluid):
    vapor_pressure = PropsSI('P', 'T', T, 'Q', 0, fluid)
    return vapor_pressure - P_target


def ox_tank(fluid, system_prev, atmospheric_pressure, time_propert, const_propert, total_system_properties):
    current_system = copy.deepcopy(system_prev)
    ox_propert = Oxidizer_Properties(system_prev['Ox_tank_temperature'], fluid)

    current_system['Oxidizer_properties'] = ox_propert
    current_system['P_oxtank'] = current_system['Oxidizer_properties']['Pressure']
    dP = ox_propert['Pressure'] - system_prev['P_chamber']
    R = PropsSI('GAS_CONSTANT', fluid) / PropsSI('M', fluid)
    #Assuming isentropic flow
    M_combust_chamb = (ox_propert['Compressibility']*system_prev['Gamma']*R*(system_prev['Ox_tank_temperature'])*(system_prev['P_chamber']/system_prev['P_oxtank'])**((system_prev['Gamma']-1)/system_prev['Gamma']))**0.5
    M_atmosphere = (ox_propert['Compressibility']*system_prev['Gamma']*R*(system_prev['Ox_tank_temperature'])*(atmospheric_pressure/system_prev['P_oxtank'])**((system_prev['Gamma']-1)/system_prev['Gamma']))**0.5
    # In a closed pipe the max mach number the fluid can reach is 1
    if M_combust_chamb > 1:
        M_combust_chamb = 1
    if M_atmosphere > 1:
        M_atmosphere = 1
    #Pressure difference should be positive
    if dP < 0:
        dP = 0

    if time_propert['end_time'] == 0 or time_propert['Current_time'] <= time_propert['end_time']:
        if system_prev['Current_liquid_oxidizer_mass'] == 0:
            current_system['Mass_Flow_Ox'] = (const_propert['Injector_Coefficient_of_Discharge']*const_propert['Number_of_Holes']*current_system['P_oxtank']/current_system['Oxidizer_properties']['Temperature']**0.5)*(system_prev['Gamma']/(current_system['Oxidizer_properties']['Compressibility']*R))**0.5*M_combust_chamb*(1+(system_prev['Gamma']-1)/2*M_combust_chamb**2)**((-system_prev['Gamma']-1)/(2*(system_prev['Gamma']-1)))
        else:
            current_system['Mass_Flow_Ox'] = (const_propert['Injector_Coefficient_of_Discharge']*math.pi*(const_propert['Injector_Hole_Size']/2)**2*const_propert['Number_of_Holes']*(2*current_system['Oxidizer_properties']['Density_liquid']*dP)**0.5)
    elif time_propert['end_time'] > 0 and time_propert['Current_time'] > time_propert['end_time']:
        current_system['Mass_Flow_Ox'] = 0
    mass_discharged = current_system['Mass_Flow_Ox']*time_propert['Change_in_time']
    current_system['Total_mass_discharged'] = system_prev['Total_mass_discharged'] - mass_discharged
    if system_prev['Current_liquid_oxidizer_mass'] < system_prev['Previous_liquid_oxidizer_mass'] and system_prev['Current_liquid_oxidizer_mass'] > 0 and current_system['Mass_Flow_Ox'] > 0:
        current_system['Previous_liquid_oxidizer_mass'] = current_system['Current_liquid_oxidizer_mass'] - mass_discharged
        current_system['Oxidizer_properties'] = Oxidizer_Properties(current_system['Ox_tank_temperature'], fluid)
        current_system['Current_liquid_oxidizer_mass'] = (constant_system_properties['Ox_tank_volume'] - current_system['Total_mass_discharged']/current_system['Oxidizer_properties']['Density_vapor'])/(1/current_system['Oxidizer_properties']['Density_liquid'] - 1/current_system['Oxidizer_properties']['Density_vapor'])
        mass_of_vapor = current_system['Previous_liquid_oxidizer_mass'] - current_system['Current_liquid_oxidizer_mass']
        dT = -mass_of_vapor*current_system['Oxidizer_properties']['Latent_heat_vaporization']/(current_system['Current_liquid_oxidizer_mass']*current_system['Oxidizer_properties']['Cp'])
        current_system['Ox_tank_temperature'] += dT
        ox_propert = Oxidizer_Properties(current_system['Ox_tank_temperature'], fluid)
        current_system['dP'] = ox_propert['Pressure'] - current_system['P_oxtank']
    elif system_prev['Current_liquid_oxidizer_mass'] >= system_prev['Previous_liquid_oxidizer_mass'] and system_prev['Current_liquid_oxidizer_mass'] > 0 and system_prev['Mass_Flow_Ox'] > 0:
        num_negatives = np.sum(total_system_properties['dP'] < 0)

        if num_negatives > 0:  
            dP_avg = np.mean(total_system_properties['dP'][:num_negatives])
        else:
            dP_avg = None

        P_new = current_system['P_oxtank'] + dP_avg
        current_system['Ox_tank_temperature'], = scipy.optimize.fsolve(find_temp_for_vapor_pressure, 280, args=(P_new, fluid))
        current_system['dP'] = current_system['Oxidizer_properties']['Pressure'] - current_system['P_oxtank']
        current_system['Oxidizer_properties'] = Oxidizer_Properties(current_system['Ox_tank_temperature'], const_propert['fluid'])
        current_system['Current_liquid_oxidizer_mass'] = (constant_system_properties['Ox_tank_volume'] - current_system['Total_mass_discharged']/current_system['Oxidizer_properties']['Density_vapor'])/(1/current_system['Oxidizer_properties']['Density_liquid'] - 1/current_system['Oxidizer_properties']['Density_vapor'])
        current_system['Previous_liquid_oxidizer_mass'] = 0
    elif system_prev['Current_liquid_oxidizer_mass'] <= 0 and current_system['Mass_Flow_Ox'] > 0:
        if current_system['Current_liquid_oxidizer_mass'] != 0:
            current_system['Current_liquid_oxidizer_mass'] = 0
            current_system['Total_mass_discharged'] = 0
            return current_system
        
        Z_old = current_system['Oxidizer_properties']['Compressibility']
        Zguess = Z_old
        epsilon = 1
        tolerance = 1e-6
        T_initial = current_system['Ox_tank_temperature']
        P_initial = current_system['P_oxtank']
        while epsilon >= tolerance:
            T_ratio = ((Zguess * current_system['Total_mass_discharged']) / (Z_old * system_prev['Total_mass_discharged'])) ** 0.3
            T_tnk = T_ratio * T_initial
            P_ratio = T_ratio ** (current_system['Gamma']/ (current_system['Gamma']-1))
            P_tnk = P_ratio * P_initial
            Z = PropsSI('Z', 'T', T_tnk, 'P', P_tnk, fluid)            
            epsilon = abs(Zguess - Z)
            Zguess = (Zguess + Z) / 2
    
    return current_system

def Const_OF(staticsystem, dynamicsystem, time):
    new_dynamic_system = copy.deepcopy(dynamicsystem)
    new_dynamic_system['Mass_flow_fuel'] = new_dynamic_system['Mass_Flow_Ox']/staticsystem['OF']
    new_dynamic_system['Regression_rate'] = new_dynamic_system['Mass_flow_fuel']/(staticsystem['Fuel_density']*math.pi*new_dynamic_system['Grain_ID']*staticsystem['Grain_length'])
    new_dynamic_system['Old_Grain_ID'] = new_dynamic_system['Grain_ID']
    new_dynamic_system['Grain_ID'] = new_dynamic_system['Grain_ID'] + 2*new_dynamic_system['Regression_rate']*time['Change_in_time']
    new_dynamic_system['Fuel_mass'] = new_dynamic_system['Fuel_mass'] - new_dynamic_system['Mass_flow_fuel']*time['Change_in_time']
    return new_dynamic_system



def comb(staticsystem, dynamicsystem, time):
    new_dynamic_system = copy.deepcopy(dynamicsystem)
    if time['Current_time'] <= time['end_time'] or time['end_time'] == 0:

        print(new_dynamic_system['Gamma'])
    return new_dynamic_system


def chamber(staticsystem, dynamicsystem, time, CEA, P_atm):
    if staticsystem['Chamber_volume'] == 0:
        V = 0.25*math.pi*dynamicsystem['Grain_ID']**2*staticsystem['Grain_length']
    else:
        V = staticsystem['Chamber_volume'] - 0.25*math.pi*(dynamicsystem['Grain_OD']**2 - dynamicsystem['Grain_ID']**2)*staticsystem['Grain_length']
    dV = 0.25*math.pi*(dynamicsystem['Grain_ID']**2 - dynamicsystem['Old_Grain_ID']**2)*staticsystem['Grain_length']/time['Change_in_time']
    dynamicsystem['Nozzle_mass_flow'] = dynamicsystem['P_chamber']*staticsystem['Nozzle_discharge_ratio']*0.25*math.pi*staticsystem['Throat_diameter']**2/dynamicsystem['Cstar']
    dm_g = dynamicsystem['Mass_flow_fuel'] + dynamicsystem['Mass_Flow_Ox'] - dynamicsystem['Nozzle_mass_flow']
    if dynamicsystem['Mass_Flow_Ox'] == 0:
        dynamicsystem['dm_g'] = -dynamicsystem['Nozzle_mass_flow']
    
    dynamicsystem['Gas_mass'] = dynamicsystem['Gas_mass'] + dm_g*time['Change_in_time']

    dP = dynamicsystem['P_chamber']*(dm_g/dynamicsystem['Gas_mass'] - dV/V)
    dynamicsystem['P_chamber'] +=dP*time['Change_in_time']
    if dynamicsystem['P_chamber'] <= P_atm:
        dynamicsystem['P_chamber'] = P_atm
        dynamicsystem['Nozzle_mass_flow'] = 0
    return dynamicsystem
    

def sim_iteration(overallsystem, staticsystem, dynamicsystem, time, iteration, CEA, P_atm):
    time['Current_time'] = time['Current_time'] + time['Change_in_time']
    cursystem = ox_tank(staticsystem['fluid'], dynamicsystem, P_atm, time, staticsystem, overallsystem)
    cursystem = Const_OF(staticsystem, cursystem, time)
    cursystem['Cstar'] = CEA.get_Cstar(cursystem['P_chamber'], constant_system_properties['OF'])
    cursystem = chamber(staticsystem, cursystem, time, CEA, P_atm)
    Isp = CEA.get_Isp(cursystem['P_chamber'], staticsystem['OF'], staticsystem['Nozzle_expansion_ratio'])
    overallsystem['time'][iteration] = time['Current_time']
    overallsystem['Total_mass_discharged'][iteration] = cursystem['Total_mass_discharged']
    overallsystem['P_oxtank'][iteration] = cursystem['P_oxtank']
    overallsystem['P_chamber'][iteration] = cursystem['P_chamber']
    overallsystem['Mass_Flow_Ox'][iteration] = cursystem['Mass_Flow_Ox']
    overallsystem['Mass_Flow_Fuel'][iteration] = cursystem['Mass_flow_fuel']
    overallsystem['OF'][iteration] = staticsystem['OF']
    overallsystem['Grain_ID'][iteration] = cursystem['Grain_ID']
    overallsystem['Nozzle_mass_flow'][iteration] = cursystem['Nozzle_mass_flow']
    overallsystem['Regression_rate'][iteration] = cursystem['Regression_rate']
    overallsystem['Fuel_mass'][iteration] = cursystem['Fuel_mass']
    overallsystem['dP'][iteration] = cursystem['dP']
    overallsystem['Isp'][iteration] = Isp*staticsystem['Nozzle_efficiency']
    return staticsystem, cursystem, overallsystem, time

def sim_loop(static_system, dynamic_system, time, overallsystem, CEA):
    i = 0
    new_static_system = copy.deepcopy(static_system)
    new_dynamic_system = copy.deepcopy(dynamic_system)
    new_overall_system = copy.deepcopy(overallsystem)
    ox_mass = new_dynamic_system['Total_mass_discharged']
    while True:
        if new_static_system['is_flying']:
            # To do: allow for flight config
            pass
        else:
            P_atm = 101325
        time['Current_time'] = i*time['Change_in_time']
        i+=1
        new_static_system, new_dynamic_system, new_overall_system, time = sim_iteration(new_overall_system, new_static_system, new_dynamic_system, time, i, CEA, P_atm)
        if new_dynamic_system['Grain_ID']>=new_dynamic_system['Grain_OD']:
            print("No fuel left")
            break
        elif new_dynamic_system['Total_mass_discharged'] <= 0:
            print("No Ox left")
            break
        elif time['Current_time'] >= time['end_time']:
            print("Max Time Reached")
            break
        elif new_dynamic_system['P_chamber'] <= P_atm:
            print("Burn Complete")
            break
    new_overall_system['Impulse'] = new_overall_system['Isp']*ox_mass*(9.8)
    new_overall_system['Mass_Flow_Ox'] = new_overall_system['Mass_Flow_Ox'][new_overall_system['Mass_Flow_Ox']!=0]
    new_overall_system['P_chamber'] = new_overall_system['P_chamber'][new_overall_system['P_chamber']!=0]
    new_overall_system['Impulse'] = new_overall_system['Impulse'][:int(time_propert['end_time']/time_propert['Change_in_time']+1)]
    new_overall_system['Thrust'] = (new_overall_system['Impulse']/time_propert['Current_time'])
    print("Max Thrust (N): ", max(new_overall_system['Thrust']))
    print("Average Thrust (N): ", np.average(new_overall_system['Thrust']))
    print("Average Mass Flow Rate (kg/s): ", np.average(new_overall_system['Mass_Flow_Ox']))
    print("Average Combustion Chamber Pressure (psi): ", np.average(new_overall_system['P_chamber'])*145/10e5)
    print("Max Combustion Chamber Pressure (psi): ", max(new_overall_system['P_chamber'])*145/10e5)
    print("Burn Time (s): ", max(new_overall_system['time']))
    visualize(new_overall_system)
    

def visualize(overallsystem, filename="output.csv"):
    # Find the maximum length of all lists in the dictionary
    max_length = max(len(overallsystem[key]) for key in overallsystem)

    # Extend all lists to the maximum length with zeros
    for key in overallsystem:
        if len(overallsystem[key]) < max_length:
            overallsystem[key] = np.append(overallsystem[key], [0] * (max_length - len(overallsystem[key])))
    
    # Find the index where time becomes 0 again
    time_array = overallsystem['time']
    zero_index = np.where(time_array == 0)[0]
    
    if len(zero_index) > 1:
        truncate_index = zero_index[1]  # Second occurrence of 0
    else:
        truncate_index = len(time_array)  # No second occurrence, no truncation needed

    # Truncate all lists in the dictionary at the truncate_index
    for key in overallsystem:
        overallsystem[key] = overallsystem[key][:truncate_index]
    
    headers = [
        'Time (s)', 
        'Total Mass Discharged (kg)', 
        'Oxidizer Tank Pressure (Pa)', 
        'Chamber Pressure (Pa)', 
        'Mass Flow Rate Oxidizer (kg/s)', 
        'Mass Flow Rate Fuel (kg/s)', 
        'O/F Ratio', 
        'Grain ID (m)', 
        'Nozzle Mass Flow Rate ()', 
        'Regression Rate (m/s)', 
        'Fuel Mass (kg)', 
        'dP (Pa)', 
        'Isp (s)'
    ]
    # Write the data to a CSV file
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for i in range(len(overallsystem['time'])):
            writer.writerow([
                overallsystem['time'][i], 
                overallsystem['Total_mass_discharged'][i], 
                overallsystem['P_oxtank'][i], 
                overallsystem['P_chamber'][i], 
                overallsystem['Mass_Flow_Ox'][i], 
                overallsystem['Mass_Flow_Fuel'][i], 
                overallsystem['OF'][i], 
                overallsystem['Grain_ID'][i], 
                overallsystem['Nozzle_mass_flow'][i], 
                overallsystem['Regression_rate'][i], 
                overallsystem['Fuel_mass'][i], 
                overallsystem['dP'][i], 
                overallsystem['Isp'][i]
            ])

    # Plot the data
    fig, axs = plt.subplots(6, 2, figsize=(15, 20))  # 6 rows, 2 columns

    axs[0, 0].plot(overallsystem['time'], overallsystem['Total_mass_discharged'])
    axs[0, 0].set_title('Total Mass Discharged (kg)')
        
    axs[0, 1].plot(overallsystem['time'], overallsystem['P_oxtank'])
    axs[0, 1].set_title('Oxidizer Tank Pressure (Pa)')

    axs[1, 0].plot(overallsystem['time'], overallsystem['P_chamber'])
    axs[1, 0].set_title('Chamber Pressure (Pa)')

    axs[1, 1].plot(overallsystem['time'], overallsystem['Mass_Flow_Ox'])
    axs[1, 1].set_title('Mass Flow Rate Oxidizer (kg/s)')

    axs[2, 0].plot(overallsystem['time'], overallsystem['Mass_Flow_Fuel'])
    axs[2, 0].set_title('Mass Flow Rate Fuel (kg/s)')

    axs[2, 1].plot(overallsystem['time'], overallsystem['OF'])
    axs[2, 1].set_title('O/F Ratio')

    axs[3, 0].plot(overallsystem['time'], overallsystem['Grain_ID'])
    axs[3, 0].set_title('Grain ID (m)')

    axs[3, 1].plot(overallsystem['time'], overallsystem['Nozzle_mass_flow'])
    axs[3, 1].set_title('Nozzle Mass Flow Rate ()')
       
    axs[4, 0].plot(overallsystem['time'], overallsystem['Regression_rate'])
    axs[4, 0].set_title('Regression Rate (m/s)')

    axs[4, 1].plot(overallsystem['time'], overallsystem['Fuel_mass'])
    axs[4, 1].set_title('Fuel Mass (kg)')

    axs[5, 0].plot(overallsystem['time'], overallsystem['dP'])
    axs[5, 0].set_title('dP (Pa)')

    axs[5, 1].plot(overallsystem['time'], overallsystem['Isp'])
    axs[5, 1].set_title('Isp (s)')

    # Adjust layout
    plt.tight_layout()
    plt.show()

def on_button_click():
    
    global time_propert, dynamic_system_propert, constant_system_properties, overall_system
    '''
    Ox_tank_vol = 0.01396612489262478177383457064491
    Ox_tank_length = 1.7526
    Ox_tank_diameter = (Ox_tank_vol/(math.pi*Ox_tank_length))**0.5
    Aluminum_weight_percent = 0
    Carbon_black_weight_precent = 10


    CC_vol = 0.4826

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
    '''
    #Assuming constant specific heat capacity for now
    gamma = 1.31
    fluid = 'N2O'

    time_step = 0.01
    # I Want to try and avoid using this
    simulation_time = 8
    OF_ratio = 4.5
    fuel_density = 1000

    card_str = """
    fuel
    fuel C C(s)       C 1.0      wt%={0}
    t(k)=298.15       h,cal=0.0

    fuel AL           AL 1.0     wt%={1}
    t(k)=298.15       h,cal=0.0

    fuel C22H46       C 22.0  H 46.0     wt%={2}
    t(k)=298.15       h,cal=-39600.0   rho=0.9
    """.format(Carbon_black_weight_precent, Aluminum_weight_percent, 100 - Carbon_black_weight_precent - Aluminum_weight_percent)
    add_new_fuel( 'Paraffin', card_str )
    C = CEA_Obj(fuelName="Paraffin", oxName='N2O')
    time_propert={
        'Current_time':0,
        'Change_in_time':time_step,
        'end_time':simulation_time
    }
    dynamic_system_propert={
        'Ox_tank_temperature':PropsSI('T', 'P', Starting_Tank_Pressure, 'Q', 1, fluid),
        'P_chamber': Starting_Chamber_Pressure,
        'P_oxtank': Starting_Tank_Pressure,
        'Gamma':gamma,
        'Total_mass_discharged':Starting_Ox_Mass,
        'Grain_ID':Grain_ID,
        'Grain_OD':Grain_OD,
        'Fuel_mass':Starting_Ox_Mass/OF_ratio,
        'Gas_mass':0,
    }
    constant_system_properties={
        'Number_of_Holes':Number_of_Injector_Holes,
        'Injector_Hole_Size':Injector_Hole_Diamter,
        'Injector_Coefficient_of_Discharge':Injector_Discharge_Coefficient,
        'Ox_tank_volume':Ox_tank_vol,
        'fluid':fluid,
        'is_flying':False,
        'OF':OF_ratio,
        'Fuel_density':fuel_density,
        'Grain_length':Grain_Length,
        'Critical_pressure':PropsSI('Pcrit', fluid),
        'Chamber_volume':CC_vol, 
        'Nozzle_discharge_ratio':Nozzle_Discharge_Ratio,
        'Throat_diameter':Nozzle_Throat_Diameter,
        'Nozzle_expansion_ratio':Nozzle_Expansion_Ratio,
        'Nozzle_efficiency':Nozzle_Efficiency
    }
    overall_system = {
        'time':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Total_mass_discharged':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'P_oxtank':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'P_chamber':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Mass_Flow_Ox':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Mass_Flow_Fuel':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'OF':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Grain_ID':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Nozzle_mass_flow':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Fuel_mass':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'dP':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Thrust':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Regression_rate':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
        'Isp':np.zeros(int(time_propert['end_time']/time_propert['Change_in_time']+1)),
    }
    dynamic_system_propert['Oxidizer_properties'] = Oxidizer_Properties(dynamic_system_propert['Ox_tank_temperature'], fluid)
    dynamic_system_propert['Current_liquid_oxidizer_mass'] = (Ox_tank_vol - Starting_Ox_Mass/dynamic_system_propert['Oxidizer_properties']['Density_vapor'])/(1/dynamic_system_propert['Oxidizer_properties']['Density_liquid'] - 1/dynamic_system_propert['Oxidizer_properties']['Density_vapor'])
    dynamic_system_propert['Previous_liquid_oxidizer_mass'] = dynamic_system_propert['Current_liquid_oxidizer_mass'] + 1
    sim_loop(constant_system_properties, dynamic_system_propert, time_propert, overall_system, C)