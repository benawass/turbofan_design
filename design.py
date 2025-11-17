# MY DESIGN

# flight conditions
altitude_m = 11_000 # Assumption, cruise altitude
mach_cruise = 0.85 # Given

# efficencys
diffuser_eff = 0.98
core_compressor_eff = 0.91
fan_eff = 0.93
turb_eff = 0.93
nozzle_eff = 0.95

# pressure ratios
fan_stag_P_ratio = 1.54

# combine Intermediate and high pressure compressor ratios
IPC_P_ratio = 9.61
HPC_P_ratio = 3.38
core_stag_P_ratio = IPC_P_ratio * HPC_P_ratio

combust_stag_P_ratio = 0.95 # NEED BETTER VALUE

# bypass ratio
overall_bypass_ratio = 10

# f value
air_mass_flow_rate = 600
m_core = air_mass_flow_rate / (1 + overall_bypass_ratio)
fuel_mass_flow_rate = 1.116
fuel_to_mass_ratio = 1.116/m_core

# air  & combustion product values (constants)
c_pa = 1.01e3
c_pp = 1.10e3
gamma_a = 1.4
gamma_c= 1.33
Q_r = 45e6 # NEED TO JUSTIFY

# nozzle exit area
core_nozzle_area_m2 = 2.01
# fan inlet diameter
fan_inlet_diameter = 2.85

from turbo_funcs import *
# --- Required Thrust Inputs ---
max_takeoff_mass_kg = 228_000
max_fuel_burn_kg = 60_000
cruise_lift_to_drag_ratio = 21
number_engines = 2

# --- Required Thrust Calculations ---
required_thrust = calculate_engine_cruise_thrust(
    max_takeoff_mass_kg=max_takeoff_mass_kg,
    max_fuel_burn_kg=max_fuel_burn_kg,
    cruise_lift_to_drag_ratio=cruise_lift_to_drag_ratio,
    number_engines=number_engines
)
# --- Display Result ---
print_header("--- Thrust Requirements ---".upper())
print(f"\t{'Required thrust per engine:':<30} {required_thrust:,.2f} N")

# Engine Performance Calcs
all_design_results, thrust_vals, final_TSFC, hpt_targets = calculate_turbofan_engine_performance(
    # flight conditions
    altitude_m,
    mach_cruise,

    # efficencys
    diffuser_eff,
    core_compressor_eff,
    fan_eff,
    turb_eff,
    nozzle_eff,

    # pressure ratios
    fan_stag_P_ratio,
    core_stag_P_ratio,
    combust_stag_P_ratio,

    # f value
    fuel_to_mass_ratio,

    # air  & combustion product values
    c_pa,
    c_pp,
    gamma_a,
    gamma_c,
    Q_r,

    # bypass and nozzle area
    overall_bypass_ratio,
    core_nozzle_area_m2,
    fan_inlet_diameter)

calculate_engine_cruise_thrust(max_takeoff_mass_kg, 
                                 max_fuel_burn_kg, 
                                 cruise_lift_to_drag_ratio, 
                                 number_engines, 
                                 g=9.81)

pretty_print_engine_stats(all_design_results)


print_header(" --- Final Calculated Thrust ---".upper())
print(f"\t{'Total Thrust Per Engine:':<30} {thrust_vals['total_thrust']:,.2f} N")
print(f"\t{'Bypass Thrust Per Engine:':<30} {thrust_vals['bypass_thrust']:,.2f} N")
print(f"\t{'Core Thrust Per Engine:':<30} {thrust_vals['core_thrust']:,.2f} N")

print_header(" --- Final Calculated TSFC ---".upper())
print(f"\t{'Final Calculated TSFC:':<30} {final_TSFC*1e6:,.2f} kg/MN.s")

pretty_print_hpt_targets(hpt_targets)