# EXAMPLE
# flight conditions
altitude_m = 0 # As used in the final get_thrust call
mach_cruise = 0.85

# efficencys
diffuser_eff = 0.9
core_compressor_eff = 0.85
fan_eff = 0.85
turb_eff = 0.88
nozzle_eff = 0.96

# pressure ratios
fan_stag_P_ratio = 2
core_stag_P_ratio = 10
combust_stag_P_ratio = 0.97

# f value
fuel_to_mass_ratio = 0.0160

# air  & combustion product values
c_pa = 1.01e3
c_pp = 1.10e3
gamma_a = 1.4
gamma_c= 1.33
Q_r = 43_000e3

# nozzle exit area
core_nozzle_area_m2 = 0.5

# bypass ratio
overall_bypass_ratio = 1

from turbo_funcs import calculate_turbofan_engine_performance, print_results

all_example_results, final_thrust, final_TSFC = calculate_turbofan_engine_performance(
    # flight conditions
    altitude_m, # As used in the final get_thrust call
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
    core_nozzle_area_m2)

print("All Station Values:")
pretty_print_engine_stats(all_example_results)
print("\nFinal Calculated Thrust:")
print(final_thrust)
print("\nFinal Calculated TSFC:")
print(final_TSFC)