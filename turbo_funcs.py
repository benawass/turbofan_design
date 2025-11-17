import numpy as np # Import numpy for sqrt

import numpy as np

def calculate_engine_cruise_thrust(max_takeoff_mass_kg, 
                                 max_fuel_burn_kg, 
                                 cruise_lift_to_drag_ratio, 
                                 number_engines, 
                                 g=9.81):
    """
    Calculates the required thrust per engine at mid-cruise conditions.

    The calculation assumes steady, level flight (Lift=Weight, Thrust=Drag)
    at the mid-cruise mass, which is the average of the takeoff
    and end-of-cruise (landing) mass.

    Args:
        max_takeoff_mass_kg (float): Maximum takeoff mass (MTOM) in kg.
        max_fuel_burn_kg (float): Total fuel mass to be burned in kg.
        cruise_lift_to_drag_ratio (float): The lift-to-drag ratio (L/D) at cruise.
        number_engines (int): The number of engines on the aircraft.
        g (float, optional): Acceleration due to gravity in m/s^2. 
                             Defaults to 9.81.

    Returns:
        float: The required thrust per engine in Newtons (N).
    """
    
    # Find design point conditions
    end_mass = max_takeoff_mass_kg - max_fuel_burn_kg
    mid_cruise_mass_kg = (max_takeoff_mass_kg + end_mass) / 2

    # Find net cruise thrust
    # In level, steady flight: Lift = Weight, Thrust = Drag
    mid_cruise_weight_N = mid_cruise_mass_kg * g

    # Use L/D to find total thrust required: T = W / (L/D)
    total_cruise_thrust_N = mid_cruise_weight_N / cruise_lift_to_drag_ratio

    # Distribute total thrust among engines
    engine_cruise_thrust_N = total_cruise_thrust_N / number_engines
    
    return engine_cruise_thrust_N

def get_atmospheric_properties(height_m):
    from ambiance import Atmosphere

    # Create an Atmosphere instance
    atm = Atmosphere(height_m)

    # Access static atmospheric properties
    P = atm.pressure[0] # Pascals
    T = atm.temperature[0] # Kelvin
    rho = atm.density[0] # kg/m^3

    # Calculate speed of sound
    gamma = 1.4  # Ratio of specific heats for air
    R = 287.05   # Gas constant for dry air (J/(kg*K))
    a = np.sqrt(gamma * R * T)

    return {'P': P, 'T': T, 'rho': rho, 'a': a}


def get_intake_values(mach_number, diffuser_eff, P_atm, T_atm, gamma_a):
  # find real and ideal stagnation temps
  T_01 = T_atm*(1 + ((gamma_a-1)/2)*mach_number**2)
  T_01s = T_atm + diffuser_eff*(T_01 - T_atm)

  # find real and ideal stagnation pressures
  P_01s = P_atm*(T_01s/T_atm)**(gamma_a/(gamma_a-1))
  P_01 = P_01s
  return {'T_01': T_01, 'P_01': P_01, 'T_01s': T_01s, 'P_01s': P_01s}


def get_fan_exit_values(fan_stag_P_ratio, fan_eff, P_01, T_01, gamma_a):
  # ideal and real stag pressures
  P_02s = fan_stag_P_ratio*P_01
  P_02 = P_02s

  # ideal and real stag temps
  T_02s = T_01 * fan_stag_P_ratio**((gamma_a-1)/gamma_a)
  T_02 = T_01 + (1/fan_eff)*(T_02s-T_01)
  return {'T_02': T_02, 'P_02': P_02, 'T_02s': T_02s, 'P_02s': P_02s}


def get_compressor_exit_values(core_stag_P_ratio, core_compressor_eff, P_02, T_02, gamma_a):
  # ideal and real stag pressures
  P_03s = core_stag_P_ratio*P_02
  P_03 = P_03s

  # ideal and real stag temps
  T_03s = T_02 * core_stag_P_ratio**((gamma_a-1)/gamma_a)
  T_03 = T_02 + (1/core_compressor_eff)*(T_03s-T_02)
  return {'T_03': T_03, 'P_03': P_03, 'T_03s': T_03s, 'P_03s': P_03s}


def get_combustor_exit_values(c_pa, c_pp, fuel_to_mass_ratio, Q_r, combust_stag_P_ratio, T_03, P_03):
  # energy balance for fuel burnt
  T_04 = ((c_pa*T_03) + (fuel_to_mass_ratio*Q_r))/((1+fuel_to_mass_ratio)*c_pp)

  # combust pressure ratio to find pressure
  P_04 = P_03 * combust_stag_P_ratio

  return {'T_04': T_04, 'P_04': P_04}

# Note: USE OF GAMMA_C as has been combusted already
def get_turbine_exit_values(bypass_ratio, c_pa, c_pp, turb_eff, all_station_values, gamma_c):
  # extract relevant values
  T_01 = all_station_values['T_01']
  T_02 = all_station_values['T_02']
  T_03 = all_station_values['T_03']
  T_04 = all_station_values['T_04']
  P_04 = all_station_values['P_04']

  # find T deltas
  delta_T_compressor = T_03 - T_02
  delta_T_fan = T_02 - T_01

  # find temp change in turbine
  delta_T_turbine = ((c_pa*delta_T_compressor) + (1+bypass_ratio)*(c_pa*delta_T_fan))/c_pp

  # find ideal and real temps at exit of turbine
  T_05 = T_04 - delta_T_turbine
  T_05s = T_04 - (delta_T_turbine/turb_eff)

  # find ideal and real pressures at exit of turbine
  P_05 = P_04 * (T_05s / T_04) ** (gamma_c / (gamma_c - 1))
  return {'T_05': T_05, 'P_05': P_05, 'T_05s': T_05s,}

def get_choked_core_nozzle_exit_values(all_station_values, gamma_c, nozzle_eff, c_pp, R_c, core_nozzle_area_m2):
  # get relevant values
  T_05 = all_station_values['T_05']
  P_05 = all_station_values['P_05']
  T_06 = T_05

  # find temps
  T_6 = T_06/((gamma_c+1)/2)
  T_6s = T_05 - ((T_05-T_6)/nozzle_eff)

  # find speed via SFEE
  u_6 = np.sqrt(2*c_pp*(T_06-T_6))

  # find pressures (isen temp/pressure relationship)
  P_6 = P_05 * (T_6s/T_05)**(gamma_c/(gamma_c-1))

  # find density at exit using ideal gas law
  rho_6 = P_6/(R_c*T_6);

  return {'T_6': T_6, 'T_06':T_06, 'u_6': u_6, 'T_6s': T_6s, 'P_6': P_6, 'rho_6': rho_6}


def get_choked_bypass_nozzle_exit_values(all_station_values, gamma_a, c_pa, R_a, bypass_ratio, nozzle_eff):
  # get relevant temps
  T_02 = all_station_values['T_02']
  T_06 = all_station_values['T_06']
  mass_flow_bypass = all_station_values['mass_flow_bypass']
  P_02 = all_station_values['P_02']

  # calculate temps
  T_07 = T_02
  T_7 = T_07/((gamma_a+1)/2) # using gamma a here as BYPASSES COMBUSTION
  T_7s = T_07 - ((T_07-T_7)/nozzle_eff)

  # calculate speed
  u_7 = np.sqrt(2*c_pa*(T_07-T_7))

  # find pressure
  P_7s = P_02 * (T_7s/T_02)**(gamma_a/(gamma_a-1))
  P_7 = P_7s # nozzle therefore no work done

  # find density
  rho_7 = P_7/(R_a*T_7)

  # find bypass area
  bypass_area_m2 = mass_flow_bypass / (rho_7 * u_7)

  # return vars
  return {'T_7': T_7, 'T_7s': T_7s, 'T_07':T_07, 'u_7': u_7, 'P_7':P_7, 'P_7s': P_7s, 'rho_7':rho_7, 'bypass_area_m2':bypass_area_m2}


def get_thrust(all_station_values, mach_cruise, altitude, core_nozzle_area_m2):
  # mass and area of core & bypass
  mass_flow_core = all_station_values['mass_flow_core']
  mass_flow_bypass = all_station_values['mass_flow_bypass']
  bypass_area_m2 = all_station_values['bypass_area_m2']

  # extract speeds out of core and bypass
  u_6 = all_station_values['u_6']
  u_7 = all_station_values['u_7']

  # find inlet speed of sound at altitude
  atm_values = get_atmospheric_properties(altitude)
  a_1 = atm_values['a']
  u_1 = a_1 * mach_cruise

  # pressure at bypass and core exits
  P_6 = all_station_values['P_6']
  P_7 = all_station_values['P_7']
  P_atm = all_station_values['P_atm']

  core_thrust = mass_flow_core*(u_6 - u_1) + (P_6 - P_atm)*core_nozzle_area_m2
  bypass_thrust = mass_flow_bypass*(u_7 - u_1) + (P_7 - P_atm)*bypass_area_m2

#   print(f'Core Thrust: {core_thrust} N')
#   print(f'Bypass Thrust: {bypass_thrust} N')
  total_thrust = core_thrust + bypass_thrust
  return {'core_thrust': core_thrust, 'bypass_thrust': bypass_thrust, 'total_thrust': total_thrust}


def get_TSFC(f, all_station_values, thrust):
  TSFC =(f * all_station_values['mass_flow_core'])/thrust
  return TSFC

def calculate_mass_flows(all_station_values, overall_bypass_ratio, fan_inlet_diameter, mach_cruise):
    # get values
    fan_inlet_area = np.pi * (fan_inlet_diameter/2)**2
    rho_atm = all_station_values['rho_atm']
    aircraft_speed = mach_cruise * all_station_values['a_atm']

    # total mass flow
    mass_flow_total = rho_atm * fan_inlet_area * aircraft_speed

    # distribute mass flows
    mass_flow_core = mass_flow_total/(1+overall_bypass_ratio)
    mass_flow_bypass = mass_flow_total * (overall_bypass_ratio/(1+overall_bypass_ratio))

    mass_flow_values = {
       'mass_flow_total': mass_flow_total,
       'mass_flow_core': mass_flow_core,
       'mass_flow_bypass': mass_flow_bypass
    }

    return mass_flow_values


def calculate_turbofan_engine_performance(
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
    core_nozzle_area_m2,
    fan_inlet_diameter):

    # Initialize an empty dictionary for all station values
    all_station_values = {}

    # 3. Call get_atmospheric_properties and store results
    atm_properties = get_atmospheric_properties(altitude_m)
    all_station_values['P_atm'] = atm_properties['P']
    all_station_values['T_atm'] = atm_properties['T']
    all_station_values['rho_atm'] = atm_properties['rho']
    all_station_values['a_atm'] = atm_properties['a']

    # 3.5 Get mass flow values
    mass_flow_vals = calculate_mass_flows(all_station_values, overall_bypass_ratio, fan_inlet_diameter, mach_cruise)
    all_station_values.update(mass_flow_vals)

    # 4. Call get_intake_values
    intake_values = get_intake_values(
        mach_cruise,
        diffuser_eff,
        all_station_values['P_atm'],
        all_station_values['T_atm'],
        gamma_a
    )
    all_station_values.update(intake_values)

    # 5. Call get_fan_exit_values
    fan_exit_values = get_fan_exit_values(
        fan_stag_P_ratio,
        fan_eff,
        all_station_values['P_01'],
        all_station_values['T_01'],
        gamma_a
    )
    all_station_values.update(fan_exit_values)

    # 6. Call compressor_exit_values
    compressor_exit_results = get_compressor_exit_values(
        core_stag_P_ratio,
        core_compressor_eff,
        all_station_values['P_02'],
        all_station_values['T_02'],
        gamma_a
    )
    all_station_values.update(compressor_exit_results)

    # 7. Call get_combustor_exit_values
    combustor_exit_results = get_combustor_exit_values(
        c_pa,
        c_pp,
        fuel_to_mass_ratio,
        Q_r,
        combust_stag_P_ratio,
        all_station_values['T_03'],
        all_station_values['P_03']
    )
    all_station_values.update(combustor_exit_results)

    # 8. Calculate R_c and R_a within the function
    R_c = c_pp * ((gamma_c - 1) / gamma_c)
    R_a = c_pa * ((gamma_a - 1) / gamma_a)
    all_station_values['R_c'] = R_c
    all_station_values['R_a'] = R_a

    # 9. Call get_turbine_exit_values
    turbine_exit_results = get_turbine_exit_values(
        overall_bypass_ratio,
        c_pa,
        c_pp,
        turb_eff,
        all_station_values,
        gamma_c
    )
    all_station_values.update(turbine_exit_results)

    # 10. Call get_choked_core_nozzle_exit_values
    core_nozzle_results = get_choked_core_nozzle_exit_values(
        all_station_values,
        gamma_c,
        nozzle_eff,
        c_pp,
        R_c,
        core_nozzle_area_m2
    )
    all_station_values.update(core_nozzle_results)

    # 11. Call get_choked_bypass_nozzle_exit_values
    bypass_nozzle_results = get_choked_bypass_nozzle_exit_values(
        all_station_values,
        gamma_a,
        c_pa,
        R_a,
        overall_bypass_ratio,
        nozzle_eff
    )
    all_station_values.update(bypass_nozzle_results)

    # 12. Call get_thrust
    thrust_values = get_thrust(
        all_station_values,
        mach_cruise,
        altitude_m,
        core_nozzle_area_m2
    )

    # 13. Call get_TSFC
    total_thrust = thrust_values['total_thrust']
    TSFC = get_TSFC(fuel_to_mass_ratio, all_station_values, total_thrust)

    # 14. Return the dictionary and final_thrust
    return all_station_values, thrust_values, TSFC

def print_header(title, char="="):
    """Helper function to print a formatted header."""
    width = 60
    print(f"\n{title.center(width, ' ')}")
    print(char * width)


def pretty_print_engine_stats(stats):
    """
    Prints a formatted report of engine station values from a dictionary.

    Args:
        stats (dict): A dictionary containing all station values, mass flows,
                      and thrust calculations.
    """
    
    print_header("--- TURBOFAN PERFORMANCE REPORT ---", char="*")

    # --- Ambient Conditions ---
    print_header("Ambient Conditions (Station 0)", char="-")
    print(f"\t{'Ambient Pressure (P_atm):':<30} {stats.get('P_atm', 0.0):,.2f} Pa")
    print(f"\t{'Ambient Temperature (T_atm):':<30} {stats.get('T_atm', 0.0):,.2f} K")
    print(f"\t{'Ambient Density (rho_atm):':<30} {stats.get('rho_atm', 0.0):,.4f} kg/m³")
    print(f"\t{'Speed of Sound (a_atm):':<30} {stats.get('a_atm', 0.0):,.2f} m/s")

    # --- Mass Flow and Areas ---
    print_header("Mass Flow & Areas", char="-")
    print(f"\t{'Total Mass Flow:':<30} {stats.get('mass_flow_total', 0.0):,.2f} kg/s")
    print(f"\t{'Core Mass Flow:':<30} {stats.get('mass_flow_core', 0.0):,.2f} kg/s")
    print(f"\t{'Bypass Mass Flow:':<30} {stats.get('mass_flow_bypass', 0.0):,.2f} kg/s")
    if 'core_area_m2' in stats:
        print(f"\t{'Core Nozzle Area:':<30} {stats.get('core_area_m2', 0.0):,.4f} m²")
    if 'bypass_area_m2' in stats:
        print(f"\t{'Bypass Nozzle Area:':<30} {stats.get('bypass_area_m2', 0.0):,.4f} m²")

    # --- Thrust ---
    if 'total_thrust' in stats or 'core_thrust' in stats:
        print_header("Thrust Output", char="-")
    if 'core_thrust' in stats:
        print(f"\t{'Core Thrust:':<30} {stats.get('core_thrust', 0.0):,.2f} N")
    if 'bypass_thrust' in stats:
        print(f"\t{'Bypass Thrust:':<30} {stats.get('bypass_thrust', 0.0):,.2f} N")
    if 'total_thrust' in stats:
        total_thrust = stats.get('total_thrust', 0.0)
        print(f"\t{'TOTAL THRUST:':<30} {total_thrust:,.2f} N  ({total_thrust/1000:,.2f} kN)")
    if 'tsfc' in stats:
         print(f"\t{'TSFC:':<30} {stats.get('tsfc', 0.0):.2e} kg/(N·s)")

    # --- Station 1: Fan Inlet ---
    print_header("Station 1 (Fan Inlet)")
    print(f"\t{'Total Pressure (P_01):':<30} {stats.get('P_01', 0.0):,.2f} Pa")
    print(f"\t{'Total Temperature (T_01):':<30} {stats.get('T_01', 0.0):,.2f} K")
    if 'P_01s' in stats:
        print(f"\t{'Isentropic Temp (T_01s):':<30} {stats.get('T_01s', 0.0):,.2f} K")

    # --- Station 2: Fan Exit / Core Inlet ---
    print_header("Station 2 (Fan Exit / Core Inlet)")
    print(f"\t{'Total Pressure (P_02):':<30} {stats.get('P_02', 0.0):,.2f} Pa")
    print(f"\t{'Total Temperature (T_02):':<30} {stats.get('T_02', 0.0):,.2f} K")
    if 'P_02s' in stats:
        print(f"\t{'Isentropic Temp (T_02s):':<30} {stats.get('T_02s', 0.0):,.2f} K")

    # --- Station 3: Core Compressor Exit ---
    print_header("Station 3 (Core Compressor Exit)")
    print(f"\t{'Total Pressure (P_03):':<30} {stats.get('P_03', 0.0):,.2f} Pa")
    print(f"\t{'Total Temperature (T_03):':<30} {stats.get('T_03', 0.0):,.2f} K")
    if 'P_03s' in stats:
        print(f"\t{'Isentropic Temp (T_03s):':<30} {stats.get('T_03s', 0.0):,.2f} K")

    # --- Station 4: Combustor Exit / Turbine Inlet ---
    print_header("Station 4 (Combustor Exit / Turbine Inlet)")
    print(f"\t{'Total Pressure (P_04):':<30} {stats.get('P_04', 0.0):,.2f} Pa")
    print(f"\t{'Total Temperature (T_04):':<30} {stats.get('T_04', 0.0):,.2f} K")

    # --- Station 5: Turbine Exit ---
    print_header("Station 5 (Turbine Exit)")
    print(f"\t{'Total Pressure (P_05):':<30} {stats.get('P_05', 0.0):,.2f} Pa")
    print(f"\t{'Total Temperature (T_05):':<30} {stats.get('T_05', 0.0):,.2f} K")
    if 'P_05s' in stats:
        print(f"\t{'Isentropic Temp (T_05s):':<30} {stats.get('T_05s', 0.0):,.2f} K")

    # --- Station 6: Core Nozzle Exit ---
    print_header("Station 6 (Core Nozzle Exit)")
    print(f"\t{'Static Pressure (P_6):':<30} {stats.get('P_6', 0.0):,.2f} Pa")
    print(f"\t{'Static Temperature (T_6):':<30} {stats.get('T_6', 0.0):,.2f} K")
    print(f"\t{'Static Density (rho_6):':<30} {stats.get('rho_6', 0.0):,.4f} kg/m³")
    print(f"\t{'Exit Velocity (u_6):':<30} {stats.get('u_6', 0.0):,.2f} m/s")
    
    # --- Station 7: Fan Nozzle Exit ---
    print_header("Station 7 (Fan Nozzle Exit)")
    print(f"\t{'Static Pressure (P_7):':<30} {stats.get('P_7', 0.0):,.2f} Pa")
    print(f"\t{'Static Temperature (T_7):':<30} {stats.get('T_7', 0.0):,.2f} K")
    print(f"\t{'Static Density (rho_7):':<30} {stats.get('rho_7', 0.0):,.4f} kg/m³")
    print(f"\t{'Exit Velocity (u_7):':<30} {stats.get('u_7', 0.0):,.2f} m/s")

    # --- Gas Properties ---
    if 'R_a' in stats or 'R_c' in stats:
        print_header("Gas Properties", char="-")
    if 'R_a' in stats:
        print(f"\t{'R (Air):':<30} {stats.get('R_a', 0.0):,.2f} J/(kg·K)")
    if 'R_c' in stats:
        print(f"\t{'R (Combustion):':<30} {stats.get('R_c', 0.0):,.2f} J/(kg·K)")

    print("\n" + "*" * 60)
