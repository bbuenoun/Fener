# -*- coding: cp1252 -*-
from __future__ import division
from scipy.optimize import fsolve
from math import exp

# import matplotlib.pyplot as plt
# import pandas as pd
import os

the_folder = os.path.dirname(__file__) + "/"

###############################################################################
# Declarations of constants:

SB = 5.6693e-8  # [SI] Stefan Boltzmann constant
T_ref = 273.15  # Conversion °C --> K
g = 9.81  # [m/s²] Gravitationnal constant
R = 8.314462  #

###############################################################################
# Definition of a class with all the inputs


class Inputs:
    def __init__(
        self,
        N_layers,
        air_r,
        argon_r,
        krypton_r,
        xenon_r,
        lambdaa,
        t,
        d,
        eps_ft,
        eps_b,
        T,
        re_ft,
        re_b,
        H,
        w,
        S_abs,
        airflow_r,
        T_exterior,
        T_interior,
        T_inlet,
        h_cv_e,
        h_cv_i,
        E_gv_exterior,
        E_gv_interior,
    ):

        # Declaration and initialization of parameters:

        self.N_layers = (
            N_layers  # The number of solid layers (including blinds if down)
        )
        # Layer 0 is the exterior, layer N+1 the interior !!!
        self.air_ratio = air_r
        self.argon_ratio = argon_r
        self.krypton_ratio = krypton_r
        self.xenon_ratio = xenon_r
        self.lambda_gv = (
            lambdaa  # [W/(m*K)] The thermal conductivity of the solid layers
        )
        self.t_gv = t  # [m] Thickness of the solid layers
        self.d_gv = d  # [m] Thickness of the air layers d_gv[i] is the thickness of the air layer between layer i and i-1
        self.Epsilon_ft = (
            eps_ft  # [-] Longwave emissivity of the front of the solid layer
        )
        self.Epsilon_b = (
            eps_b  # [-] Longwave emissivity of the back of the solid layers
        )
        self.Tau = T  # [-] Longwave transmissivity of solid layers
        self.r_ft = re_ft  # [-] Longwave reflectivity of the front of the solid layers
        self.r_b = re_b  # [-] Longwave reflectivity of the front of the solid layers
        self.height_glazed_area = H  # [m] the height of the glazed area
        self.width_glazed_area = w  # [m] the width of the glazed area.

        self.S = S_abs  # [W/mÂ²] Energy absorbed by each layer
        self.airflow_rates = airflow_r  # [m³/h] Airflow rate of the corresponding  air gaps. airflow_rates[i] correspond to the airflow between layers i and i-1.

        self.T_ai_ex = T_exterior  # [K] Temperature of the outside air
        self.T_ai_int = T_interior  # [K] Temperature of the inside air
        self.T_air_inlet = T_inlet  # [K] Temperature of the air enterring in the gaps

        self.h_cv_ex = h_cv_e  # [W/(m²*K)] Exterior convective heat transfer coeficient
        self.h_cv_int = (
            h_cv_i  # [W/(m²*K)] Interior convective heat transfer coeficient
        )

        self.E_gv_ex = E_gv_exterior  # [W/m²] Longwave irradiance incoming at external glass portion from the exterior
        self.E_gv_int = E_gv_interior  # [W/m²] Longwave irradiance incoming at internal glass portion from the interior


###############################################################################
# Definition of a class with all the outputs


class Outputs:
    def __init__(self):

        # Declaration and initialization of outputs:

        self.T_ft = (
            []
        )  # [K] The front temperatures of the layers, 0 =exterior, 1=layer 1 , 2= layer 2,..., n+1= last layer, n+2=interior
        self.T_b = (
            []
        )  # [K] The back temperatures of the layers, 0 =exterior, 1=layer 1 , 2= layer 2,..., n+1= last layer, n+2=interior
        self.J_ft = (
            []
        )  # [W/m²] The total longwave radiation going from the front of layer i towards layer i-1
        self.J_b = (
            []
        )  # [W/m²] The total longwave radiation going from the back of layer i towards layer i+1
        self.T_gap = (
            []
        )  # [K]  In case of ventilated cavities, the average air temperature in the cavities
        self.T_gap_out = (
            []
        )  # [K] In case of ventilated cavities, the air outlet temperature of the cavities.
        self.T_gap_middle = (
            []
        )  # [K] In case of ventilated cavities, the air temperature at middle height
        self.T_gap_inlet = []  # [K] Just to check
        self.h_cv_ft = (
            []
        )  # [W/(m²*K)] In case of closed cavity, it is the CHTC between the layer i and the layer i-1 or the exterior. In case of ventilated cavity, it is the CHTC between the layer i and the air between i and i-1.
        self.h_cv_b = (
            []
        )  # [W/(m²*K)] In case of closed cavity, it is the CHTC between the layer i and the layer i+11 or the interior. In case of ventilated cavity, it is the CHTC between the layer i and the air between i and i+1.
        self.q_int = (
            -77777
        )  # [W/m²] The total (convection+longwave radiation) heat flux going from the room towards the window
        self.U_gv = (
            -77777
        )  # [W/(m²*K)]. The U-value as defined by the standard. Only valid without irradiance.
        self.nat_airflow_rate = (
            -77777
        )  # [kg/s] In case of natural convection through holes between the inner and outer cavities

    def reset(self):
        self.T_ft = (
            []
        )  # [K] The front temperatures of the layers, 0 =exterior, 1=layer 1 , 2= layer 2,..., n+1= last layer, n+2=interior
        self.T_b = (
            []
        )  # [K] The back temperatures of the layers, 0 =exterior, 1=layer 1 , 2= layer 2,..., n+1= last layer, n+2=interior
        self.J_ft = (
            []
        )  # [W/m²] The total longwave radiation going from the front of layer i towards layer i-1
        self.J_b = (
            []
        )  # [W/m²] The total longwave radiation going from the back of layer i towards layer i+1
        self.T_gap = (
            []
        )  # [K]  In case of ventilated cavities, the average air temperature in the cavities
        self.T_gap_out = (
            []
        )  # [K] In case of ventilated cavities, the air outlet temperature of the cavities.
        self.T_gap_middle = (
            []
        )  # [K] In case of ventilated cavities, the air temperature at middle height
        self.T_gap_inlet = []  # [K] Just to check
        self.h_cv_ft = (
            []
        )  # [W/(m²*K)] In case of closed cavity, it is the CHTC between the layer i and the layer i-1 or the exterior. In case of ventilated cavity, it is the CHTC between the layer i and the air between i and i-1.
        self.h_cv_b = (
            []
        )  # [W/(m²*K)] In case of closed cavity, it is the CHTC between the layer i and the layer i+11 or the interior. In case of ventilated cavity, it is the CHTC between the layer i and the air between i and i+1.
        self.q_int = (
            -77777
        )  # [W/m²] The total (convection+longwave radiation) heat flux going from the room towards the window
        self.U_gv = (
            -77777
        )  # [W/(m²*K)]. The U-value as defined by the standard. Only valid without irradiance.
        self.nat_airflow_rate = (
            -77777
        )  # [kg/s] In case of natural convection through holes between the inner and outer cavities


###############################################################################
# Various functions:


def gas_mixture_properties(
    T_gas, air_ratio, argon_ratio, krypton_ratio, xenon_ratio, P=101325
):
    # Calculates the thermal properties of a gas mixture
    # From ISO15099 p. 18/19
    # Attention: T_gas in K
    # 0=Air, 1=Argon, 2=Krypton, 3=Xenon
    lambda_air, mu_air, Cp_air, rho_air, M_air = gas_properties(T_gas, 0)
    lambda_argon, mu_argon, Cp_argon, rho_argon, M_argon = gas_properties(T_gas, 1)
    lambda_krypton, mu_krypton, Cp_krypton, rho_krypton, M_krypton = gas_properties(
        T_gas, 2
    )
    lambda_xenon, mu_xenon, Cp_xenon, rho_xenon, M_xenon = gas_properties(T_gas, 3)

    # Molecular mass
    molecular_mass_mix = (
        (air_ratio * M_air)
        + (argon_ratio * M_argon)
        + (krypton_ratio * M_krypton)
        + (xenon_ratio * M_xenon)
    )

    # Density
    rho_mix = P * molecular_mass_mix / (R * T_gas * 1000)

    # Specifiec heat
    cpi_air = Cp_air * M_air
    cpi_argon = Cp_argon * M_argon
    cpi_krypton = Cp_krypton * M_krypton
    cpi_xenon = Cp_xenon * M_xenon
    cp_mix = (
        (air_ratio * cpi_air)
        + (argon_ratio * cpi_argon)
        + (krypton_ratio * cpi_krypton)
        + (xenon_ratio * cpi_xenon)
    )
    Cp_mix = cp_mix / molecular_mass_mix

    # Viscosity
    mu_mix = (
        (air_ratio * mu_air)
        + (argon_ratio * mu_argon)
        + (krypton_ratio * mu_krypton)
        + (xenon_ratio * mu_xenon)
    )

    # Thermal conductivity
    lambda_mix = (
        (air_ratio * lambda_air)
        + (argon_ratio * lambda_argon)
        + (krypton_ratio * lambda_krypton)
        + (xenon_ratio * lambda_xenon)
    )

    return lambda_mix, mu_mix, Cp_mix, rho_mix


def gas_properties(T_gas, gas=0, P=101325):
    # Calculates the thermal properties of different gases, based on Annex B of ISO15099
    # ThP 2015
    # Attention: T_gas in K
    # 0=Air, 1=Argon, 2=Krypton, 3=Xenon

    lambda_a = [2.873e-3, 2.285e-3, 9.443e-4, 4.538e-4]
    lambda_b = [7.760e-5, 5.149e-5, 2.826e-5, 1.723e-5]
    mu_a = [3.723e-6, 3.379e-6, 2.213e-6, 1.069e-6]
    mu_b = [4.94e-8, 6.451e-8, 7.777e-8, 7.414e-8]
    Cp_a = [1002.7370, 521.9285, 248.0907, 158.3397]
    Cp_b = [1.2324e-2, 0, 0, 0]
    molecular_masses = [
        28.97,
        39.948,
        83.80,
        131.30,
    ]  # [kg/kgmol] Molecular masses of the different gases

    lambda_gas = lambda_a[gas] + lambda_b[gas] * T_gas  # [W/(m*K)] Thermal conductivity
    mu_gas = mu_a[gas] + mu_b[gas] * T_gas  # [Pa*s] Dynamic viscosity
    Cp_gas = (
        Cp_a[gas] + Cp_b[gas] * T_gas
    )  # [J/(kg*K)]  Specific heat capacity at constant pressure
    rho_gas = (
        P * molecular_masses[gas] / (R * T_gas * 1000)
    )  # [kg/m³] Density of the gas
    M_gas = molecular_masses[gas]
    return lambda_gas, mu_gas, Cp_gas, rho_gas, M_gas


def h_cv_closed_cavity(T1, T2, i, param):
    # Calculates the convective heat transfer coefficient of vertical closed cavities
    # From ISO15099 p. 17
    # ThP 2015
    # d_gv [m] thickness of the gas layer
    # gas type: 0=Air, 1=Argon, 2=Krypton, 3=Xenon

    Tm = (T1 + T2) / 2  # Mean surface temperatures [K]

    # We get the gas properties:
    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        Tm,
        param.air_ratio[i],
        param.argon_ratio[i],
        param.krypton_ratio[i],
        param.xenon_ratio[i],
    )

    Ra = ((rho_gas ** 2) * ((param.d_gv[i]) ** 3) * g * Cp_gas * abs(T1 - T2)) / (
        mu_gas * lambda_gas * Tm
    )  # [-] Rayleigh number of the air in the cavity

    A_gv = param.height_glazed_area / param.d_gv[i]  # [-] Aspect ratio of the cavity

    Nu_1 = 0

    if 0 < Ra <= 1e4:
        Nu_1 = 1 + ((1.7596678e-10) * (Ra ** 2.2984755))
    elif 1e4 < Ra <= 5e4:
        Nu_1 = 0.028154 * (Ra ** 0.4134)
    elif Ra > 5e4:
        Nu_1 = 0.0673838 * (Ra ** (1 / 3))

    Nu_2 = 0.242 * ((Ra / A_gv) ** 0.272)

    Nu = max(Nu_1, Nu_2)

    h_cv = Nu * (
        lambda_gas / param.d_gv[i]
    )  # [W/(m²*K)] Convective heat transfer coefficient from one surface to the other

    return h_cv


def h_cv_vent_cavity(T1, T2, i, param, V):
    # h_cv in [W/(m²*K)]
    # This function calculates the convective heat transfer coefficient in ventilated gaps
    # It is the HTC between one surface and the average air temperature in the ventilated cavity

    h_cv_vent = 2 * h_cv_closed_cavity(T1, T2, i, param) + (4 * V)

    return h_cv_vent


def ventilated_case(
    T1, T2, i, T_gap_1_2, param, nat_mflow_r=-77777, nat_T_gap_inlet=-77777
):
    # The following calculations of values are needed in the calculation of q_i and q_i_plus_1
    # The outputs are used in the calcultaions of q_i, q_i_plus_1 and T_gap[i]

    T_av_1_2 = 0.5 * (T2 + T1)

    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        T_gap_1_2,
        param.air_ratio[i],
        param.argon_ratio[i],
        param.krypton_ratio[i],
        param.xenon_ratio[i],
    )
    if nat_mflow_r > 0:  # [kg/s] Natural ventilation, holes model
        V = nat_mflow_r / (
            rho_gas * param.d_gv[i] * param.width_glazed_area
        )  # [m/s] The mean velocity in the cavity
    else:
        V = param.airflow_rates[i] / (
            3600 * param.d_gv[i] * param.width_glazed_area
        )  # [m/s] The mean velocity in the cavity
    h_cv_1_2 = h_cv_vent_cavity(T1, T2, i, param, V)
    H_0_1_2 = (
        rho_gas * Cp_gas * param.d_gv[i] * V / (2 * h_cv_1_2)
    )  # [m] The characteristic height of the temperatur profile
    if nat_T_gap_inlet > 0:
        T_gap_inl_1_2 = nat_T_gap_inlet
    else:
        T_gap_inl_1_2 = param.T_air_inlet[
            i
        ]  # [K] We assume the entering temperature for all cavities is the exterior temperature TO BE CHANGED
    T_gap_out_1_2 = T_av_1_2 - (
        (T_av_1_2 - T_gap_inl_1_2) * exp(-param.height_glazed_area / H_0_1_2)
    )  # [K] The gap outlet temperature
    T_gap_middle = T_av_1_2 - (
        (T_av_1_2 - T_gap_inl_1_2) * exp(-param.height_glazed_area / (2 * H_0_1_2))
    )
    return h_cv_1_2, T_av_1_2, H_0_1_2, T_gap_inl_1_2, T_gap_out_1_2, T_gap_middle


def h_cv_internal(T_bn, T_ai_int, height_glazed_area):
    # Calculates the convective heat transfer coefficient for the internal side of the window
    # in case of natural convection which usually occurs on the internal side
    # From ISO15099 p.51/52

    # T_int [K] internal temperature
    # T_bn [K] temperature of the internal glazing surface
    # height_glazed_area [m] is the height of the fenestration system

    T_mf = T_ai_int + 0.25 * (T_bn - T_ai_int)  # [K] Mean film temperature
    # print("T_mf "+str(T_mf))
    # We get the gas properties for this mean film temperature :
    lambda_gas, mu_gas, Cp_gas, rho_gas, M_gas = gas_properties(T_mf, 0)

    Ra_H = (
        (rho_gas ** 2) * (height_glazed_area ** 3) * g * Cp_gas * abs(T_bn - T_ai_int)
    ) / (
        mu_gas * lambda_gas * T_mf
    )  # [-] Rayleigh number of the air inside
    # print("Ra_H "+str(Ra_H))

    Ra_cv = 2.5 * (10 ** 5) * ((exp(0.72 * 90) / 1) ** (1 / 5))
    # print("Ra_cv "+str(Ra_cv))

    # Nusselt number
    if Ra_H <= Ra_cv:
        Nu_int = 0.56 * ((Ra_H * 1) ** (1 / 4))
    else:  # Ra_H>Ra_cv :
        Nu_int = 0.13 * ((Ra_H ** (1 / 3)) - (Ra_cv ** (1 / 3))) + 0.56 * (
            (Ra_cv * 1) ** (1 / 4)
        )
    # print("Nu "+str(Nu_int))

    h_cv_int = Nu_int * (
        lambda_gas / height_glazed_area
    )  # [W/(m²*K)] Convective heat transfer coefficient for the internal side

    return h_cv_int


def h_cv_external_building(gamma_az, gamma_wN, V_wind):
    # Calculates the convective heat transfer coefficient for the external side of the window
    # in case of forced convection which usally occurs on the external side because of the wind
    # For real building fenestration component annual energy analysis
    # From ISO15099 p.53/54

    # gamma_az [°] Wall azimuth (positive degrees westward from south and negative eastward)
    # gamma_wN [°] Wind direction (angle measured clockwise from north)
    # V_wind [m/s] wind velocity measured at a height of 10m above ground level
    # Vs [m/s] free stream velocity near the fenestration surface

    # To determine whether the surface is windward or leeward :
    # We determine gamma_w [°] the wind direction relative to the wall surface

    gamma_w = gamma_az + 180 - gamma_wN  # [°] Wind direction

    if abs(gamma_w) > 180:
        gamma_w = 360 - (abs(gamma_w))

    if (-45 <= abs(gamma_w)) and (abs(gamma_w) <= 45):
        windward = 1
        leeward = 0
    else:
        windward = 0
        leeward = 1

    if windward == 1 and leeward == 0:  # If the surface is windward
        if V_wind > 2:
            Vs = 0.25 * V_wind
        elif V_wind <= 2:
            Vs = 0.5
    elif windward == 0 and leeward == 1:  # If the surface is leeward
        Vs = 0.3 + (0.05 * V_wind)

    h_cv_ext = 4.7 + (
        7.6 * Vs
    )  # [W/(m²*K)] Convective heat transfer coefficient for the external side

    return h_cv_ext


def h_cv_simple(param):
    h_cv_simple = 4 + (4 * param.V_wind)
    return h_cv_simple


def print_outputs(param, otpts):
    empty_strg = "\n\t"
    first_strg = "\n\t"
    mid_strg = "\n{0:5.1f}   ".format(otpts.T_b[0] - T_ref)
    last_strg = "\n\t"
    for i in range(0, param.N_layers):
        empty_strg = empty_strg + "\t|\t"
        mid_strg = mid_strg + "   {0:5.1f}|{1:5.1f}  ".format(
            (otpts.T_ft[i + 1] - T_ref), (otpts.T_b[i + 1] - T_ref)
        )
        if param.airflow_rates[i + 1] > 0:
            first_strg = first_strg + "{0:5.1f}   |\t".format(
                (otpts.T_gap_out[i + 1] - T_ref)
            )
            last_strg = last_strg + "{0:5.1f}   |\t".format(
                (param.T_air_inlet[i + 1] - T_ref)
            )
        else:
            first_strg = first_strg + "\t|\t"
            last_strg = last_strg + "\t|\t"
        # mid_strg=mid_strg+str(T_ft_final[i+1]-T_ref)+"|"+str(T_b_final[i+1]-T_ref)
    mid_strg = mid_strg + "   {0:5.1f}".format(otpts.T_ft[param.N_layers + 1] - T_ref)
    #        print("\n################################################################\n")
    #        print("Convergence achieved !")

    print(
        first_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + mid_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + empty_strg
        + last_strg
    )
    # print ("\nResults: \n\t"+str(result.tolist()))
    # print("\nEquations: \n\t"+str(equations(result,param)))

    print(
        "\nT_ai_ex: "
        + str(param.T_ai_ex - T_ref)
        + " / T_air_int: "
        + str(param.T_ai_int - T_ref)
        + " / E_gv_ex: "
        + str(param.E_gv_ex)
        + " / E_gv_int: "
        + str(param.E_gv_int)
    )
    print(
        "\nTemperature of the air entering in the cavity: "
        + str(param.T_air_inlet[2] - T_ref)
    )
    print("\nT_ft: " + str([T - T_ref for T in otpts.T_ft]))
    print("\nT_b: " + str([T - T_ref for T in otpts.T_b]))
    print("\nJ_ft: " + str(otpts.J_ft))
    print("\nJ_b: " + str(otpts.J_b))
    print("\nT_gap: " + str([T - T_ref for T in otpts.T_gap]))
    print("\nT_gap_inl: " + str([T - T_ref for T in otpts.T_gap_inlet]))
    print("\nT_gap_out: " + str([T - T_ref for T in otpts.T_gap_out]))
    print("\nT_gap_middle: " + str([T - T_ref for T in otpts.T_gap_middle]))
    print("\nh_cv_ft: " + str(otpts.h_cv_ft))
    print("\nh_cv_b: " + str(otpts.h_cv_b))
    print("\nq_int: " + str(otpts.q_int))
    print("\nU_gv: " + str(otpts.U_gv))
    print(
        "\nNatural airflow rate: "
        + str(otpts.nat_airflow_rate)
        + " (between outer and inner air gap in case of top and bottom holes)"
    )
    print("\n################################################################\n")


def dP_gap_i(mflow, T_gap_i, T_gap_k, param, i, A_holes):
    # This function returns the sum of the pressure drop (Pa) in gap i if air is circulating between gap i and gap k
    sum_delta_P_gap = 0
    gap_thickness_i = param.d_gv[i]  # [m] the thickness of the air gap
    gap_width_i = param.width_glazed_area  # [m] the width of the air gap
    A_s_i = gap_thickness_i * gap_width_i  # [m²] Cross section of the air gap
    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        T_gap_i,
        param.air_ratio[i],
        param.argon_ratio[i],
        param.krypton_ratio[i],
        param.xenon_ratio[i],
    )
    V = mflow / (rho_gas * A_s_i)  # [m/s] Mean air speed in the air gap
    delta_P_B_i = (
        rho_gas * (V ** 2) / 2
    )  # [Pa] Acceleration of the air to the velocity V
    delta_P_HP_i = (
        12 * mu_gas * param.height_glazed_area * V / (gap_thickness_i ** 2)
    )  # [Pa] Steady laminar flow (Haben-Poiseuille law)
    if T_gap_i > T_gap_k:
        A_equ_inl_i = A_holes
        A_equ_out_i = A_holes
    else:
        A_equ_inl_i = A_holes
        A_equ_out_i = A_holes
    Z_inl_i = (A_s_i / (0.6 * A_equ_inl_i) - 1) ** 2
    Z_out_i = (A_s_i / (0.6 * A_equ_out_i) - 1) ** 2
    delta_P_Z_i = (
        rho_gas * (V ** 2) * (Z_inl_i + Z_out_i) / 2
    )  # [Pa] Pressure loss in the inlet and outlet opening
    # print("i: "+str(i))
    # print("delta_P_B_i="+str(delta_P_B_i))
    # print("delta_P_HP_i="+str(delta_P_HP_i))
    # print("delta_P_Z_i="+str(delta_P_Z_i))
    sum_delta_P_gap = delta_P_B_i + delta_P_HP_i + delta_P_Z_i
    return sum_delta_P_gap


def dP_driving_i_k(T_gap_i, T_gap_k, param, i):
    # This function calculates the driving pressure difference for natural convection through holes between two cavities
    T_0 = 283  # [K] Reference temperature
    rho_0 = gas_mixture_properties(
        T_0,
        param.air_ratio[i],
        param.argon_ratio[i],
        param.krypton_ratio[i],
        param.xenon_ratio[i],
    )[
        -1
    ]  # [kg/m³] Density at the reference temperature
    dP_driving = (
        rho_0
        * T_0
        * g
        * param.height_glazed_area
        * abs(T_gap_i - T_gap_k)
        / (T_gap_i * T_gap_k)
    )
    # print("dP_driving="+str(dP_driving))
    return dP_driving


###############################################################################
# Set up of the equations


def equations(var, param2):
    # In case of non ventilated cavities, we have 4 unknown variable for each layer: T_ft,T_b,J_ft,J_b
    # In case of ventilated cavity, we add a variable T_gap and an equation per ventilated cavity

    param = param2[0]  # the parameters
    outpts = param2[1]  # the outputs passed through the calling of the function
    outpts.reset()  # to erase the old outputs by the new outputs at each iteration
    A_holes = param2[2]
    the_equations = []

    T_ft = []  # [K] Front temperatures of layers
    T_b = []  # [K] Back temperatures of layers
    J_ft = (
        []
    )  # [W/m²] The amount of longwave radiation flux going from the front of one layer to the back of the next layer in the outter direction
    J_b = (
        []
    )  # [W/m²] The amount of longwave radiation flux going from the back of one layer to the front of the next layer in the inner direction
    T_gap = (
        []
    )  # [K] The equivalent temperature of the air gaps in case of ventilated cavities
    T_gap_outlet = (
        []
    )  # [K] The gap outlet temperature in case of natural convection between the inner and outer gap

    # List to contain the outputs:
    T_gap_middle_tab = (
        []
    )  # [K] The temperature at a height H from the bottom of the cavity, here we will use the height H/2
    T_gap_inl_tab = (
        []
    )  # [K] Equivalent temperature of the air on the inlet of the gap in case of ventilated cavity
    T_gap_out_tab = (
        []
    )  # [K] Equivalent temperature of the air on the outlet of the gap in case of ventilated cavity
    h_cv_ft_tab = (
        []
    )  # [W/(m²*K)] Convective heat transfer coefficient on the front of the layer
    h_cv_b_tab = (
        []
    )  # [W/(m²*K)] Convective heat transfer coefficient on the back of the layer
    # h_cv_ft_tab[i] corresponds to the convective heat transfer coefficient between layers i and i-1
    # h_cv_b_tab[i] corresponds to the convective heat transfer coefficient between layers i and i+1

    # Transmission of boundaries and variables:

    # Layer 0 = exterior
    T_ft.append(-77777)  # Only for consistency of index
    T_b.append(param.T_ai_ex)
    J_ft.append(-77777)  # Only for consistency of index
    J_b.append(param.E_gv_ex)
    T_gap.append(-77777)  # Only useful for ventilated gaps
    T_gap_outlet.append(
        -77777
    )  # Only in case of natural convection between inner and outer gap
    # Extraction of variables for layers 1 to N:
    # The variable list has following structure: T_ft_1,_T_b_1,J_ft_1,J_b_1,...,T_ft_N,_T_b_N,J_ft_N,J_b_N

    shifter = 0

    for i2 in range(0, param.N_layers):

        T_ft.append(var[i2 * 4 + shifter])
        T_b.append(var[i2 * 4 + 1 + shifter])
        J_ft.append(var[i2 * 4 + 2 + shifter])
        J_b.append(var[i2 * 4 + 3 + shifter])
        if param.airflow_rates[i2 + 1] > 0 or (
            (i2 == 1 or i2 == (param.N_layers - 1)) and (A_holes > 0)
        ):
            shifter = shifter + 1
            T_gap.append(var[i2 * 4 + 3 + shifter])
            if (i2 == 1 or i2 == (param.N_layers - 1)) and (
                A_holes > 0
            ):  # We add the gap outlet temperatur in case of convection between the inner and outer gap
                shifter = shifter + 1
                T_gap_outlet.append(var[i2 * 4 + 3 + shifter])
        else:
            T_gap.append(-77777)
            T_gap_outlet.append(-77777)
    # Layer N+1=interior
    T_ft.append(param.T_ai_int)
    T_b.append(-77777)  # Not used, only for consistency of index
    J_ft.append(param.E_gv_int)
    J_b.append(-77777)  # Not used, only for consistency of index
    T_gap.append(-77777)  # Not used, only for consistency of index
    T_gap_outlet.append(-77777)  # Not used, only for consistency of index

    # In case of natural ventilation between the outer and inner gap, "holes model":
    if A_holes > 0:
        nat_mflow_rate = var[-1]  # [kg/s]
        # print("nat_mflow_rate: "+str(nat_mflow_rate))

    # We build the system of equations:

    for i in range(1, param.N_layers + 1):

        # Calculation of q_i, heat flux leaving the front of layer i:

        if param.airflow_rates[i] <= 0 and not (
            (i == 2 or i == param.N_layers) and (A_holes > 0)
        ):  # Non-ventilated case:
            if i == 1:
                h_cv_ft = param.h_cv_ex
                h_cv_ft_tab.extend(
                    [-77777, h_cv_ft]
                )  # First value of convective heat transfer coefficient
                # for the front side is the external coefficient
                h_cv_b_tab.append(h_cv_ft)
                T_gap_inl_tab.extend(
                    [-77777, -77777]
                )  # Not used, only for consistency of index
                T_gap_out_tab.extend(
                    [-77777, -77777]
                )  # Not used, only for consistency of index
                T_gap_middle_tab.extend([-77777, -77777])
            elif 1 < i:
                h_cv_ft = h_cv_closed_cavity(T_b[i - 1], T_ft[i], i, param)
                h_cv_ft_tab.append(h_cv_ft)

                T_gap_inl_tab.append(-77777)  # Not used, only for consistency of index
                T_gap_out_tab.append(-77777)  # Not used, only for consistency of index
                T_gap_middle_tab.append(-77777)

            q_i = (
                h_cv_ft * (T_ft[i] - T_b[i - 1]) + J_ft[i] - J_b[i - 1]
            )  # [W/m²] Heat fluxes from one layer to the other, positive from interior to exterior

        else:  # Ventilated case:
            if (i == 2 or i == param.N_layers) and (
                A_holes > 0
            ):  # Natural ventilation with holes
                if i == 2:
                    T_gap_in = T_gap_outlet[param.N_layers]
                if i == param.N_layers:
                    T_gap_in = T_gap_outlet[2]
                (
                    h_cv_ft,
                    T_av_i,
                    H_0_i,
                    T_gap_inl_i,
                    T_gap_out_i,
                    T_gap_middle_i,
                ) = ventilated_case(
                    T_b[i - 1], T_ft[i], i, T_gap[i], param, nat_mflow_rate, T_gap_in
                )
                h_cv_ft_tab.append(h_cv_ft)
                T_gap_inl_tab.append(T_gap_inl_i)
                T_gap_out_tab.append(T_gap_out_i)
                T_gap_middle_tab.append(T_gap_middle_i)
                q_i = h_cv_ft * (T_ft[i] - T_gap[i]) + J_ft[i] - J_b[i - 1]
            else:
                (
                    h_cv_ft,
                    T_av_i,
                    H_0_i,
                    T_gap_inl_i,
                    T_gap_out_i,
                    T_gap_middle_i,
                ) = ventilated_case(T_b[i - 1], T_ft[i], i, T_gap[i], param)
                h_cv_ft_tab.append(h_cv_ft)
                T_gap_inl_tab.append(T_gap_inl_i)
                T_gap_out_tab.append(T_gap_out_i)
                T_gap_middle_tab.append(T_gap_middle_i)
                q_i = h_cv_ft * (T_ft[i] - T_gap[i]) + J_ft[i] - J_b[i - 1]

        # Calculation of q_i_plus_1, heat flux entering the back of layer i:

        if param.airflow_rates[i + 1] <= 0 and not (
            (i == 1 or i == (param.N_layers - 1)) and (A_holes > 0)
        ):  # Non-ventilated case:
            if i == param.N_layers:
                # h_cv_b=h_cv_int
                h_cv_b = param.h_cv_int
                h_cv_b_tab.extend(
                    [h_cv_b, -77777]
                )  # Last value of convective heat transfer coefficient
                # for the back side is the internal coefficient
                h_cv_ft_tab.append(h_cv_b)
                T_gap_inl_tab.append(-77777)
                T_gap_out_tab.append(-77777)
                T_gap_middle_tab.append(-77777)
            elif i < param.N_layers:
                h_cv_b = h_cv_closed_cavity(T_b[i], T_ft[i + 1], i + 1, param)
                h_cv_b_tab.append(h_cv_b)

            q_i_plus_1 = (
                h_cv_b * (T_ft[i + 1] - T_b[i]) + J_ft[i + 1] - J_b[i]
            )  # [W/m²] Heat fluxes from one layer to the other, positive from interior to exterior

            if i == param.N_layers:
                outpts.q_int = q_i_plus_1

        else:  # Ventilated case:
            if (i == 1 or i == (param.N_layers - 1)) and (
                A_holes > 0
            ):  # Natural ventilation with holes
                if (
                    i == 1
                ):  # We have to calculate the outlet gap temperature of the inner gap which is the inlet temperature of this gap !
                    T_gap_in_plus_1 = T_gap_outlet[param.N_layers]
                if i == (
                    param.N_layers - 1
                ):  # We have to calculate the outlet gap temperature of the outer gap which is the inlet temperature of this gap !
                    T_gap_in_plus_1 = T_gap_outlet[2]
                (
                    h_cv_b,
                    T_av_i_plus_1,
                    H_0_i_plus_1,
                    T_gap_inl_i_plus_1,
                    T_gap_out_i_plus_1,
                    T_gap_middle_plus_1,
                ) = ventilated_case(
                    T_b[i],
                    T_ft[i + 1],
                    i + 1,
                    T_gap[i + 1],
                    param,
                    nat_mflow_rate,
                    T_gap_in_plus_1,
                )
            else:
                (
                    h_cv_b,
                    T_av_i_plus_1,
                    H_0_i_plus_1,
                    T_gap_inl_i_plus_1,
                    T_gap_out_i_plus_1,
                    T_gap_middle_plus_1,
                ) = ventilated_case(T_b[i], T_ft[i + 1], i + 1, T_gap[i + 1], param)

            h_cv_b_tab.append(h_cv_b)
            q_i_plus_1 = h_cv_b * (T_gap[i + 1] - T_b[i]) + J_ft[i + 1] - J_b[i]

            if i == param.N_layers:
                outpts.q_int = q_i_plus_1

        # Then, for each layer, we have 4 equations for 4 unknown variables T_ft, T_b, J_ft and J_b.
        the_equations.extend(
            [
                -q_i + param.S[i] + q_i_plus_1,
                -J_ft[i]
                + param.Epsilon_ft[i] * SB * ((T_ft[i]) ** 4)
                + param.Tau[i] * J_ft[i + 1]
                + param.r_ft[i] * J_b[i - 1],
                -J_b[i]
                + param.Epsilon_b[i] * SB * ((T_b[i]) ** 4)
                + param.Tau[i] * J_b[i - 1]
                + param.r_b[i] * J_ft[i + 1],
                -T_b[i]
                + T_ft[i]
                + param.t_gv[i]
                * (2 * q_i_plus_1 + param.S[i])
                / (2 * param.lambda_gv[i]),
            ]
        )
        # and one equation more if the gap is ventilated
        if param.airflow_rates[i] > 0 or (
            (i == 2 or i == (param.N_layers)) and (A_holes > 0)
        ):
            the_equations.extend(
                [
                    T_gap[i]
                    - (
                        T_av_i
                        - H_0_i * (T_gap_out_i - T_gap_inl_i) / param.height_glazed_area
                    )
                ]
            )
            # and one equation more i for the gap outlet temperatures if the "hole model" is activated, see ISO p.46-48
            if (i == 2 or i == (param.N_layers)) and (A_holes > 0):
                the_equations.extend([T_gap_outlet[i] - T_gap_out_i])
    # and one equation more i for the natural mass flow rate if the "hole model" is activated, see ISO p.46-48
    if A_holes > 0:
        the_equations.extend(
            [
                dP_driving_i_k(T_gap[2], T_gap[param.N_layers], param, 2)
                - dP_gap_i(
                    nat_mflow_rate, T_gap[2], T_gap[param.N_layers], param, 2, A_holes
                )
                - dP_gap_i(
                    nat_mflow_rate,
                    T_gap[param.N_layers],
                    T_gap[2],
                    param,
                    param.N_layers,
                    A_holes,
                )
            ]
        )
        # print(dP_driving_i_k(T_gap[2],T_gap[param.N_layers],param,2))
        # print(dP_gap_i(nat_mflow_rate,T_gap[2],T_gap[param.N_layers],param,2)+dP_gap_i(nat_mflow_rate,T_gap[param.N_layers],T_gap[2],param,param.N_layers))

    outpts.T_ft.extend(T_ft)
    outpts.T_b.extend(T_b)
    outpts.J_ft.extend(J_ft)
    outpts.J_b.extend(J_b)
    outpts.T_gap.extend(T_gap)
    outpts.T_gap_out.extend(T_gap_out_tab)
    outpts.T_gap_inlet.extend(T_gap_inl_tab)
    outpts.T_gap_middle.extend(T_gap_middle_tab)
    outpts.h_cv_ft.extend(h_cv_ft_tab)
    outpts.h_cv_b.extend(h_cv_b_tab)
    if A_holes > 0:
        outpts.nat_airflow_rate = nat_mflow_rate
    return the_equations


###############################################################################
# Solving of the equations


def equ_first_guess_holes(var, param2):
    nat_mflow_rate = var[0]
    param, T_outer_gap, T_inner_gap = param2[0], param2[1], param2[2]
    # print(dP_driving_i_k(T_outer_gap,T_inner_gap,param,2))
    return (
        dP_driving_i_k(T_outer_gap, T_inner_gap, param, 2)
        - dP_gap_i(nat_mflow_rate, T_outer_gap, T_inner_gap, param, 2)
        - dP_gap_i(nat_mflow_rate, T_inner_gap, T_outer_gap, param, param.N_layers)
    )


def calculate_ISO_15099(
    param, calc_outputs, A_holes=-77777, h_int=-77777, h_ext=-77777
):

    # A_holes [m²] the equivalent holes at the top and bottom, connecting the outer gap to the inner gap, according to IS 15099 p. 49. Index 2 is between cavity 1 and 2, index 2 between cavity 2 and 3, etc.
    # If A_holes=0.05 m², there is a holes of 0.05 m² at the top of the element connecting the inner cavity to the outer cavity and a 0.05 m² hole at the bottom allowing large scale convection

    # If we want to use h_int and h_ext instead of splitting convective and radiative part:
    if h_int > 0 and h_ext > 0:
        E_ft_1 = param.Epsilon_ft[1]
        param.Epsilon_ft[1] = 0
        r_ft_1 = param.r_ft[1]
        param.r_ft[1] = 0
        E_b_N = param.Epsilon_b[param.N_layers]
        param.Epsilon_b[param.N_layers] = 0
        r_b_N = param.r_b[param.N_layers]
        param.r_b[param.N_layers] = 0
        h_cv_e = param.h_cv_ex
        param.h_cv_ex = h_ext
        h_cv_i = param.h_cv_int
        param.h_cv_int = h_int
        E_gv_e = param.E_gv_ex
        param.E_gv_ex = 0
        E_gv_i = param.E_gv_int
        param.E_gv_int = 0
    #    print("\nCheck:\n"+str([param.Epsilon_ft,param.h_cv_ex,param.r_b]))

    # Starting values:
    starting_values = []

    # The variables' list has following structure: T_ft_1,_T_b_1,J_ft_1,J_b_1,...,,
    # T_ft_N,_T_b_N,J_ft_N,J_b_N
    for i1 in range(1, param.N_layers + 1):
        starting_values.append(
            param.T_ai_ex + (param.T_ai_int - param.T_ai_ex) / param.N_layers * i1 - 0.2
        )  # [°C] T_ft: Front temperatures of layers
        starting_values.append(
            param.T_ai_ex + (param.T_ai_int - param.T_ai_ex) / param.N_layers * i1 + 0.2
        )  # [°C] T_b: Back temperatures of layers
        starting_values.append(
            SB * (starting_values[-2]) ** 4
        )  # [W/m²] J_ft: The amount of longwave radiation flux going from the front of one layer to the back of the next layer in the outter direction
        starting_values.append(
            SB * (starting_values[-2]) ** 4
        )  # [W/m²] J_b: The amount of longwave radiation flux going from the back of one layer to the front of the next layer in the inner direction
        if param.airflow_rates[i1] > 0 and 2 <= i1 <= param.N_layers:
            starting_values.append(
                param.T_air_inlet[i1]
            )  # [K] In case of ventilated cavity, the average temperature of the cavity air
        if i1 == 2 and (A_holes > 0):
            starting_values.append(
                param.T_ai_ex
            )  # [K] The average temperature of the cavity air in case of natural convection between the inner air gap and the outer air gap, "holes" model
            starting_values.append(
                param.T_ai_ex
            )  # [K] The outlet temperature of the cavity air in case of natural convection between the inner air gap and the outer air gap, "holes" model
        if i1 == (param.N_layers) and (A_holes > 0):
            starting_values.append(
                param.T_ai_int
            )  # [K] The average temperature of the cavity air in case of natural convection between the inner air gap and the outer air gap, "holes" model
            starting_values.append(
                param.T_ai_int
            )  # [K] The outlet temperature of the cavity air in case of natural convection between the inner air gap and the outer air gap, "holes" model
    if A_holes > 0:
        # First guess of the massflow in case of natural convetion through holes between inner and outer gap
        m_flow_first_guess = fsolve(
            equ_first_guess_holes, 0.001, args=[param, param.T_ai_ex, param.T_ai_int]
        )[0]
        # print("Mass flow first guess: "+str(m_flow_first_guess)+" kg/s")
        starting_values.append(
            m_flow_first_guess
        )  # [kg/s] Natural massflow rate between the inner air gap and the outer air gap, "holes" model

    # The calculation
    result = fsolve(equations, starting_values, args=[param, calc_outputs, A_holes])

    if (
        max(equations(result, [param, calc_outputs, A_holes])) > 0.0005
    ):  # In case of NO convergence:
        print("\n################################################################\n")
        print("ALERT: NO CONVERGENCE")
        print("\n################################################################\n")
        calc_outputs.reset()
        return 0

    else:  # In case of convergence:
        #        print("Convergence achieved.")
        # If we used h_int and h_ext, we set back the parameters:
        if h_int > 0 and h_ext > 0:
            param.Epsilon_ft[1] = E_ft_1
            param.r_ft[1] = r_ft_1
            param.Epsilon_b[param.N_layers] = E_b_N
            param.r_b[param.N_layers] = r_b_N
            param.h_cv_ex = h_cv_e
            param.h_cv_int = h_cv_i
            param.E_gv_ex = E_gv_e
            param.E_gv_int = E_gv_i
        # Calculation of U-value if no irradiance:
        if max(param.S) <= 0:
            calc_outputs.U_gv = 1 / (
                (param.T_ai_int - param.T_ai_ex) / calc_outputs.q_int
            )  # [W*m-2*K-1] The U-value of the element
        return 1


###############################################################################
