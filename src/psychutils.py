# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 15:39:39 2014

@author: mstreet
"""
import numpy as np
import warnings


def cp_water(temp):
    """Function to determine the specific heat of water at
    constant pressure between 0 and 100 degrees Celsius.
    Interpolates a data set from TA Instruments that 
    references Osborne, Stimson and Ginnings 1939 in Handbook of
    Chemistry and Physics.

    param: temp: Water temperatures for the cp calculations.

    Calls the function numpy.interp()

    Returns the specific heat in [J/kgC].
    """

    t = np.array([0, 4, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

    cp = np.array(
        [
            4217.7,
            4204.8,
            4192.2,
            4181.9,
            4178.5,
            4178.6,
            4180.7,
            4184.4,
            4189.6,
            4196.4,
            4205.1,
            4216.0,
        ]
    )

    return np.interp(temp, t, cp)


def densityWater(temp):
    """ Funciton to determine the density of water between 0 and 
    100 degrees Celsius.  Interpolates a data set copied from
    Engineering Toolbox.

    param: temp: Water temperatures for the density calculations.

    Calls the function numpy.interp()
    """

    t = np.array([0, 4, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

    rho = np.array(
        [
            999.8,
            1000,
            999.7,
            998.2,
            995.7,
            992.2,
            988.1,
            983.2,
            977.8,
            971.8,
            965.3,
            958.4,
        ]
    )

    result = np.interp(temp, t, rho)

    return result


def spliceFunction(pos, neg, x, deltax):
    """ Function to provide a once continuously differentiable
    transition between two arguments.
    
    Equation from:
    Buildings.Utilities.Math.Functions.spliceFunction
    
    :param pos: Argument of x > 0.
    :param neg: Argument of x < 0.
    :param x:  Independent value.
    :param deltax: Half width of transition interval.
    """
    scaledx1 = x / deltax
    asin1 = np.arcsin(1)
    if scaledx1 <= -0.99999999999:
        out = neg
    elif scaledx1 >= 0.99999999999:
        out = pos
    else:
        y = (np.tanh(np.tan(scaledx1 * asin1)) + 1.0) / 2.0
        out = pos * y + (1.0 - y) * neg

    return out


def saturationPressureLiquid(Tsat):
    """ Return saturation pressure of water as a function of
    temperature T in the range of 273.16 to 373.16 K.
    
    Equation from:
    Buildings.Media.PerfectGases.MoistAirUnsaturated.saturationPressureLiquid
    
    
    :param Tsat: The temperature of liquid water in Kelvin
    """
    if Tsat < 223.16 or Tsat > 373.16:
        raise ValueError(
            "Tsat is "
            + str(Tsat)
            + " which is not in"
            + " the range 223.16 <= Tsat <= 373.16 K."
        )

    pSat = 611.657 * np.exp(17.2799 - 4102.99 / (Tsat - 35.719))

    return pSat


def sublimationPressureIce(Tsat):
    """Return sublimation pressure of water as a function of temperature
    between 223.16 and 273.16K
    
    Equation from:
    Modelica.Media.Air.MoistAirUnsaturated.sublimationPressureIce
    
    :param Tsat: The temperature of liquid water in Kelvin
    """
    if Tsat < 223.16 or Tsat > 373.16:
        raise ValueError(
            "Tsat is "
            + str(Tsat)
            + " which is not in"
            + " the range 223.16 <= Tsat <= 273.16 K for ice."
        )

    pSat = 611.657 * np.exp(22.5159 * (1.0 - 273.16 / Tsat))

    return pSat


def saturationPressure(Tsat):
    """ Return saturation pressure of moist air valid for temperatures
    223.16 <= Tsat <= 373.16 K.
    
    Equation from:
    Buildings.Media.PerfectGases.MoistAirUnsaturated.saturationPressure
    
    :param Tsat: saturation temperature [K]
    """

    if Tsat < 223.16 or Tsat > 373.16:
        raise ValueError(
            "Tsat is "
            + str(Tsat)
            + " which is not in"
            + " the range 223.16 <= Tsat <= 373.16 K."
        )

    pSat = spliceFunction(
        saturationPressureLiquid(Tsat), sublimationPressureIce(Tsat), Tsat - 273.16, 1.0
    )
    return pSat


def X_pSatpphi(pSat, p, phi):
    """ Return humidity ratio for given water vapor pressure.
    
    Equation from:
    Buildings.Utilities.Psychometrics.Functions.X_pSatpphi
    
    :param pSat: Saturation pressure [Pa]
    :param p: Fluid pressure [Pa]
    :param phi: Relative humidity [0 <= phi <= 1]
    """
    if phi > 1.01 and phi < 100.01:
        warnings.warn(
            "Assuming that the phi values are percentages.  " + "Dividing by 100.",
            RuntimeWarning,
        )
        phi = phi * 0.01
    elif phi < 0:
        raise ValueError("Relative humidity out of acceptable range.")

    k = 0.621964713077499  # Ratio of molar masses

    # Calculate water vapor concentration per total mass of air.
    X_w = phi * k / (k * phi + p / pSat - phi)

    return X_w


def X_pTphi(p_in, T, phi):
    """ Return steam mass fraction as a function of relative
    humidity, temperature and pressure.  Only for fluids with 
    a liquid water phase.
    
    Equation from:
    Buildings.Utilities.Psychometrics.X_pTphi
    
    :param p_in: Atmospheric pressure [Pa]
    :param T:  Temperature [K]
    :param phi: Relative humidity [0 <= phi <= 1]
    """
    if phi > 1.01 and phi < 100.01:
        warnings.warn(
            "Assuming that the phi values are percentages.  " + "Dividing by 100.",
            RuntimeWarning,
        )
        phi = phi * 0.01
    elif phi < 0:
        raise ValueError("Relative humidity out of acceptable range.")
    pSat = saturationPressure(T)
    # X = np.array([0,0])
    # X[0] = X_pSatpphi(pSat, p_in, phi)
    # X[1] = 1 - X[0]
    X = X_pSatpphi(pSat, p_in, phi)
    return X


def density_moist_air(T, phi, p=101325):

    # equations found on http://www.engineeringtoolbox.com/density-air-d_680.html

    # T air temperature in K
    # phi relative humidity between 0 and 1
    # p air pressure in Pa

    R_a = 286.9  # [J/(kg*K)]   ....of dry air
    R_w = 461.5  # [J/(kg*K)]   ....of saturated air

    density_da = p / (R_a * T)  # density of dry air

    x = X_pTphi(p, T, phi)  # [kg water/kg gas] humidity ratio

    density = density_da * (1 + x) / (1 + x * R_w / R_a)

    return density


def phi_pTX(p, T, X):
    """ Return Relative humidity [0 <= phi <= 1] as a function of water fraction
    , temperature and pressure.  Only for fluids with 

    a liquid water phase.
    
    Equation from:
    Buildings.Utilities.Psychometrics.X_pTphi
    
    :param p_in: Atmospheric pressure [Pa]

    :param T:  Temperature [K]
    :param phi: Relative humidity [0 <= phi <= 1]
    """

    pw = saturationPressure(T)

    Xs = X_pSatpphi(pw, p, 1)

    mu = X / Xs

    phi = mu / (1 - (1 - mu) * (pw / p))

    return phi
