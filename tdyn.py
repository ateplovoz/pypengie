# -*- coding: utf-8 -*-
"""
tdyn.py

Module with thermodynamics functions

---
Created: Vadim Pribylov, 2018-10-19
"""

import math

AIR_R = 287.06 # [J/(kg K)]

def cp_methane(t_start, t_end, afr, use_new=False):
    """
    Calculates specific heat of combustion products of natural gas.

    Args:
        t_start:float - lowest temperature of gas, °C
        t_end:float - highest temperature of gas, °C
        afr:float - air to fuel ratio. AFR=0 calculates specific heat of air,
            AFR in [0 ÷ 1]
        use_new:bool - use new table for polynom coefficients

    Returns:
        float - specific isobaric heat of gas mixture
    """
    # TODO: update style and give more clarity that this function can be used
    # to calculate specific heat of air.
    # TODO: change calculations to integrate specific heat instead of averaging
    # between Cp(t_start) and Cp(t_end).

    gas_coef = {
        'air': [.252192, -.593306, 11.2026, -76.8453, 276.44, -515.04, 392.2],
        'co2': [.104706, 2.11718, -13.1785, 56.2368, -154.60, 243.69, -166.7],
        'h2o': [.448938, -.544201, 13.4255, -65.9598, 159.88, -192.86, 89.17],
        'o2':  [.208363, -.056140, 7.45289, -68.3167, 292.27, -614.50,  512.0],
        'h':   [3.07088, 11.1537, -163.633, 1330.410, -5513.1, 11419., -9424.],
        }
    gas_coef_new = {
        'air': [1.004117, -2.754852*10**-6, 6.890151*10**-7, -8.358055*10**-10,
                3.343216*10**-13, 7.965546*10**-17, -1.175786*10**-19,
                3.856751*10**-23, -4.344547*10**-27],
        'co2': [0.8148627, 1.104033*10**-3, -1.301332*10**-6, 1.323860*10**-9,
                -1.118083*10**-12, 6.735382*10**-16, -2.561333*10**-19,
                5.420749*10**-23, -4.844547*10**-27],
        'h2o': [1.858979, 2.066147*10**-4, 1.409027*10**-6, -2.616702*10**-9,
                3.558973*10**-12, -3.276883*10**-15, 1.857165*10**-18,
                -6.186433*10**-22, 1.112322*10**-25, -8.334864*10**-30],
        'o2':  [0.9146970, 1.026171*10**-4, 1.142046*10**-6, -2.773659*10**-9,
                3.127974*10**-12, -1.990618*10**-15, 7.322529*10**-19,
                -1.451865*10**-22, 1.200398*10**-26],
        'h':   [14.19732, 3.926051*10**-3, -1.924455*10**-5, 4.817688*10**-8,
                -6.327742*10**-11, 5.072448*10**-14, -2.57018*10**-17,
                8.017741*10**-21, -1.402408*10**-24, 1.050463*10**-28],
            }
    g_coef = {
        'air': 0,
        'co2': 0.75,
        'h2o': 0.25,
        'o2': 0,
        'h': 0
        }

    poly_coef = []
    if use_new:
        t1 = t_start
        t2 = t_end
        ts = (t_start + t_end)/2
        gas_coef = gas_coef_new
    else:
        t1 = (t_start + 273.15)*0.0001
        t2 = (t_end + 273.15)*0.0001
        ts = (t1 + t2)/2.

    g_vozd = 1 - afr  # Количество воздуха в газах
    g_co2 = g_coef['co2']*afr*3.66666666  # Количество СО2 в газах
    g_h2o = g_coef['h2o']*afr*9  # Количество водяного пара в газах

    poly_coef = [
        gas_coef['air'][0]*g_vozd + gas_coef['co2'][0]*g_co2 +
        gas_coef['h2o'][0]*g_h2o,
        gas_coef['air'][1]*g_vozd + gas_coef['co2'][1]*g_co2 +
        gas_coef['h2o'][1]*g_h2o,
        gas_coef['air'][2]*g_vozd + gas_coef['co2'][2]*g_co2 +
        gas_coef['h2o'][2]*g_h2o,
        gas_coef['air'][3]*g_vozd + gas_coef['co2'][3]*g_co2 +
        gas_coef['h2o'][3]*g_h2o,
        gas_coef['air'][4]*g_vozd + gas_coef['co2'][4]*g_co2 +
        gas_coef['h2o'][4]*g_h2o,
        gas_coef['air'][5]*g_vozd + gas_coef['co2'][5]*g_co2 +
        gas_coef['h2o'][5]*g_h2o,
        gas_coef['air'][6]*g_vozd + gas_coef['co2'][6]*g_co2 +
        gas_coef['h2o'][6]*g_h2o
        ]

    if use_new:
        cp = 0
        i = 0
        for poly in poly_coef:
            cp += poly*pow(ts, i)*(i+1)
            i += 1
        return cp*1000
    else:
        cp = poly_coef[0] + poly_coef[1]*(t1 + t2)
        i = 2
        for poly in poly_coef[2:]:
            cp += poly*pow(ts, i)*(i+1)
            i += 1
        return cp*4186.8


def air_h(temperature, **kwargs):
    """
    Calculates specific enthalpy (enthalpy from 0 °C to `temperature`) of air
    based on temperature and pressure

    Source:
        Бурцев С. И. and Цветков, Ю. Н.
        Влажный воздух. Состав и свойства. Учебное пособие.
        Санкт-Петербург: СПбГАХПТ, 1998
        ISBN 5-89565-005-8

    Arguments
    ------
    temperature: float [K]
        temperature of air

    Kwargs
    ------
    use_new: bool
        if True then uses new tables to calculate specific isobaric heat of air

    Returns
    ------
    specific enthalpy [J/kg]
    """
    # TODO: generalize function so it can take any two of three arguments
    # (temperature, pressure, enthalpy) and calculate the unknown
    if 'use_new' in kwargs:
        air_cp = cp_methane(0, temperature - 273.15, 0, use_new=True)
    else:
        air_cp = cp_methane(0, temperature - 273.15, 0)
    return temperature * air_cp

def air_stp(temperature, pressure, **kwargs):
    """
    Calculates specific enthalpy (enthalpy from 0 °C to `temperature`) of air
    based on temperature and pressure

    Arguments
    ------
    temperature: float [K]
        temperature of air
    pressure: float [Pa]
        pressure of air

    Kwargs
    ------
    use_new: bool
        if True then uses new tables to calculate specific isobaric heat of air

    Returns
    ------
    specific entropy [J/(kg K)]
    """

    if 'use_new' in kwargs:
        air_cp = cp_methane(0, temperature, 0, use_new=True)
    else:
        air_cp = cp_methane(0, temperature, 0)
    return (
        air_cp*math.log(temperature/273.15) - AIR_R*math.log(pressure/101325)
    )


def air_tsp(entropy, pressure, init_temp, **kwargs):
    """
    Calculates temperature of air based on specific entropy and pressure

    Arguments
    ------
    entropy: float [J/(kg K)]
        specific entropy of air
    pressure: float [Pa]
        pressure of air
    init_temp: float [K]
        initial temperature, for specific heat calculation

    Kwargs
    ------
    use_new: bool
        if True then uses new tables to calculate specific isobaric heat of air
    use_kj: bool
        if True then uses kilojoules to work around floating point limit

    Returns
    ------
    temperature [K]
        temperature of the air
    """
    # TODO: write a paper about this

    if 'use_new' in kwargs:
        air_cp = cp_methane(0, init_temp, 0, use_new=True)
    else:
        air_cp = cp_methane(0, init_temp, 0)
    if 'use_kj' in kwargs:
        air_cp = air_cp/1000
        return (
            (
                math.e**(entropy/1000) * (pressure/101325)**(AIR_R/1000)
            )**(1/air_cp) * 273.15
        )
    else:
        return (
            (math.e**entropy * (pressure/101325)**AIR_R)**(1/air_cp) * 273.15
        )


def air_pst(entropy, temperature, **kwargs):
    """
    Calculates pressure of air based on specific entropy and temperature

    Arguments
    ------
    entropy: float [J/(kg K)]
        specific entropy of air
    temperature: float [K]
        temperature of air

    Kwargs
    ------
    use_new: bool
        if True then uses new tables to calculate specific isobaric heat of air

    Returns
    ------
    pressure [Pa]
        pressure of the air
    """
    # TODO: write a paper about this

    if 'use_new' in kwargs:
        air_cp = cp_methane(0, temperature, 0, use_new=True)
    else:
        air_cp = cp_methane(0, temperature, 0)
    return (
        (math.e**-entropy * (temperature/273.15)**air_cp)**(1/AIR_R) * 101325
    )
