'''
tcconvert.py

Module for thermocouple EMF and temperature conversion

---
Created: Vadim Pribylov, 2018-10-19
'''


def tc_emf(temp, tctype='k'):
    '''
    Calculates thermocouple EMF based on hot junction temperature
    (according to GOST 8.585–2001)

    Args:
        temp: temperature of thermocouple hot junction,
        tctype: type of the thermocouple (only `k` supported yet)

    Returns:
        emf: EMF at thermocouple cold junction

    ====================

    Thermocouple types explanation

    K-type:

    for temperature range from -270 °C to 0 °C
    $$
    \begin{align}
    E &= \sum\limits_{i=0}^{10} A_i t^i, \\
    A_0 &= 0, \\
    A_1 &= 3.9450128025 \cdot 10^{-2}, \\
    A_2 &= 2.622373598 \cdot 10^{-5}, \\
    A_3 &= -3.2858906784 \cdot 10^{-7}, \\
    A_4 &= -4.9904828777 \cdot 10^{-9}, \\
    A_5 &= -6.7509059173 \cdot 10^{-11}, \\
    A_6 &= -5.7410327428 \cdot 10^{-13}, \\
    A_7 &= -3.1088872894 \cdot 10^{-15}, \\
    A_8 &= -1.0451609365 \cdot 10^{-17}, \\
    A_9 &= -1.9889266878 \cdot 10^{-20}, \\
    A_{10} &= -1.6322697486 \cdot 10^{-23}; \\
    \end{align}
    $$

    for temperature range from 0 °C to 1372 °C:
    $$
    \begin{align}
    E &= \sum\limits_{i=0}^{9} A_i t^i + C_0 e^{C_1 (t - 126.9686)^2}
    A_0 &= -1.7600413686 \cdot 10^{-2},\\
    A_1 &=  3.8921204975 \cdot 10^{02},\\
    A_2 &=  1.8558770032 \cdot 10^{-5},\\
    A_3 &= -9.9457592874 \cdot 10^{-8},\\
    A_4 &=  3.1840945719 \cdot 10^{-10},\\
    A_5 &= -5.6072844889 \cdot 10^{-13},\\
    A_6 &=  5.6075059059 \cdot 10^{-16},\\
    A_7 &= -3.2020720003 \cdot 10^{-19},\\
    A_8 &=  9.7151147152 \cdot 10^{-23},\\
    A_9 &= -1.2104721275 \cdot 10^{-26},\\
    C_0 &=  1.185976 \cdot 10^{-1},\\
    C_1 &= -1.183432 \cdot 10^{-4}.
    \end{align}
    $$
    '''

    if tctype == 'k':
        if temp <= 0:
            coef = [
                3.9450128025 * 10 ** (-2),
                2.3622373598 * 10 ** (-5),
                -3.2858906784 * 10 ** (-7),
                -4.9904828777 * 10 ** (-9),
                -6.7509059173 * 10 ** (-11),
                -5.7410327428 * 10 ** (-13),
                -3.1088872894 * 10 ** (-15),
                -1.0451609365 * 10 ** (-17),
                -1.9889266878 * 10 ** (-20),
                -1.6322697486 * 10 ** (-23),
            ]
            const = None
        else:
            coef = [
                -1.7600413686 * 10 ** (-2),
                3.8921204975 * 10 ** (-2),
                1.8558770032 * 10 ** (-5),
                -9.9457592874 * 10 ** (-8),
                3.1840945719 * 10 ** (-10),
                -5.6072844889 * 10 ** (-13),
                5.6075059059 * 10 ** (-16),
                -3.2020720003 * 10 ** (-19),
                9.7151147152 * 10 ** (-23),
                -1.2104721275 * 10 ** (-26),
            ]
            const = [
                1.185976 * 10 ** (-1),
                -1.183432 * 10 ** (-4),
            ]
    else:
        raise AttributeError('thermocouple type not supported')

    emf = np.array(
        [ai * (temp ** order) for ai, order in zip(coef, range(len(coef)))]
    )
    if const:
        emf = (
            emf + const[0] * np.exp(const[1] * (temp - 126.9686) ** 2)
        )
    return emf.sum()


def tc_temp(emf, cj=0, tctype='k'):
    '''
    Calculates thermocouple hot junction temp based on emf
    (according to GOST 8.585–2001)

    for temperature range from -200 °C to 0 °C, TEMF from -5.891 to 0 mV:

    $$
    \begin{align}
    t &= \sum\limits_{i=0}^{8} C_i E^i,\\
    C_0 &= 0,\\
    C_1 &=  2.5173462 \cdot 10,\\
    C_2 &= -1.1662878,\\
    C_3 &= -1.0833638,\\
    C_4 &= -8.9773540 \cdot 10^{-1},\\
    C_5 &= -3.7342377 \cdot 10^{-1},\\
    C_6 &= -8.6632643 \cdot 10^{-2},\\
    C_7 &= -1.0450598 \cdot 10^{-2},\\
    C_8 &= -5.1920577 \cdot 10^{-4};\\
    \end{align}
    $$

    from 0 °C to 500 °C, TEMF from 0 to 20.644 mV:

    $$
    \begin{align}
    t &= \sum\limits_{i=0}^{9} C_i E^i,\\
    C_0 &= 0,\\
    C_1 &=  2.508355 \cdot 10,\\
    C_2 &=  7.860106 \cdot 10^{-2},\\
    C_3 &= -2.503131 \cdot 10^{-1},\\
    C_4 &=  8.315270 \cdot 10^{-2},\\
    C_5 &= -1.228034 \cdot 10^{-2},\\
    C_6 &=  9.804036 \cdot 10^{-4},\\
    C_7 &= -4.413030 \cdot 10^{-5},\\
    C_8 &=  1.057734 \cdot 10^{-6},\\
    C_9 &= -1.052755 \cdot 10^{-8};\\
    \end{align}
    $$

    from 500 °C to 1372 °C, TEMF from 20.644 to 54.886 mV:

    $$
    \begin{align}
    t &= \sum\limits_{i=0}^{6} C_i E^i,\\
    C_0 &= -1.318058 \cdot 10^{2},\\
    C_1 &=  4.830222 \cdot 10,\\
    C_2 &= -1.646031,\\
    C_3 &=  5.464731 \cdot 10^{-2},\\
    C_4 &= -9.650715 \cdot 10^{-4},\\
    C_5 &=  8.802193 \cdot 10^{-6},\\
    C_6 &= -3.110810 \cdot 10^{-8}.
    \end{align}
    $$
    '''

    if tctype == 'k':
        if emf < 0:
            coef = [
                0,
                2.5173462 * 10,
                -1.1662878,
                -1.0833638,
                -8.9773540 * 10 ** (-1),
                -3.7342377 * 10 ** (-1),
                -8.6632643 * 10 ** (-2),
                -1.0450598 * 10 ** (-2),
                -5.1920577 * 10 ** (-4),
            ]
        elif emf < 20.644:
            coef = [
                0,
                2.508355 * 10,
                7.860106 * 10 ** (-2),
                -2.503131 * 10 ** (-1),
                8.315270 * 10 ** (-2),
                -1.228034 * 10 ** (-2),
                9.804036 * 10 ** (-4),
                -4.413030 * 10 ** (-5),
                1.057734 * 10 ** (-6),
                -1.052755 * 10 ** (-8),
            ]
        else:
            coef = [
                -1.318058 * 10 ** (2),
                4.830222 * 10,
                -1.646031,
                5.464731 * 10 ** (-2),
                -9.650715 * 10 ** (-4),
                8.802193 * 10 ** (-6),
                -3.110810 * 10 ** (-8),
            ]
    else:
        raise AttributeError('thermocouple type not supported')

    temp = np.array(
        [ai * (emf ** order) for ai, order in zip(coef, range(len(coef)))]
    )
    return temp.sum() + cj
