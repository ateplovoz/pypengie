'''
utils.py

Module for fluid dynamics and fluid parameters calculation

---
Created: Vadim Pribylov, 2018-10-19
'''

KAIR = 1.4
MAIR = 0.0405


def fluidflow_eqs(lang='en'):
    '''
    Simple function to get theoretical info about fluid flow dimensionless
    parameters.

    Basically copypasta for technical reports.

    Args:
        lang: str - language of output. Supports `en` and `ru` for now.

    Returns:
        str - handful of text with equations and parameter explanation
    '''
    eq = (
        '''
        $$\n
        \pi(\lambda) = \dfrac{p}{p_0} =
        \left(1-\dfrac{\kappa-1}{\kappa+1}\lambda^2\right)
        ^{\dfrac{\kappa}{\kappa-1}},\\\\\n
        q(\lambda) = \dfrac{f}{f_0} =
        \left(\dfrac{\kappa+1}{2}\right)^{\dfrac{1}{\kappa-1}}\lambda
        \left(1-\dfrac{\kappa-1}{\kappa+1}\lambda^2\right)
        ^{\dfrac{1}{\kappa-1}},\\\\\n
        G_\text{возд} =
        \dfrac{mq(\lambda_1)p_0F_0}{\sqrt{T_0}} \left[\dfrac{кг}{с}\right],
        '''
    )
    if lang == 'en':
        expl = (
            '''
            $p$ — static pressure at smallest diameter, kPa; $p_0$ — total
            pressure (atmospheric/ambient); $T_0$ — full air temperature, K;
            $\kappa=1.4$ — a constant; $\lambda$ — dimensionless flow speed;
            $\pi$ — dimensionless pressure, $q$ — dimensionless flow rate;
            $m=0.0405$ — a constant; $F_0$ — smallest area of nozzle
            '''
        )
    elif lang == 'ru':
        expl = (
            '''
            где $p$ — статическое давление в узком сечении мерного цилиндра,
            $p_0$ — полное давление (атмосферное), $T_0$ — полная температура
            воздуха, $\kappa=1{,}4$ — постоянная, $\lambda$ — безразмерная
            скорость, $\pi$ — безразмерное давление, $q$ — безразмерный расход,
            $m=0.0405$ — постоянная, $F_0$ — характерное сечение сопла.
            '''
        )
    return eq + '\n' + expl


def pi_lambda(g_lambda):
    '''
    Returns dimensionless pressure based on dimensionless speed
    '''
    return (
        pow(
            (1 - (KAIR - 1)/(KAIR + 1) * pow(g_lambda, 2)),
            KAIR/(KAIR - 1)
        )
    )


def lambda_pi(g_pi):
    '''
    Returns dimensionless flow speed based on dimensionless pressure
    '''
    return (
        pow(
            (1 - pow(g_pi, (KAIR - 1) / KAIR))
            * (KAIR + 1)/(KAIR - 1),
            0.5
            )
    )


def q_lambda(g_lambda):
    '''
    Returns dimensionless flow rate based on dimensionless speed
    '''
    return (
        pow(
            (KAIR + 1)/2,
            1 / (KAIR - 1)
        ) * g_lambda * pow(
            (1 - (KAIR - 1)/(KAIR + 1) * pow(g_lambda, 2)),
            1 / (KAIR - 1)
        )
    )


def g_q(g_q, p0, f0, t0):
    '''
    Returns flow rate based on dimensionless flow rate

    Args:
        g_q - input data
        p0 - ambient pressure
        f0 - smallest nozzle area
        t0 - ambient temperature
    '''
    return (
        (MAIR * g_q * p0 * f0) / pow(t0, 0.5)
    )
