import numpy as np

def layer_temp(
        t_drop:float,
        t_gas:float,
        alpha:float=1/3
    ) -> float:
    """
    Calculate the layer temperature according to the droplet temperature and the gas temperature

    :param t_drop: Droplet temperature
    :param t_gas: Gas temperature
    :param alpha: Alpha parameter
    :return: Layer temperature

    """
    return t_drop + alpha*(t_gas - t_drop)


def sat_pressure_g(
        t_gas,
        a=77.3450,
        b=0.0057,
        c=7235,
        d=8.2
    ) -> float:
    """
    Calculate the gas saturation pressure according to the gas temperature

    :param t_gas: Gas temperature
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :param d: coefficient d
    :return: Saturation pressure

    """
    return (np.exp(a+b*t_gas-c/t_gas))/(t_gas**d)


def vap_pressure_g(
        ppm:float,
        press:float,
        mm_gas:float,
        mm_drop:float
    ) -> float:
    """
    Evaluate the gas vapour pressure according to the gas pressure, the gas molecular mass, the droplet molecular mass and the water content

    :param ppm: Water content
    :param press: Gas pressure
    :param mm_gas: Gas molecular mass
    :param mm_drop: Droplet molecular mass
    :return: Vapour pressure

    """
    return (ppm*mm_gas*press) / (mm_drop*1.0e6 + mm_gas*ppm)


def eval_omega(
        kbtepsilon12:float,
        a:float=0.49924,
        b:float=574.81152,
        c=1,
        d:float=11267.71941,
        e:float=1,
        f:float=1.45781
    ) -> float:
    """
    Evaluate the omega parameter according to the Boltzmann constant, the epsilon parameter and the temperature

    :param kbtepsilon12: Boltzmann constant * temperature / epsilon
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :param d: coefficient d
    :param e: coefficient e
    :param f: coefficient f
    :return: Omega parameter

    """
    return a+b/(c+kbtepsilon12*d)**(e/f)


def eval_diff_coeff(
        t_layer:float,
        omega:float,
        sigma12:float,
        p:float,
        mm_g:float,
        mm_d:float,
        a:float=1.86e-3,
        b:float=3/2,
        c:float=1/2,
        d:float=9.86923e-6,
        e:float=2
    ) -> float:
    """
    Evaluate the diffusion coefficient according to the layer temperature, the omega parameter, the surface tension, the gas pressure, the gas molecular mass and the droplet molecular mass

    :param t_layer: Layer temperature
    :param omega: Omega parameter
    :param sigma12: Surface tension
    :param p: Gas pressure
    :param mm_g: Gas molecular mass
    :param mm_d: Droplet molecular mass
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :param d: coefficient d
    :param e: coefficient e
    :return: Diffusion coefficient

    """
    return a*(t_layer**(b))*(1/mm_d+1/mm_g)**(c)/(p*d*(sigma12**e)*omega)


def eval_vap_heat(
        t_d:float,
        a:float=-1.66282e-7,
        b:float=2.53533e-4,
        c:float=-0.14657,
        d:float=35.23021,
        e:float=-434.96806
    ) ->float:
    """
    Evaluate the vapour heat according to the droplet temperature

    :param t_d: Droplet temperature
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :param d: coefficient d 
    :param e: coefficient e
    :return: Vapour heat

    """
    return a*(t_d)**4+b*(t_d)**3+c*(t_d)**2+d*(t_d)+e


def eval_viscosity(
        t_layer:float,
        a:float,
        b:float
    ) -> float:
    """
    Evaluate the viscosity according to the layer temperature

    :param t_layer: Layer temperature
    :param a: coefficient a
    :param b: coefficient b
    :return: Viscosity

    """
    return t_layer*a+b


def eval_term(
        mu_g_g:float,
        mu_l_g:float,
        chi_vs:float,
        mm_g:float,
        mm_d:float
    ) -> float:
    """
    Evaluate the term according to the gas viscosity, the liquid viscosity, the chi parameter, the gas molecular mass and the droplet molecular mass

    :param mu_g_g: Gas viscosity
    :param mu_l_g: Liquid viscosity
    :param chi_vs: Chi parameter
    :param mm_g: Gas molecular mass
    :param mm_d: Droplet molecular mass
    :return: Term
    
    """
    return mu_l_g/(1+(1-chi_vs)/chi_vs*((1+((mu_l_g/mu_g_g)**0.5)*((mm_d/mm_g)**0.25))**2)/(2.83*(1+mm_d/mm_g)**0.5)**0.5)


def eval_m_viscosity(
        mu_g_g:float,
        mu_l_g:float,
        chi_vs:float,
        mm_g:float,
        mm_d:float
    ) -> float:
    """
    Evaluate the viscosity according to the gas viscosity, the liquid viscosity, the chi parameter, the gas molecular mass and the droplet molecular mass

    :param mu_g_g: Gas viscosity
    :param mu_l_g: Liquid viscosity
    :param chi_vs: Chi parameter
    :param mm_g: Gas molecular mass
    :param mm_d: Droplet molecular mass
    :return: Viscosity

    """
    term_a = eval_term(mu_g_g, mu_l_g, chi_vs, mm_g, mm_d)
    term_b = eval_term(mu_l_g, mu_g_g, 1-chi_vs, mm_d, mm_g)
    return term_a+term_b


def eval_cp(
        t_m:float,
        a:float,
        b:float,
        c:float,
        d:float,
        e:float,
        mm_x:float
    ) -> float:
    """
    Evaluate the heat capacity according to the layer temperature, the gas molecular mass and the droplet molecular mass

    :param t_m: Layer temperature
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :param d: coefficient d
    :param e: coefficient e
    :param mm_x: Molecular mass
    :return: Heat capacity

    """
    return (a+b*(t_m)/1000+c*((t_m)/1000)**2+d*(t_m/1000)**3+e/((t_m/1000)**2))/mm_x


def eval_k(
        t_m:float,
        a:float,
        b:float,
        c:float
    ) -> float:
    """
    Evaluate the thermal conductivity according to the layer temperature

    :param t_m: Layer temperature
    :param a: coefficient a
    :param b: coefficient b
    :param c: coefficient c
    :return: Thermal conductivity

    """
    return a*t_m**2+b*t_m+c


def eval_phi(
        mu_g_g:float,
        mu_l_g:float,
        mm_g:float,
        mm_d:float
    ) -> float:
    """
    Evaluate the phi parameter according to the gas viscosity, the liquid viscosity, the gas molecular mass and the droplet molecular mass

    :param mu_g_g: Gas viscosity
    :param mu_l_g: Liquid viscosity
    :param mm_g: Gas molecular mass
    :param mm_d: Droplet molecular mass
    :return: Phi parameter
    
    """
    return (1+(mu_g_g/mu_l_g)**0.5*(mm_d/mm_g)**0.25)**2/(np.sqrt(8)*(1+(mm_g/mm_d))**0.5)
