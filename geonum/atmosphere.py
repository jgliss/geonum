# -*- coding: utf-8 -*-
#
# Geonum is a Python library for geographical calculations in 3D
# Copyright (C) 2017 Jonas Gliss (jonasgliss@gmail.com)
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License a
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Module for atmospheric calculations relevant for geonum
"""

from numpy import deg2rad, cos, ndarray, meshgrid, pi, exp, log

#: Pressure at sea level [Pa] 
p0 = 101325.0

#: Avogadro constant [mol-1]
NA = 6.022140857e+23

#: Standard sea level temperature [K]
T0_STD = 288.15 

#: Average molar mass of air, Unit: [g mol-1]
M_AIR_AVG = 28.9645 

#: Gas constant (U.S. standard atmosphere) [Nm mol-1 K-1]
R_STD = 8.31432  

#: Atmospheric lapse rate (Standard atmosphere) [K m-1]
L_STD_ATM = -6.5e-3
#: Atmospheric lapse rate (dry atmosphere) [K m-1]
L_DRY_AIR = -9.8e-3

def _g0(phi=0.0):
    """Gravitational accelaration at sea level (function of latitude)
    
    From Bodhaine et al., 1999, *On Rayleigh Optical Depth Calculations*

    
    Parameters
    ----------
    phi : float or ndarray
        latitude(s) in rads
    
    Returns
    -------
    float or ndarray
        value(s) of Earth gravitational accelaration at sea level for input 
        latitude(s)
    """
    return 9.806160 * (1 - 0.0026373 * cos(2 * phi) 
                         + 0.0000059 * cos(2 * phi)**2)

def _g_acc(phi, alt):
    return (_g0(phi) - (3.085e-4  + 2.27e-7 * cos(2 * phi)) * alt
                     + (7.254e-11 + 1e-13   * cos(2 * phi)) * alt ** 2
                     + (1.517e-17 + 6e-20   * cos(2 * phi)) * alt ** 3)
                         
def g(lat=45.0, alt=0.0):
    """Gravitational accelaration as function of latitude and altitude
    
    From Bodhaine et al., 1999, *On Rayleigh Optical Depth Calculations* who 
    have it from `List et al. 1968 <http://agris.fao.org/agris-search/search.
    do?recordID=US201300011727>`_.

    Parameters
    ----------
    lat : float or ndarray
        latitude(s) in degrees
    alt : float or ndarray
        altitude(s) in m
    
    Returns
    -------
    float or ndarray
        value(s) of Earth gravitational accelaration at sea level for input 
        latitude(s) and altitude(s). If both inputs are arrays, then a 2D array
        is returned with the first axis (index 0) corresponding to altitudes 
        and second axis to latitude
        
    Examples
    --------
    
    >>> import geonum
    >>> lats = np.linspace(0, 90, 200)
    >>> alts = np.linspace(0, 1000, 100)
    >>> dat = geonum.atmosphere.g(lats, alts)
    >>> dat.shape
    (100, 200)      
    """ 
    phi = deg2rad(lat)
    if all([isinstance(x, ndarray) for x in [phi, alt]]):
        phi, alt = meshgrid(phi, alt)
    
    return _g_acc(phi, alt)

def g0(lat=45.0):
    """Gravitational accelaration at sealevel as function of latitude
    
    See also :func:`g` for details
    
    Parameters
    ----------
    lat : float or ndarray
        latitude(s) in degrees
    alt : float or ndarray
        altitude(s) in m
    
    Returns
    -------
    float or ndarray
        value(s) of Earth gravitational accelaration at sea level for input 
        latitude(s) and altitude(s). If both inputs are arrays, then a 2D array
        is returned with the first axis (index 0) corresponding to altitudes 
        and second axis to latitude
    """
    return g(lat, alt=0.0)
                    
def temperature(alt=0.0, ref_temp=T0_STD, ref_alt=0.0, lapse_rate=L_STD_ATM):
    """Get temperature for a certain altitude (or altitude array)
    
    **Formula:**
        
    .. math::
        
        T = T_{ref} + L \\cdot (h-h_{ref}) \\quad [K]
        
    where:
        
        - :math:`$T_{ref}$` is a reference temperature
        - :math:`$L$` is the atmospheric lapse-rate
        - :math:`$h$` is the altitude in m
        - :math:`$h_{ref}$` is a reference altitude
        
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m 
    ref_temp : float, optional
        reference temperature in K (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to reference temperature (in m)
    lapse_rate : float, optional
        atmospheric lapse rate (in K m-1)
    
    Returns
    -------
    float or ndarray
        temperature(s) in K corresponding to altitudes 
    """    
    return float(ref_temp) + float(lapse_rate) * (alt - float(ref_alt))

def beta_exp(mol_mass=M_AIR_AVG, lapse_rate=L_STD_ATM, lat=45.0):
    """Exponent for conversion of atmosperic temperatures and pressures
    
    Based on barometric height formula (see e.g. `wikipedia <https://en.
    wikipedia.org/wiki/Barometric_formula>`__) for details.
    
    **Formula: **
        
    .. math:: 
        
        \\beta\,=\,\\frac{g_0 M_{air}}{R L}
        
    where:
        
        - :math:`$g_{0}$` is the gravitational constant at sea level
        - :math:`$M_{Air}$` is the molar mass of air (defaults to standard atm)
        - :math:`$R$` is the gas constant (:attr:`R_STD`)
        - :math:`$L$` is the atmospheric lapse-rate (cf. :attr:`L_STD_ATM`, :attr:`L_DRY_AIR`)
    
    Parameters
    ----------
    mol_mass : float, optional
        molar mass of air M 
    lapse_rate : float, optional
        atmospheric lapse rate L (in K m-1)
    lat : float, optional
        latitude for calculation of gravitational constant
        
    Returns
    -------
    float 
        Beta Exponent for formula 
    """
    return g0(lat) * mol_mass / (1000 * R_STD * lapse_rate)

def pressure(alt=0.0, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
             lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """Get atmospheric pressure in units of Pascal
    
    **Formula:**
    
    .. math:: 
        
        p = p_{ref} \\left[ \\frac{T_{ref}}{T_{ref} + L_{ref} (h-h_{ref})} \\right] ^\\beta
    
    where:
        
        - :math:`$\\beta$` is computed using :func:`beta_exp`
        - :math:`$p_{ref}$` is a reference pressure 
        - :math:`$T$` is the temperature  (cf. :func:`temperature`)
        - :math:`$T_{ref}$` is a reference temperature
        - :math:`$L_{ref}$` is the atmospheric lapse-rate 
        - :math:`$h$` is the altitude in m
        - :math:`$h_{ref}$` is a reference altitude
        
        
    
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m
    ref_p : float, optional
        reference pressure (default is std atm sea level)
    ref_temp : float, optional
        reference temperature in K corresponding to `ref_p`
        (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to `ref_p`
    lapse_rate : float, optional
        atmospheric lapse rate (in K / m)
    mol_mass : float, optional
        molar mass of air 
    lat : float, optional
        latitude for calculation of gravitational constant
        
    Returns
    -------
    float or ndarray
        pressure(s) corresponding to input altitude(s)
    """
    temp = temperature(alt, ref_temp, ref_alt, lapse_rate)
    exp = beta_exp(mol_mass, lapse_rate, lat)
    return ref_p * (ref_temp / temp) ** exp

def pressure_hPa(alt=0.0, *args, **kwargs):
    """Calls :func:`pressure` using provided input and returns in unit of hPa
    
    Defaults to standard atmosphere. This can be changed using further input
    parameters, for details see :func:`pressure`.
    
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m
    ref_p : float, optional
        reference pressure (default is std atm sea level)
    ref_temp : float, optional
        reference temperature in K corresponding to `ref_p`
        (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to `ref_p`
    lapse_rate : float, optional
        atmospheric lapse rate (in K / m)
    mol_mass : float, optional
        molar mass of air 
    lat : float, optional
        latitude for calculation of gravitational constant
    temp : float, optional
        if unspecified (default), the temperature is calculated using
        :func:`temperature`
        
    Returns
    -------
    float or ndarray
        pressure(s) in units of hPa corresponding to input altitude(s)
    """
    return pressure(alt, *args, **kwargs) / 100

def pressure2altitude(p, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
                      lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """General formula to convert atm. pressure to altitude
    
    **Formula:**
        
    .. math::
        
        h = h_{ref} + \\frac{T_{ref}}{L} \\left(\\exp\\left[-\\frac{\\ln\\left(\\frac{p}{p_{ref}} \\right)}{\\beta}\\right] - 1\\right) \\quad [m]
    
    where:
        
        - :math:`$h_{ref}$` is a reference altitude         
        - :math:`$T_{ref}$` is a reference temperature
        - :math:`$L$` is the atmospheric lapse-rate 
        - :math:`$p$` is the pressure (cf. :func:`pressure`)
        - :math:`$p_{ref}$` is a reference pressure
        - :math:`$\\beta$` is computed using :func:`beta_exp`
    
    Parameters
    ----------
    p : float
        pressure in Pa
    ref_p : float, optional
        reference pressure (default is std atm sea level)
    ref_temp : float, optional
        reference temperature in K corresponding to `ref_p`
        (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to `ref_p`
    lapse_rate : float, optional
        atmospheric lapse rate (in K / m)
    mol_mass : float, optional
        molar mass of air 
    lat : float, optional
        latitude for calculation of gravitational constant
        
    Returns
    -------
    float
        altitude corresponding to pressure level in defined atmosphere
    
        
    Examples
    --------
    >>> import geonum.atmosphere as atm
    >>> p = atm.pressure(2000) # pressure at 2000 m standard atm
    >>> int(atm.pressure2altitude(p))
    2000
    """
    
    beta = beta_exp(mol_mass, lapse_rate, lat=lat)
    return (ref_temp / lapse_rate * (exp(-log(p / ref_p) / beta) - 1) + ref_alt)

def density(alt=0.0, temp=None, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
            lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """Get atmospheric density in units of :math:`$g\,m^-3$`
    
    **Formula:**
        
    .. math::
        
        \\rho = p \\cdot \\frac{M_{Air}} {RT} \\quad \\left[ \\frac{g}{m^-3} \\right]
        
    where:
        
        - :math:`$p$` is the atm. pressure (cf. :func:`pressure`)
        - :math:`$M_{Air}$` is the molar mass of air (defaults to standard atm)
        - :math:`$R$` is the gas constant (:attr:`R_STD`)
        - :math:`$T$` is the temperature  (cf. :func:`temperature`)
    
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m
    temp : float, optional
        temperature in K, if None, :func:`temperature` is used to compute 
        temperature.
    ref_p : float, optional
        reference pressure (default is std atm sea level)
    ref_temp : float, optional
        reference temperature in K corresponding to `ref_p`
        (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to `ref_p`
    lapse_rate : float, optional
        atmospheric lapse rate (in K / m)
    mol_mass : float, optional
        molar mass of air 
    lat : float, optional
        latitude for calculation of gravitational constant
    
    Returns
    -------
    float or ndarray
        density corresponding to altitudes (in g m-3)
    """
    if temp is None:
        temp = temperature(alt, ref_temp, ref_alt, lapse_rate)
    p = pressure(alt, mol_mass=mol_mass, lapse_rate=lapse_rate, lat=lat)
    return p * mol_mass / (R_STD * temp)
    
def number_density(alt=0.0, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
                   lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """Get atmospheric number density in m-3
    
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m
    ref_p : float, optional
        reference pressure (default is std atm sea level)
    ref_temp : float, optional
        reference temperature in K corresponding to `ref_p`
        (default is std atm sea level)
    ref_alt : float, optional
        altitude corresponding to `ref_p`
    lapse_rate : float, optional
        atmospheric lapse rate (in K m-1)
    mol_mass : float, optional
        molar mass of air 
    lat : float, optional
        latitude for calculation of gravitational constant
        
    Returns
    -------
    float or ndarray
        number density in m-3 corresponding to altitudes
    """
    rho = density(alt, ref_p, ref_temp, ref_alt, lapse_rate, mol_mass, lat)
    return rho * NA / mol_mass

def refr_idx_300ppm_co2(lbda_mu=0.300):
    """Get refractive index of air using eq. of Peck and Reeder, 1972
    
    Uses either of two parametrisations for wavelengths above or below 0.23 
    microns using parametrisations from `Peck and Reeder, 1972, Dispersion of 
    air <https://www.osapublishing.org/josa/abstract.cfm?uri=josa-62-8-958>`_.
    
    Parameters 
    ----------
    lbda_mu : float or ndarray
        wavelength in microns
    
    Returns
    -------
    float
        refractive index
    """
    if lbda_mu >= 0.23:
        return 1 + (5791817 / (238.0185 - (1 / lbda_mu) ** 2) +
                    167909  / (57.362   - (1 / lbda_mu) ** 2)
                   ) * 10 ** (-8)
    else:
        return 1 + (8060.51 +
                    2480990 / (132.274  - (1 / lbda_mu) ** 2) +
                    17455.7 / (39.32957 - (1 / lbda_mu) ** 2)
                   ) * 10 ** (-8)
     
def refr_idx(lbda_mu=0.300, co2_ppm=400.0):
    """Calculate refr. index of atm for varying CO2 amount
    
    Uses either of two parametrisations for wavelengths above or below 0.23 
    microns from `Peck and Reeder, 1972, Dispersion of 
    air <https://www.osapublishing.org/josa/abstract.cfm?uri=josa-62-8-958>`_ 
    to retrieve the refractive index for CO2 concentrations at 300ppm and 
    the parametrisation for varying CO2 amounts from Bodhaine et al., 1999, 
    *On Rayleigh Optical Depth Calculations*.
    
    Parameters
    ----------
    lbda_mu : float or ndarray
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    float or ndarray
        refractive index
    """
    n_300 = refr_idx_300ppm_co2(lbda_mu)
    return 1 + (n_300 - 1) * (1 + 0.54 * (co2_ppm * 10**(-6) - 0.0003))
    
def _F_N2(lbda_mu=0.300):
    """Depolarisation factor of N2 for Rayleigh scattering cross section
    
    Taken from `Bates 1984 <http://www.sciencedirect.com/science/article/pii/
    0032063384901028>`_
    
    Parameters
    ----------
    lbda_mu : float or ndarray
        wavelength in micrometers
    
    Returns
    -------
    float or ndarray
        depolarisation factor of N2  
    """
    return 1.034 + 3.17e-4 / lbda_mu ** 2

def _F_O2(lbda_mu=0.300):
    """Depolarisation factor of O2 for Rayleigh scattering cross section
    
    Taken from `Bates 1984 <http://www.sciencedirect.com/science/article/pii/
    0032063384901028>`_
    
    Parameters
    ----------
    lbda_mu : float or ndarray
        wavelength in micrometers
    
    Returns
    -------
    float or ndarray
        depolarisation factor of O2
        
    """
    return 1.096 + 1.385e-3 / lbda_mu ** 2 + 1.448e-4 / lbda_mu**4
    
def F_AIR(lbda_mu=0.300, co2_ppm=400.0):
    """Depolarisation factor of air for Rayleigh scattering cross section
    
    Taken from Bodhaine et al., 1999, *On Rayleigh Optical Depth Calculations*
    
    Parameters
    ----------
    lbda_mu : float or ndarray
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    float or ndarray
        depolarisation factor of O2  
    """
    co2 = co2_ppm * 10**(-6)
    return (78.084 * _F_N2(lbda_mu) + 
            20.946 * _F_O2(lbda_mu) + 
            0.934 + co2 * 1.15) / (78.084 + 20.946 + 0.934 + co2)
            
def sigma_rayleigh(lbda_mu=0.300, co2_ppm=400.0):
    """Rayleigh scattering cross section
    
    Parameters
    ----------
    lbda_mu : float
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    float
        value of Rayleigh scattering cross section in cm-2    
    """
    n = refr_idx(lbda_mu, co2_ppm)
    lbda = lbda_mu * 10**(-6)
    num_dens = number_density()
    refr_fac = (n**2 - 1)**2 / (n**2 + 2)**2
    F = F_AIR(lbda_mu, co2_ppm)
    return 24 * pi**3 / (lbda**4 * num_dens**2) * refr_fac * F * 100 ** 2
           
    
def rayleigh_vol_sc_coeff(alt=0.0, lbda_mu=0.300, co2_ppm=400.0, **kwargs):
    """Rayleigh volume scattering coefficient :math:`\\beta(z,\\lambda)`
    
    Parameters
    ----------
    alt : float or ndarray
        altitude(s) in m
    lbda_mu : float
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    float or ndarray
        vol. scattering coeff. in cm-1 corresponding to altitudes
    """
    num_dens = number_density(alt, **kwargs) * 100**(-3) # cm^-3
    return num_dens * sigma_rayleigh(lbda_mu, co2_ppm)

if __name__ == '__main__':
    print('No input temp')
    print(density(0))
    
    T1 = 273
    print('input temp', T1)
    print(density(0, temp=T1))
    
    T2 = 253
    print('input temp', T2)
    print(density(0, temp=T2))