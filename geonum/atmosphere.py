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

from numpy import deg2rad, cos, ndarray, meshgrid, pi
#from .constants 
#: Pressure at sea level, Unit: Pa (N/m^2)
p0 = 101325.0

#: Avogadro constant (g / mol)
NA = 6.022140857e+23

#: Standard sea level temperature, Unit: K
T0_STD = 288.15 

#: Average molar mass of air, Unit: g/mol
M_AIR_AVG = 28.9645 

#: Gas constant (U.S. standard atmosphere), Unit: N m mol^-1 K^-1
R_STD = 8.31432  

#: Atmospheric lapse rates, Unit: K / m
L_STD_ATM = -6.5e-3
L_DRY_AIR = -9.8e-3

def _g0(phi=0.0):
    """Gravitational accelaration at sea level (function of latitude)
    
    From Bodhaine et al., 1999, *On Rayleigh Optical Depth Calculations*

    
    Parameters
    ----------
    phi : :obj:`float`, :obj:`array`
        latitude(s) in rads
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
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
    lat : :obj:`float`, :obj:`array`
        latitude(s) in degrees
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
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


                    
def temperature(alt=0.0, ref_temp=T0_STD, ref_alt=0.0, lapse_rate=L_STD_ATM):
    """Get temperature for a certain altitude (or altitude array)
    
    Parameters
    ----------
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    ref_temp : :obj:`float`, optional
        reference temperature in K (default is std atm sea level)
    ref_alt : :obj:`float`, optional
        altitude corresponding to reference temperature (in m)
    lapse_rate : :obj:`float`, optional
        atmospheric lapse rate (in K / m)
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        temperature(s) in K corresponding to altitudes 
    """    
    return float(ref_temp) + float(lapse_rate) * (alt - float(ref_alt))
    
def pressure(alt=0.0, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
             lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0, temp=None):
    """Get atmospheric pressure in units of Pascal
    
    Parameters
    ----------
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    ref_p : :obj:`float`, optional
        reference pressure (default is std atm sea level)
    ref_temp : :obj:`float`, optional
        reference temperature in K corresponding to :param:`ref_pressure`
        (default is std atm sea level)
    ref_alt : :obj:`float`, optional
        altitude corresponding to :param:`ref_pressure`
    lapse_rate : :obj:`float`, optional
        atmospheric lapse rate (in K / m)
    mol_mass : :obj:`float`, optional
        molar mass of air 
    lat : :obj:`float`, optional
        latitude for calculation of gravitational constant
    temp : :obj:`float`, optional
        if unspecified (default), the temperature is calculated using
        :func:`temperature`
        
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        pressure(s) corresponding to input altitude(s)
    """
    if temp is None:
        temp = temperature(alt, ref_temp, ref_alt, lapse_rate)
    exp = (g(lat, alt) * mol_mass / (1000 * R_STD * lapse_rate))
    return ref_p * (ref_temp / temp) ** exp
    
def density(alt=0.0, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
            lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """Get atmospheric density in units of :math:`$g\,m^-3$`
    
    
    Parameters
    ----------
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    ref_p : :obj:`float`, optional
        reference pressure (default is std atm sea level)
    ref_temp : :obj:`float`, optional
        reference temperature in K corresponding to :param:`ref_pressure`
        (default is std atm sea level)
    ref_alt : :obj:`float`, optional
        altitude corresponding to :param:`ref_pressure`
    lapse_rate : :obj:`float`, optional
        atmospheric lapse rate (in K / m)
    mol_mass : :obj:`float`, optional
        molar mass of air 
    lat : :obj:`float`, optional
        latitude for calculation of gravitational constant
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        density corresponding to altitudes (in :math:`$g\,m^-3$`)
    """
    temp = temperature(alt, ref_temp, ref_alt, lapse_rate)
    p = pressure(alt, mol_mass=mol_mass, lapse_rate=lapse_rate, temp=temp,
                 lat=lat)
    return p * mol_mass / (R_STD * temp)
    
def number_density(alt=0.0, ref_p=p0, ref_temp=T0_STD, ref_alt=0.0, 
                   lapse_rate=L_STD_ATM, mol_mass=M_AIR_AVG, lat=45.0):
    """Get atmospheric number density :math:`$[m^-3]$`
    
    Parameters
    ----------
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    ref_p : :obj:`float`, optional
        reference pressure (default is std atm sea level)
    ref_temp : :obj:`float`, optional
        reference temperature in K corresponding to :param:`ref_pressure`
        (default is std atm sea level)
    ref_alt : :obj:`float`, optional
        altitude corresponding to :param:`ref_pressure`
    lapse_rate : :obj:`float`, optional
        atmospheric lapse rate (in K / m)
    mol_mass : :obj:`float`, optional
        molar mass of air 
    lat : :obj:`float`, optional
        latitude for calculation of gravitational constant
        
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        number density in m^-3 corresponding to altitudes
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
    lbda_mu : (:obj:`float`, :obj:`array`)
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
    lbda_mu : (:obj:`float`, :obj:`array`)
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    float
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
    lbda_mu : (:obj:`float`, :obj:`array`)
        wavelength in micrometers
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        depolarisation factor of N2
        
    """
    return 1.034 + 3.17e-4 / lbda_mu ** 2

def _F_O2(lbda_mu=0.300):
    """Depolarisation factor of O2 for Rayleigh scattering cross section
    
    Taken from `Bates 1984 <http://www.sciencedirect.com/science/article/pii/
    0032063384901028>`_
    
    Parameters
    ----------
    lbda_mu : (:obj:`float`, :obj:`array`)
        wavelength in micrometers
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        depolarisation factor of O2
        
    """
    return 1.096 + 1.385e-3 / lbda_mu ** 2 + 1.448e-4 / lbda_mu**4
    
def F_AIR(lbda_mu=0.300, co2_ppm=400.0):
    """Depolarisation factor of air for Rayleigh scattering cross section
    
    Taken from Bodhaine et al., 1999, *On Rayleigh Optical Depth Calculations*
    
    Parameters
    ----------
    lbda_mu : (:obj:`float`, :obj:`array`)
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
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
    alt : :obj:`float`, :obj:`array`
        altitude(s) in m
    lbda_mu : float
        wavelength in micrometers
    co2_ppm : float
        atmospheric CO2 volume mixing ratio in ppm
    
    Returns
    -------
    (:obj:`float` or :obj:`array`)
        vol. scattering coeff. in cm-1 corresponding to altitudes
    """
    num_dens = number_density(alt, **kwargs) * 100**(-3) # cm^-3
    return num_dens * sigma_rayleigh(lbda_mu, co2_ppm)