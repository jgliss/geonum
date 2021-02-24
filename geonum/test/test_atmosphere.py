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
import pytest
import numpy.testing as npt
from geonum import atmosphere as atm

def test_p0():
#: Pressure at sea level [Pa]
    assert atm.p0 == 101325.0

def test_NA():
    #: Avogadro constant [mol-1]
    assert atm.NA == 6.022140857e+23

def test_T0_STD():
    #: Standard sea level temperature [K]
    assert atm.T0_STD == 288.15

def test_M_AIR_AVG():
    #: Average molar mass of air, Unit: [g mol-1]
    assert atm.M_AIR_AVG == 28.9645

def test_R_STD():
    #: Gas constant (U.S. standard atmosphere) [Nm mol-1 K-1]
    assert atm.R_STD == 8.31432

def test_L_STD_ATM():
    #: Atmospheric lapse rate (Standard atmosphere) [K m-1]
    assert atm.L_STD_ATM == -6.5e-3

def test_L_DRY_AIR():
    #: Atmospheric lapse rate (dry atmosphere) [K m-1]
    assert atm.L_DRY_AIR == -9.8e-3

@pytest.mark.parametrize('phi,should_be',[
    (0, 9.780356),
    (0.79, 9.806398), #45deg
    (1.57, 9.83208), #90deg
    ])
def test__g0(phi,should_be):
    val = atm._g0(phi)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('phi,alt,should_be',[
    (0, 0, 9.780356),
    (0.79, 0, 9.806398), #45deg
    (1.57, 0, 9.83208), #90deg
    (0, 1000, 9.471702),
    (0.79, 1000, 9.497973), #45deg
    (1.57, 1000, 9.523879), #90deg
    ])
def test__g_acc(phi,alt,should_be):
    val = atm._g_acc(phi,alt)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('lat,alt,should_be',[
    (0, 0, 9.780356),
    (45, 0, 9.80616), #45deg
    (90, 0, 9.83208), #90deg
    (0, 1000, 9.471702),
    (45, 1000, 9.497733), #45deg
    (90, 1000, 9.523879), #90deg
    ])
def test_g(lat,alt,should_be):
    val = atm.g(lat,alt)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('lat,should_be',[
    (0, 9.780356),
    (45, 9.80616), #45deg
    (90, 9.83208), #90deg
    ])
def test_g0(lat,should_be):
    val = atm.g0(lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('alt,ref_temp,ref_alt,lapse_rate,should_be',[
    (0,atm.T0_STD,0,atm.L_STD_ATM,atm.T0_STD),
    (0,atm.T0_STD,0,atm.L_STD_ATM,288.15),
    (0,atm.T0_STD,0,atm.L_DRY_AIR,288.15),
    (1000,atm.T0_STD,0,atm.L_DRY_AIR,278.35),
    (1000,atm.T0_STD,0,atm.L_STD_ATM,281.65),
    (1000,atm.T0_STD,100,atm.L_STD_ATM,282.3),
    (1000,atm.T0_STD,1000,atm.L_STD_ATM,288.15),
    ])
def test_temperature(alt,ref_temp,ref_alt,lapse_rate,should_be):
    val = atm.temperature(alt,ref_temp,ref_alt,lapse_rate)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('mol_mass, lapse_rate, lat, should_be',[
    (atm.M_AIR_AVG, atm.L_STD_ATM, 45, -5.255632),
    (atm.M_AIR_AVG, atm.L_STD_ATM, 0, -5.241802),
    (atm.M_AIR_AVG, atm.L_STD_ATM, 90, -5.269523),
    (atm.M_AIR_AVG, atm.L_DRY_AIR, 45, -3.485878),
    ])
def test_beta_exp(mol_mass, lapse_rate, lat, should_be):
    val = atm.beta_exp(mol_mass, lapse_rate,lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be', [
    (0.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,101325),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,100129.493856),
    (1000.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,89875.07181),

    (0.0, 42, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,42),
    (0.0, 42, 273.15, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,42),
    (100.0, 42, 273.15, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,41.477378),

    (0.0, atm.p0, 273.15, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,101325),
    (100.0, atm.p0, 273.15, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,100064.17512),

    (0.0, atm.p0, 273.15, 100.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,102598.657973),

    (0.0, atm.p0, 273.15, 0,atm.L_DRY_AIR,atm.M_AIR_AVG,45.0,101325),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_DRY_AIR,atm.M_AIR_AVG,45.0,100128.81154),

    (0.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,20,45.0,101325),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,20,45.0,100497.98743),

    (0.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,0,101325),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,0,100132.621128),
    (0.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,90,101325),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,90,100126.352659),
    ])
def test_pressure(alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be):
    val = atm.pressure(alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)


@pytest.mark.parametrize('alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be', [
    (0.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,101325/100),
    (100.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,100129.493856/100),
    (1000.0, atm.p0, atm.T0_STD, 0.0,atm.L_STD_ATM,atm.M_AIR_AVG,45.0,89875.07181/100),
    ])
def test_pressure_hPa(alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be):
    val = atm.pressure_hPa(alt,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('p, ref_p, ref_temp, ref_alt, lapse_rate, mol_mass, lat, should_be', [
    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,0),
    (100000, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,110.889658),
    (50000, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,5574.679741),

    (atm.p0, atm.p0, 42,0,42,42,42,0),
    (50000, 50000, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,0),
    (50000, 80000, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,3792.338559),

    (atm.p0, atm.p0, 273.15,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,0),
    (50000, atm.p0, 273.15,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,5284.482982),
    (50000, atm.p0, 300,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,5803.935181),

    (atm.p0, atm.p0, atm.T0_STD,100,atm.L_STD_ATM,atm.M_AIR_AVG,45,100),

    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_DRY_AIR,atm.M_AIR_AVG,45,0),
    (50000, atm.p0, atm.T0_STD,0,atm.L_DRY_AIR,atm.M_AIR_AVG,45,5392.870582),

    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,20,45,0),
    (50000, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,20,45,7840.324557),

    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,0),
    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,0,0),
    (atm.p0, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,90,0),

    (50000, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,0,5588.419043),
    (50000, atm.p0, atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,90,5560.946563),

    ])
def test_pressure2altitude(p, ref_p, ref_temp, ref_alt, lapse_rate, mol_mass, lat, should_be):
    val = atm.pressure2altitude(p, ref_p, ref_temp, ref_alt, lapse_rate,
                                mol_mass, lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)

@pytest.mark.parametrize('alt,temp,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be', [

    (0,None,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1225.0033852145957),
    (100,None,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1213.286799),
    (1000,None,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1111.65185),
    (10000,None,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,412.733471),

    (0,273.15,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1292.274301),
    (0,300,atm.p0,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1176.615751),

    (0,None,50000,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,604.492171),
    (-3000,None,50000,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,798.755002),


    #(1000,None,50000,atm.T0_STD,0,atm.L_STD_ATM,atm.M_AIR_AVG,45,1225.0033852145957),
    ])
def test_density(alt,temp,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat,should_be):
    val = atm.density(alt,temp,ref_p,ref_temp,ref_alt,lapse_rate,mol_mass,lat)
    npt.assert_allclose(val, should_be, rtol=1e-7)

def test_number_density():
    pass

def test_refr_idx_300ppm_co2():
    pass

def test_refr_idx():
    pass

def test__F_O2():
    pass

def test_F_AIR():
    pass

def test_sigma_rayleigh():
    pass

def test_rayleigh_vol_sc_coeff():
    pass

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
