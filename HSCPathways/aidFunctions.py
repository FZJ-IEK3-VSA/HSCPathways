# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 13:08:41 2016

@author: m.reuss
"""

import numpy as np
import CoolProp.CoolProp as CP
import pandas as pd

CP.set_config_string(
    CP.ALTERNATIVE_REFPROP_PATH,
    'C:\\Program Files (x86)\\REFPROP\\')
np.seterr(divide='ignore', invalid='ignore')

#%%H2 Constant Values at Normal Conditions
class H2Values (object):

    def __init__(self):
        self.M = 2.01588  # [kg/kmol], molare Masse von Wasserstoff
        # [J/kg K], spezifische Gaskonstante von Wasserstoff
        self.R_i = 4124.48269490247
        # [kg/m^3], Dichte von Wasserstoff im Normzustand
        self.roh_n = 0.089882
        # [-] Realgasfaktor von Wasserstoff im Normzustand
        self.Z_n = 1.00062387922965
        self.LHV_n = 119.833493175241  # [MJ/kg]
        self.g = 9.81  # [m/s^2], Erdbeschleunigung
        self.T_n = 273.15  # [K]
        self.p_n = 1.01325e5  # [Pa]

#%% Supporting Functions
def parabel(para, p):
    return (para[0] / 1e6) * (p / 1e5)**para[1] + \
        8.79676122460001e-06  # Parabelgleichung

def square(para, p):
    return para[0] * p**para[1] + para[2]  # squaregleichung
    

def getDiaPress(demArr, distArr, p_1, p_min):
    '''
    Calculation of Pipeline diameter and end pressure:
    Input Parameter:
    demArr=demand Array in kg/day
    distArr= distance Array in km
    p_1=Input Pressure at start of pipeline in bar
    p_min=minimal output pressure in bar
    '''
    #             Initialization              #
    V_para_parabel_20 = np.array([0.000125571318762396, 1.50162559878953])
    D_para_square_20 = np.array(
        [3.24859458677547e-06, 0.912591206027628, -0.166716162511868])
    Z_para_square_20 = np.array(
        [3.23101813258933e-09, 1.03880932425032, 1.00048097412768])
    T_m = np.array(20 + 273.15)  # K
    k = 0.02  # mm

    # Less diameter variances
    DK = np.linspace(0.1, 1.0, 901)  # Average class of diameter

    propH2 = H2Values()

    demHourly = demArr / 24 / 3600  # kg/day to kg/s

    distMeter = distArr * 1000  # km to m

    p_1 = p_1 * 1e5  # bar to Pa
    ###             Calculation                 ###
    res1 = len(distArr)
    res2 = demArr.shape[1]

    p_2 = np.zeros((res1, res2))
    w_1 = np.zeros((res1, res2))
    Re_1 = np.zeros((res1, res2))
    diameter = np.ones((res1, res2)) / 1000
    x = np.zeros((res1, res2))
    for i1 in range(demArr.shape[1]):
        for i2 in range(len(distArr)):
            while p_2[i2, i1] <= p_min * 1e5 or np.isnan(p_2[i2, i1]):
                # Calculation of Norm Volume Flow
                V_n = demHourly[0, i1] / propH2.roh_n  # m^3/s (i.N.)
                # Startwerte
                # Calculation of input density
                roh_1 = square(D_para_square_20, p_1[i2, i1])  # kg/m3
                # Volume flow at entrance
                V_1 = demHourly[0, i1] / roh_1  # m^3/s
                # inner diameter of the Pipeline
                diameter[i2, i1] = DK[x[i2, i1]]      # m
                # Velocity Entrance
                w_1[i2, i1] = V_1 / (np.pi * diameter[i2, i1]**2 / 4)
                # Berechnung der dynamischen Viskosität am Eintritt
                eta_1 = parabel(V_para_parabel_20, p_1[i2, i1])  # Pa*s
                # Berechnung der kinematischen Viskosität
                nu_1 = eta_1 / roh_1  # m^2/s
                # Berechnung der Reynoldszahl am Eintritt
                Re_1[i2, i1] = w_1[i2, i1] * diameter[i2, i1] / nu_1  # -
                # Berechnung der Rohrreibungszahl nach Zanke bei Re_1 für
                # Startwert
                alpha = np.e**(-1 * np.e**(6.75 - 0.0025 * Re_1[i2, i1]))
                lambda_1 = (64 / Re_1[i2, i1]) * (1 - alpha) + alpha * (-2 * np.log10((2.7 * (np.log10(
                    Re_1[i2, i1]))**1.2 / Re_1[i2, i1]) + (k / (3.71 * 1000 * diameter[i2, i1]))))**(-2)  # -
                # Simplification: Re_1 = Re_m --> lambda_m = lambda_1
                lambda_m = lambda_1
                # Berechnung der Leitungscharakteristik C_1
                # kg/(m s^2)=Pa
                C_1 = (lambda_1 * distMeter[i2] * roh_1 *
                       w_1[i2, i1]**2) / (diameter[i2, i1] * 2)
                # Ausgangsdruck bei raumbeständiger Fortleitung
                p_20 = p_1[i2, i1] - C_1  # Pa
                # Annahme: mittlere Druckes entspricht Ausgangsdruck bei
                # raumbeständiger Fortleitung
                p_m0 = p_20  # [Pa)
                # Annahme: mittlerer Realgasfaktor wird immer bei p_m0 bestimmt
                Z_m = square(Z_para_square_20, p_m0)
                # Berechnung der mittleren Kompressibilitätszahl
                K_m = Z_m / propH2.Z_n
                # Berechnung der Leitungscharakteristik C
                C = (lambda_m * 16 * propH2.roh_n * T_m * propH2.p_n *
                     K_m) / (np.pi**2 * propH2.T_n)  # kg Pa/m^3
                # Berechnung des Ausgangsdruckes
                p_2[i2, i1] = (p_1[i2, i1]**2 - (C * distMeter[i2]
                                                 * V_n**2) / diameter[i2, i1]**5)**0.5  # Pa

                if x[i2, i1] == len(DK):
                    break
                if p_2[i2, i1] <= p_min * 1e5 or np.isnan(p_2[i2, i1]):
                    x[i2, i1] += 1
                    x[i2:, i1:] = x[i2, i1]

    p_2 = p_2 * 1e-5
    diameter = diameter * 1000
    return diameter, p_2, w_1  # Diameter in mm and outlet pressure in bar

# %% Compressor Energy Demand per Stage (with isentropic coefficient)
# direct Method from Tietze


def getCompressionEnergyStage(p_1, p_2, T_1, eta_is_S):
    '''
    calculation of specific hydrogen compression energy in every compression stage
    Input:
    p_1=Inlet Pressure
    p_2=outlet Pressure
    T_1 = Inlet Temperature
    eta_is_S = isentropic efficiency
    '''
    fluid = 'HYDROGEN'
    # fluid='REFPROP::HYDROGEN'
    # Entropy
    s = CP.PropsSI('S', 'T', T_1, 'P', p_1 *
                   100000, fluid)  # [T]=K, [P]=kPa, [h]=J/kg
    # Enthalpy before 
    h_1 = CP.PropsSI('H', 'P', p_1 * 100000, 'S', s, fluid)
    # isentrope Enthalpie nach Verdichtung
    h_2_is = CP.PropsSI('H', 'P', p_2 * 100000, 'S', s,
                        fluid)  # [T]=K, [P]=kPa, [h]=J/kg
    # isentrope Endtemperatur
    # T_2_is = CP.PropsSI('T','P',p_2*100,'S',s,fluid); # [T]=K, [P]=kPa, [h]=J/kg
    # spez. isentrope Verdichterarbeit für reales Gas
    w_is = (h_2_is - h_1) / 1000  # [kJ/kg], massenspez. Verdichterarbeit
    # spezifische Verdichterarbeit für reales Gas
    w = w_is / eta_is_S  # [kJ/kg], massenspez. Verdichterarbeit
    w_spec = w / 3600
    # tatsächliche Enthalpie nach Verdichtung
    h_2 = w * 1000 + h_1  # [h]=J/kg
    # tatsächliche Temperatur nach Verdichtung
    T_2 = CP.PropsSI('T', 'P', p_2 * 100000, 'H', h_2,
                     fluid)  # [T]=K, [P]=kPa, [h]=J/kg
    return [w_spec, T_2]

# %% CompressionDemand


def getCompressionEnergy(
        p_1,
        p_2,
        demand,
        T_1=20,
        eta_isen=0.88,
        eta_mech=0.95,
        p_highlow_max=2.1,
        max_stages=2):
    '''
    calculation of specific hydrogen compression energy
    Input:
    p_1=Inlet Pressure in bar
    p_2=outlet Pressure in bar
    demand = hydrogen demand in kg/day
    T_1 = Inlet Temperature
    eta_is_S = isentropic efficiency
    '''
    # eta_isen=0.92-p_2/880*(0.24)
    if p_2 > p_1:

        compressorStages = np.log(p_2 / p_1) / np.log(p_highlow_max)
        compressorStages = np.ceil(compressorStages).astype(int)

        if compressorStages > max_stages:
            compressorStages = max_stages
        p_highlow = (p_2 / p_1)**(1 / compressorStages)
        # Initialize
        p_in = np.zeros(compressorStages)
        p_out = np.zeros(compressorStages)
        T_in = np.zeros(compressorStages)
        T_out = np.zeros(compressorStages)
        w_stage = np.zeros(compressorStages)
        # Stagedependent Calculation
        for i in range(compressorStages):
            if i == 0:
                p_in[i] = p_1
                T_in[i] = 273.15 + T_1
            else:
                p_in[i] = p_out[i - 1]
                T_in[i] = 273.15 + 40.
            p_out[i] = p_in[i] * p_highlow
            w_stage[i], T_out[i] = getCompressionEnergyStage(p_in[i],
                                                             p_out[i],
                                                             T_in[i],
                                                             eta_isen)
        T_out = T_out - 273.15
        w_mech = np.sum(w_stage) / eta_mech
        P_shaft = demand * w_mech / 24
        eta_motor = 8e-5 * np.log(P_shaft)**4 - 0.0015 * np.log(P_shaft)**3 + \
            0.0061 * np.log(P_shaft)**2 + 0.0311 * np.log(P_shaft) + 0.7617
        P_el = P_shaft / eta_motor
        w_total = w_mech / eta_motor
    else:
        w_total = 0
        P_el = 0

    return w_total, P_el
