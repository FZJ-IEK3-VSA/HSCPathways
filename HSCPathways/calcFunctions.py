# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 09:10:33 2016

@author: m.reuss
"""
from HSCPathways import classes as HSC
import numpy as np
import CoolProp.CoolProp as CP
import pandas as pd
from HSCPathways import aidFunctions as aFun

CP.set_config_string(
    CP.ALTERNATIVE_REFPROP_PATH,
    'C:\\Program Files (x86)\\REFPROP\\')
np.seterr(divide='ignore', invalid='ignore')


    #%% Calculation of Modules

def calcHSC(demArr, distArr, dfHSC, dfTable, i):

    a = HSC.Production(
        demArr,
        distArr,
        dfHSC['Production'][i],
        dfTable,
        i,
        dfHSC)

    # Connector1 Calculation
    b = HSC.Connector(demArr, distArr, dfTable, a.getTotalCost(), i, dfHSC)

    # Storage Calculation
    c = HSC.Storage(
        demArr,
        distArr,
        dfHSC['Storage'][i],
        dfTable,
        b.getTotalCost(),
        i,
        dfHSC)

    # Connector 2
    d = HSC.Connector2(demArr, distArr, dfTable, c.getTotalCost(), i, dfHSC)

    # Transport Calculation
    e = HSC.Transport(
        demArr,
        distArr,
        dfHSC['Transport'][i],
        dfTable,
        d.getTotalCost(),
        i,
        dfHSC)

    # Station Calculation
    f = HSC.Station(
        demArr,
        distArr,
        dfHSC['Station'][i],
        dfTable,
        e.getTotalCost(),
        i,
        dfHSC)

    return a, b, c, d, e, f


# %% Built up Hydrogen Supply Chains
def BuiltHSC(dfTable,
             conn=0,
             conversionStorageTransport=True,
             cavern=True,
             gTanks=False,             
             LOHCStationVariation=False            
             ):
    '''Production of different Hydrogen Supply Chains
    dfTable = DataFrame imported from excel
    conn==0 --> \n between storage and transport
    conn==1 --> space between storage and transport
    '''
    # Get the Names
    listProduction = list(dfTable['Production'])
    listTransport = list(dfTable['Transport'])
    listStorage = list(dfTable['Storage'])
    listStation = list(dfTable['Station'])
    listConnector = list(dfTable['Connector'])

    # Get the Form/State of the Chain Parts
    formProduction = dfTable['Production'].ix['form'].values
    formStorage = dfTable['Storage'].ix['form'].values
    formTransport = dfTable['Transport'].ix['form'].values
    systemConnector = dfTable['Connector'].ix['system'].values

    # Initialize the Pandas DataFrame for the Hydrogen Supply Chains
    dfHSC = pd.DataFrame()
    x = 0
    if conn == 0:
        co = '\n'
    else:
        co = ' + '
    # Name of the Columns
    col = ([['Production', 'Connector1', 'Storage', 'Connector2',
             'Transport', 'Station', 'General', 'StationTech']])
    arrHSC = np.array(col)

    # Do the loop for the Supply Chains: If the forms are compatible its fine
    for i1 in range(len(formProduction)):
        for i2 in range(len(formStorage)):
            for i3 in range(len(formTransport)):
                # Static Modules
                x1 = listProduction[i1]
                y1 = formProduction[i1]
                x3 = listStorage[i2]
                y3 = formStorage[i2]
                x5 = listTransport[i3]
                y5 = formTransport[i3]
                x6 = listStation[i3]
                # Connect Modules
                y2 = 10 * y1 + y3
                y4 = 10 * y3 + y5
                x7 = x3 + co + x5  # x1
                x8 = x5
                if i3 == 3 and LOHCStationVariation==True:
                    x7 = 'LOHC\n(nat. gas heated)'
                    x8 = 'LOHC\n(nat. gas heated)'
                # Kick out useless options
                if gTanks==False:
                    if x3 == 'GH2-Tank':  # and y5!=1:#Gtank as seasonal storage and LOHC/Liquid-Transport==boo
                        continue
                #Kick out Cavern+
                if cavern==False:
                    if x3 == 'GH2-Cavern':#Gtank as seasonal storage and LOHC/Liquid-Transport==boo
                        continue
                
                #LOHC comparison
                if LOHCStationVariation==True:
                    conversionStorageTransport==True
#                    if y5 != 3 and y3!=3:
#                        continue
                #kick out transport != storage
                if conversionStorageTransport==False:
                    if y5 != y3:
                        continue

                if y2 and y4 in systemConnector:
                    x2 = listConnector[np.where(systemConnector == y2)[0][0]]
                    x4 = listConnector[np.where(systemConnector == y4)[0][0]]
                    arrLoop = np.array([[x1, x2, x3, x4, x5, x6, x7, x8]])
                    arrHSC = np.vstack([arrHSC, arrLoop])
                    x = x + 1
                    if LOHCStationVariation==True:
                        if i3==3:
                            #LOHC Station Implementation (H2)
                            x6=listStation[4]
                            x7 = 'LOHC\n(hydrogen heated)' #x1
                            x8 = 'LOHC\n(hydrogen heated)'
                            x2 = listConnector[np.where(systemConnector==y2)[0][0]]
                            x4 = listConnector[np.where(systemConnector==y4)[0][0]]
                            arrLoop = np.array([[x1,x2,x3,x4,x5,x6,x7,x8]])
                            arrHSC=np.vstack([arrHSC,arrLoop])
                            x=x+1
                            # LOHC Station Electricity
                            x6 = listStation[5]
                            x7 = 'LOHC\n(electricity heated)'
                            x8 = 'LOHC\n(electricity heated)'
                            x2 = listConnector[np.where(systemConnector==y2)[0][0]]
                            x4 = listConnector[np.where(systemConnector==y4)[0][0]]
                            arrLoop = np.array([[x1,x2,x3,x4,x5,x6,x7,x8]])
                            arrHSC = np.vstack([arrHSC,arrLoop])
                            x = x + 1
                            #LOHC Station Diesel
                            x6 = listStation[6]
                            x7 = 'LOHC\n(diesel heated)' #x1
                            x8 = 'LOHC\n(diesel heated)'
                            x2 = listConnector[np.where(systemConnector==y2)[0][0]]
                            x4 = listConnector[np.where(systemConnector==y4)[0][0]]
                            arrLoop = np.array([[x1,x2,x3,x4,x5,x6,x7,x8]])
                            arrHSC = np.vstack([arrHSC,arrLoop])
                            x = x+1

                        
    dfHSC = pd.DataFrame(arrHSC[1:, :], columns=col)

    return dfHSC

# %% built arrays of Demand and Distance
def builtDemDistArray(demMax, distMax, res):
    distDamp = distMax / res
    distArr = np.linspace(distDamp, distMax, res)
    distArr = np.array([distArr]).T
    # Resolution Demand
    demDamp = demMax / res
    demArr = np.linspace(demDamp, demMax, res)
    demArr = np.array([demArr])
    return demArr, distArr