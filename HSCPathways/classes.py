# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:54:01 2016

@author: m.reuss
"""
import numpy as np
import pandas as pd
from HSCPathways import aidFunctions as aFun
import CoolProp.CoolProp as CP
#%% Motherclass
class Module(object):
    '''
    Generic mother class of all defined modules. 
    Inputs:
    form: Storage form like gaseous, liquid, lohc, metalhydride...
    pressure: in bar
    electricityCost: €/kWh prize for bought electricity
    NGCost: €/kWh prize for bought heat
    WACC: interestrate for the whole system
    demand: daily demand for hydrogen (normally vector)
    distance: distance between storage and distribution
    dfTable: Dataframe for Parameter
    '''
    
    def __init__(self, demand, distance,dfTable, i, dfHSC):
        #Initialize the Class with General Parameters
        self.GeneralTab=dfTable['General']
        #self.electricityCost = self.GeneralTab['General']['electricityCostAvg']
        self.electricityCostRES = self.GeneralTab['General']['electricityCostRES']
        self.NGCost = self.GeneralTab['General']['NGCost']
        self.heatGain = self.GeneralTab['General']['heatGain']
        self.WACC = self.GeneralTab['General']['WACC']
        self.storageDays=self.GeneralTab['General']['storageDays']
        self.driverCost=self.GeneralTab['General']['driverCost']
        self.demand=demand
        self.distance=distance
        self.waterCost=self.GeneralTab['General']['waterCost']
        self.dieselCost=self.GeneralTab['General']['dieselCost']
        self.hours=self.GeneralTab['General']['op. Hours RES']
        self.overDimension=1/(self.hours/8760)
        self.utilization=self.GeneralTab['General']['utilization Station']    
        self.electricityDemandPrecooling=self.GeneralTab['General']['electricityDemandPrecooling']   
        self.storagePart=self.GeneralTab['General']['storagePart']  
        self.distributionDistance=self.GeneralTab['General']['distributionDistance']
        
        self.residentialTime=self.storageDays/self.storagePart
        self.speed=self.GeneralTab['General']['truckSpeed']
                #Names of whole systems
        self.nameProduction=dfHSC['Production'][i]
        self.nameStorage=dfHSC['Storage'][i]
        self.nameStation=dfHSC['Station'][i]
        self.nameC1=dfHSC['Connector1'][i]
        self.nameC2=dfHSC['Connector2'][i]
        self.nameTransport=dfHSC['Transport'][i]
        
        if dfTable['Transport'][self.nameTransport]['form'] == dfTable['Production'][self.nameProduction]['form']:
            self.storagePart=self.storagePart
        else:
            self.storagePart=1.
        
    def getAnnuity(self,WACC,lifetime):
        '''
        Function to return the Annuity
        '''
        annuity=(WACC*(1+WACC)**lifetime)/((1+WACC)**lifetime-1)
        return annuity
        
    def getInvestScale(self,demand, base, compare, scale):
        '''
        Function for calculating the investcost based on a scaling function:
        y=base*(demand/compare)^scale
        '''
        invest=base*(self,demand/compare)**scale
        return invest
        
    def dayToSec(self,demand):
        '''
        mutation from kg/day to kg/s
        '''
        massflow=demand/24/3600
        return massflow

    def getElectricityCost(self,annualDemand):
        '''
        get electricity Price based on annual production
        values from eurostat 2016 http://appsso.eurostat.ec.europa.eu/nui/show.do
        '''
        kind=1          #0='nothing'
                        #1='noTaxes'
                        #2='all'
        
        eCost=annualDemand/1000
        if kind==0:
            eCost[eCost>70000.]=0.0555
            eCost[eCost>20000.]=0.0617
            eCost[eCost>2000.]=0.0696
            eCost[eCost>500.]=0.0811
            eCost[eCost>20.]=0.1020
            eCost[eCost>1.]=0.1342
        
        elif kind==1: 
            eCost[eCost>70000.]=0.0976
            eCost[eCost>20000.]=0.1118
            eCost[eCost>2000.]=0.1308
            eCost[eCost>500.]=0.1501
            eCost[eCost>20.]=0.1738
            eCost[eCost>1.]=0.219
        
        
        
        
        elif kind==2:
            eCost[eCost>70000.]=0.1344
            eCost[eCost>20000.]=0.1513
            eCost[eCost>2000.]=0.1740
            eCost[eCost>500.]=0.1970
            eCost[eCost>20.]=0.2251
            eCost[eCost>1.]=0.2789
        
        return eCost


#%%Production Instance
class Production(Module):
    def __init__(self, demand, distance, name, dfTable, i, dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):
        
        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Transport Technology
        '''
         #Initialize of Motherclass
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        self.electricityCost=self.electricityCostRES
        #Import Transport Table
        self.ProductionTab=dfTable['Production']
        
        self.form =self.ProductionTab[name]['form']
        self.pressureOut =self.ProductionTab[name]['pressureOut']
        self.investBase =self.ProductionTab[name]['investBase']
        self.investCompare =self.ProductionTab[name]['investCompare']
        self.investScale =self.ProductionTab[name]['investScale']
        self.investLifetime =self.ProductionTab[name]['investLifetime']
        self.boilOff =self.ProductionTab[name]['boilOff']        
        self.investOM =self.ProductionTab[name]['investOM']
        self.electricityDemand =self.ProductionTab[name]['electricityDemand']
        self.waterDemand =self.ProductionTab[name]['waterDemand']
        
        #
    def getTOTEX(self):
        
        self.annuity=self.getAnnuity(self.WACC, self.investLifetime)
        self.CAPEX=(self.annuity*self.investBase*(self.demand*self.overDimension/self.investCompare)**self.investScale)/(self.demand*365)
        self.fixOPEX=(self.investOM*self.investBase*(self.demand*self.overDimension/self.investCompare)**self.investScale)/(self.demand*365)
        self.varOPEX=self.electricityDemand*self.electricityCost+self.waterCost*self.waterDemand
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        self.electricityDemandAnnual=self.electricityDemand*self.demand*365
        
        return self.TOTEX
        
    def getTotalCost(self):
        return self.getTOTEX()
        
    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen, ElectricityRES, electricityGrid, Heat, Diesel
        return 0., self.electricityDemand,0., 0., 0.
    
    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX

        

#%% Storage Instance
class Storage(Module):
    def __init__(self, demand, distance, name,dfTable, costH2In, i, dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):
        
        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Storage Technology
        '''
        
        #Initialize of Motherclass
            #Getting the General Informations
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        
        #Import Storage Table
        self.StorageTab=dfTable['Storage']
        #Initialize of Parameter
        self.pressureOut =self.StorageTab[name]['pressureOut']
        self.pressureIn =self.StorageTab[name]['pressureIn']
        self.form = self.StorageTab[name]['form']
        self.investBase= self.StorageTab[name]['investBase']
        self.investCompare= self.StorageTab[name]['investCompare']
        self.investScale= self.StorageTab[name]['investScale']
        self.investLifetime= self.StorageTab[name]['investLifetime']
        self.boilOff= self.StorageTab[name]['boilOff']
        self.investOM= self.StorageTab[name]['investOM']
        self.costH2In=costH2In
        
        
    def getTOTEX(self):
        self.annuity=self.getAnnuity(self.WACC,self.investLifetime)
        self.capacity=self.storageDays*self.demand
        if self.nameStorage=='GH2-Cavern' or self.nameStorage=="GH2-Kaverne":
            self.density = CP.PropsSI('D','T',298,'P',self.pressureIn*1e5,'hydrogen')-CP.PropsSI('D','T',298,'P',self.pressureOut*1e5,'hydrogen'); # [T]=K, [P]=kPa, [h]=J/kg
            self.invest=self.investBase*(self.capacity/self.density/self.investCompare)**self.investScale
        else:
            self.invest=self.investBase*(self.capacity/self.investCompare)**self.investScale
        self.CAPEX=self.invest*self.annuity/(self.demand*365)
        self.fixOPEX=self.invest*self.investOM/(self.demand*365)
        self.varOPEX=self.storagePart*self.boilOff*self.costH2In*self.residentialTime
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        
        return self.TOTEX
        
    def getTotalCost(self):
        
        self.costH2Out=self.costH2In+self.getTOTEX()
        return self.costH2Out        
        
    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen Electricity, Heat, Diesel
        return 1+self.storagePart*self.boilOff*self.storageDays, 0., 0., 0., 0.
        
    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX        
        

#%%Storage Instance        
class Transport(Module):
    def __init__(self, demand, distance,name,dfTable, costH2In, i,dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):
        
        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Transport Technology
        '''
         #Initialize of Motherclass
            #Getting the General Informations
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        
        #Import Transport Table
        self.TransportTab=dfTable['Transport']       
        self.form =self.TransportTab[name]['form']
        self.pressureIn =np.ones((len(distance),demand.shape[1]))*self.TransportTab[name]['pressureIn']
        self.pipeInvestA=self.TransportTab[name]['pipeInvestA']
        self.pipeInvestB=self.TransportTab[name]['pipeInvestB']
        self.pipeInvestC=self.TransportTab[name]['pipeInvestC']
        self.pipeLifetime=self.TransportTab[name]['pipeLifetime']
        self.pipeHours=self.TransportTab[name]['pipeHours']
        self.pipeOM=self.TransportTab[name]['pipeOM']
        #self.pipeElectricityDemand=self.TransportTab[name]['pipeElectricityDemand']
        self.truckInvest=self.TransportTab[name]['truckInvest']
        self.truckLifetime=self.TransportTab[name]['truckLifetime']
        self.truckHours=self.TransportTab[name]['truckHours']
        self.pipeSystem=self.TransportTab[name]['pipeSystem']
        self.truckOMfix=self.TransportTab[name]['truckOMfix']
        self.truckDriver=self.TransportTab[name]['truckDriver']
        self.truckFuelDemand=self.TransportTab[name]['truckFuelDemand']
        self.trailerInvest=self.TransportTab[name]['trailerInvest']
        self.trailerLifetime=self.TransportTab[name]['trailerLifetime']
        self.trailerHours=self.TransportTab[name]['trailerHours']
        self.trailerOM=self.TransportTab[name]['trailerOM']
        self.trailerCapacity=self.TransportTab[name]['trailerCapacity']
        #self.speed=self.TransportTab[name]['speed']
        self.loadingtime=self.TransportTab[name]['loadingtime']
        self.costH2In=costH2In
        self.boilOffHourly=self.TransportTab[name]['boilOffHourly']  
        self.pressureStationMin=self.TransportTab[name]['pipePressureStation'] 
        self.pressureHubMin=self.TransportTab[name]['pipePressureHub']
        nameStation=dfHSC['Station'][i]
        self.stationDemand=np.ones((len(distance),demand.shape[1]))*dfTable['Station'][nameStation]['stationDemand']*self.utilization
        self.stationDistance=np.ones((len(distance)))*self.distributionDistance
        
    def getTOTEX(self):
        if self.pipeSystem==1: ###Pipeline
                    #Transmission Pipeline Storage - Hub
            #Diameter calculation for the pipeline--> minimal diameter 100
            self.pipeDiameterTrans,self.pressureHub, self.speedTrans=aFun.getDiaPress(self.demand, self.distance, self.pressureIn, self.pressureHubMin)
            #Calculation of Pipeline investcost (€/km)
            self.pipeInvestTransSpec=(self.pipeInvestA*self.pipeDiameterTrans**2+self.pipeInvestB*self.pipeDiameterTrans+self.pipeInvestC)*1000
            #Investkosten Transmission
            self.pipeTransInvest=self.distance*self.pipeInvestTransSpec
            
                    ###Distribution Pipeline Hub - Station ###
            #Diameter calculation for the pipeline--> minimal diameter 100
            self.pipeDiameterDist,self.pressureOut, self.speedDist=aFun.getDiaPress(self.stationDemand, self.stationDistance, self.pressureHub, self.pressureStationMin)
            #Calculation of Pipeline investcost (€/km)
            self.pipeInvestDistSpec=(self.pipeInvestA*self.pipeDiameterDist**2+self.pipeInvestB*self.pipeDiameterDist+self.pipeInvestC)*1000     
            self.pipeDistInvest=self.stationDistance*self.pipeInvestTransSpec
            
            #annuity calculation for pipeline
            self.pipeAnnuity=self.getAnnuity(self.WACC,self.pipeLifetime)
            self.pipeCAPEXTrans=self.pipeAnnuity*self.pipeTransInvest/(self.pipeHours/24*self.demand)           
            self.pipeCAPEXDist=self.pipeAnnuity*self.pipeDistInvest/(self.pipeHours/24*self.stationDemand)            
            
            self.pipeOPEXTrans=self.pipeOM*self.pipeTransInvest/(self.pipeHours/24*self.demand)
            self.pipeOPEXDist=self.pipeOM*self.pipeDistInvest/(self.pipeHours/24*self.stationDemand)
            #CAPEX
            self.CAPEX=self.pipeCAPEXTrans+self.pipeCAPEXDist
            self.fixOPEX=self.pipeOPEXDist+self.pipeOPEXTrans
            self.varOPEX=0
            print()
#            self.TOTEXTrans=self.pipeCAPEXTrans+self.pipeOPEXTrans
#            self.TOTEXDist=self.pipeCAPEXDist+self.pipeOPEXDist
            self.boilOff=0
            self.truckFuelUse=0
            
        else:     ###Truck
            #Annuity for the trailer truck
            self.truckAnnuity=self.getAnnuity(self.WACC,self.truckLifetime)
            #annuity calculation for the trailer
            self.trailerAnnuity=self.getAnnuity(self.WACC,self.trailerLifetime)        
            #time for traveling to destination and back with loading
            self.truckTime=(self.distance/self.speed+self.loadingtime)*2
            #Used fuel for traveling to destination and back
            self.truckFuelUse=self.truckFuelDemand*self.distance*2/100
        
            #Additional Calculations
            self.truckHourly=(self.truckAnnuity+self.truckOMfix)*self.truckInvest/self.truckHours
            self.trailerHourly=(self.trailerAnnuity+self.trailerOM)*self.trailerInvest/self.trailerHours
            #OPEX calculation in €/kg
            self.truckCAPEX=self.truckAnnuity*self.truckInvest/self.truckHours*self.truckTime/self.trailerCapacity
            self.trailerCAPEX=self.trailerAnnuity*self.trailerInvest/self.trailerHours*self.truckTime/self.trailerCapacity
            #OPEX calculation in €/kg        
            self.truckOPEXHours=(self.truckOMfix*self.truckInvest/self.truckHours+self.driverCost*self.truckDriver)*self.truckTime/self.trailerCapacity
            self.truckOPEXFuel=self.truckFuelUse*self.dieselCost/self.trailerCapacity
            self.truckOPEX=self.truckOPEXHours+self.truckOPEXFuel
            self.trailerOPEX=self.trailerOM*self.trailerInvest/self.trailerHours*self.truckTime/self.trailerCapacity
            self.boilOff=self.boilOffHourly*self.truckTime/2
        
        #Total CAPEX OPEX and TOTEX
            self.varOPEX=self.truckOPEXFuel+self.costH2In*self.boilOff
            self.CAPEX=self.truckCAPEX+self.trailerCAPEX
            self.fixOPEX=self.truckOPEXHours
            self.OPEX=self.truckOPEX+self.trailerOPEX+self.costH2In*self.boilOff
            
        
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        self.electricityDemandAnnual=0
        ### Total Cost H2
        self.costH2Out=self.costH2In+self.TOTEX
        
        return self.TOTEX

    def getTotalCost(self):
        self.costH2Out=self.costH2In+self.getTOTEX()
        return self.costH2Out

    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen, ElectricityRES, electricityGrid, Heat, Diesel
        return 1+self.boilOff, 0, 0., 0., self.truckFuelUse/self.trailerCapacity*11#kWh/l       

    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX
#%% Station Instance
class Station(Module):
    def __init__(self, demand, distance,name,dfTable, costH2In, i, dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):
        
        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Storage Technology
        '''
        #Implementation of Motherclass
            #Getting the General Informations
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        #Import Transport Table
        self.StationTab=dfTable['Station']
        
        self.form = self.StationTab[name]['form']
        #self.pressureIn = self.StationTab[name]['pressureIn']
        self.stationDemand=self.StationTab[name]['stationDemand']
        self.stationInvest=self.StationTab[name]['stationInvest']
        self.stationLifetime=self.StationTab[name]['stationLifetime']
        self.stationOM=self.StationTab[name]['stationOM']
        self.electricityDemand=self.StationTab[name]['electricityDemand']
        self.heatDemand=self.StationTab[name]['heatDemand']
        self.boilOffEff=self.StationTab[name]['boilOffEff']
        self.dieselDemand=self.StationTab[name]['dieselDemand']        
        self.costH2In=costH2In
        
        if name=='GStation':
            self.trailerInvest=dfTable['Transport'][self.nameTransport]['trailerInvest']
            self.trailerLifetime=dfTable['Transport'][self.nameTransport]['trailerLifetime']
            self.trailerOM=dfTable['Transport'][self.nameTransport]['trailerOM']
            self.trailerAnnuity=self.getAnnuity(self.WACC,self.trailerLifetime)
        else:
            self.trailerAnnuity=0
            self.trailerInvest=0
            self.trailerOM=0
                              
            
            
    def getTOTEX(self):         
        #Annuity Calculation
        self.annuity=self.getAnnuity(self.WACC,self.stationLifetime)
        #CAPEX calculation
        self.CAPEX=(self.annuity*self.stationInvest+self.trailerAnnuity*self.trailerInvest)/(self.stationDemand*self.utilization)/365
        #OPEX calculation
        self.fixOPEX=(self.stationOM*self.stationInvest+self.trailerOM*self.trailerInvest)/(self.stationDemand*self.utilization)/365
        self.electricityDemandAnnual=np.array([self.electricityDemand*self.stationDemand*self.utilization*365])
        #variable electricity price - uncomment for predefined
        self.electricityCost=self.getElectricityCost(self.electricityDemandAnnual)
        
        self.varOPEX=self.electricityDemand*self.electricityCost+self.heatDemand*self.NGCost++self.boilOffEff*self.costH2In+self.dieselDemand*self.dieselCost/11
        #TOTEX calculation
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        

        return self.TOTEX
        
    def getTotalCost(self):
        self.costH2Out=self.costH2In+self.getTOTEX()
        return self.costH2Out   
        
    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen, ElectricityRES, electricityGrid, Heat, Diesel
        return 1+self.boilOffEff,0. , self.electricityDemand, self.heatDemand, self.dieselDemand      

    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX
#%%Connector Class        
class Connector(Module):
    def __init__(self, demand, distance, dfTable, costH2In, i, dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):

        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Connection Technology
        position: Connected Modules
        '''
        #Implementation of Motherclass
        #Getting the General Informations
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        self.costH2In=costH2In
        self.name=dfHSC['Connector1'][i]
        self.dfTable=dfTable

        
        #Import from Inlet
        self.formIn = dfTable['Production'][self.nameProduction]['form']
        self.pressureIn = dfTable['Production'][self.nameProduction]['pressureOut']
        
        #Import from Outlet
        self.formOut = dfTable['Storage'][self.nameStorage]['form']
        self.pressureOut = dfTable['Storage'][self.nameStorage]['pressureIn']   
        
        
    def getValues(self, name):
        '''
        Import the Values of the specific Connector
        '''
                
        #Import Connector Table
        self.ConnectorTab=self.dfTable['Connector'][name]
                #Parameter Import
        self.investBase=self.ConnectorTab['investBase']
        self.investCompare=self.ConnectorTab['investCompare']
        self.investScale=self.ConnectorTab['investScale']
        self.investLifetime=self.ConnectorTab['investLifetime']
        self.investOM=self.ConnectorTab['investOM']
        self.electricityDemandBase=self.ConnectorTab['electricityDemandBase']
        self.electricityDemandScale=self.ConnectorTab['electricityDemandScale']
        self.electricityDemandCompare=self.ConnectorTab['electricityDemandCompare']
        self.heatDemand=self.ConnectorTab['heatDemand']
        self.heatSupply=self.ConnectorTab['heatSupply']
        self.boilOffEff=self.ConnectorTab['boilOffEff']
        self.system=self.ConnectorTab['system']

        
    def getTOTEX(self):
      
        self.getValues(self.name)
        
        self.overCapacity=((1+self.boilOffEff)
            *(1+self.dfTable['Connector'][self.nameC2]['boilOffEff'])
            *(1+self.dfTable['Station'][self.nameStation]['boilOffEff'])
            *(1+self.dfTable['Storage'][self.nameStorage]['boilOff']*self.storageDays))
        
        #annual production
        self.annualProduction=self.demand*self.overCapacity*365
        #annuity calculation
        self.annuity=self.getAnnuity(self.WACC, self.investLifetime)                
        
        # If Connector == Compressor then the compression systematic is used
        if self.name=='Compressor':
            self.electricityDemand, self.electricityPower=aFun.getCompressionEnergy(self.pressureIn, self.pressureOut, self.demand)            
            #CAPEX
            self.CAPEX=(self.annuity*self.investBase/self.annualProduction
                *(self.electricityPower*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
            #OPEX
            self.fixOPEX=(self.investOM*self.investBase/self.annualProduction*
                (self.electricityPower*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
        else:#if no compression --> liquid or LOHC is used
            self.electricityDemand=self.electricityDemandBase*(self.demand/self.electricityDemandCompare)**self.electricityDemandScale
            
            #CAPEX
            self.CAPEX=(self.annuity*self.investBase/self.annualProduction
                *(self.demand*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
            #OPEX
            self.fixOPEX=(self.investOM*self.investBase/self.annualProduction*
                (self.demand*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
            
            if self.system==31: #in case of dehydrogenation we need additional compression!!!
                self.getValues('Compressor')
                self.electricityDemandComp, self.electricityPowerComp=aFun.getCompressionEnergy(self.pressureIn, self.pressureOut, self.demand)            
                self.electricityDemand+=self.electricityDemandComp
                #CAPEX
                self.CAPEXComp=(self.annuity*self.investBase/self.annualProduction
                    *(self.electricityPowerComp*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
                self.CAPEX+=self.CAPEXComp
                #OPEX
                self.OPEXComp=(self.investOM*self.investBase/self.annualProduction*
                    (self.electricityPowerComp*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
                self.fixOPEX+=self.OPEXComp
        

        self.varOPEX=self.storagePart*(self.electricityDemand*self.electricityCostRES+self.heatDemand*self.NGCost-self.heatSupply*self.heatGain+self.boilOffEff*self.costH2In)
        
        #TOTEX
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        self.costH2Out=self.costH2In+self.TOTEX
        
        return self.TOTEX
        
    def getTotalCost(self):
        self.costH2Out=self.costH2In+self.getTOTEX()
        return self.costH2Out        

    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen, ElectricityRES, electricityGrid, Heat, Diesel
        return 1+self.boilOffEff, self.electricityDemand,0. , self.heatDemand, 0.             
     
    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX
        
#%%
class Connector2(Module):
    def __init__(self, demand, distance, dfTable, costH2In, i, dfHSC):
        '''
        transport Module:
        Inputs (like Motherclass):

        demand: demand of hydrogen per day
        distance: distance between hydrogen source and sink
        name: Kind of Connection Technology
        position: Connected Modules
        '''
        #Implementation of Motherclass
        #Getting the General Informations
        Module.__init__(self, demand, distance,dfTable, i, dfHSC)
        self.costH2In=costH2In
        self.name=dfHSC['Connector2'][i]
        self.dfTable=dfTable
        
        #Import from Inlet
        self.formIn = dfTable['Storage'][self.nameStorage]['form']
        self.pressureIn = dfTable['Storage'][self.nameStorage]['pressureOut']
        
        #Import from Outlet
        self.formOut = dfTable['Transport'][self.nameTransport]['form']
        self.pressureOut = dfTable['Transport'][self.nameTransport]['pressureIn']   
            
    def getValues(self, name, kind):
        '''
        Import the Values of the specific Connector
        '''
                
        #Import Connector Table
        self.ConnectorTab=self.dfTable['Connector'][name]
                #Parameter Import
        self.investBase=self.ConnectorTab['investBase']
        self.investCompare=self.ConnectorTab['investCompare']
        self.investScale=self.ConnectorTab['investScale']
        self.investLifetime=self.ConnectorTab['investLifetime']
        self.investOM=self.ConnectorTab['investOM']
        self.electricityDemandBase=self.ConnectorTab['electricityDemandBase']
        self.electricityDemandScale=self.ConnectorTab['electricityDemandScale']
        self.electricityDemandCompare=self.ConnectorTab['electricityDemandCompare']
        self.system=self.ConnectorTab['system']
        if kind==0:        
            self.heatDemand=self.ConnectorTab['heatDemand']
            self.heatSupply=self.ConnectorTab['heatSupply']
            self.boilOffEff=self.ConnectorTab['boilOffEff']

        
    def getTOTEX(self):
        self.getValues(self.name, kind=0)
        
        self.overCapacity=((1+self.boilOffEff)
            *(1+self.dfTable['Connector'][self.nameC2]['boilOffEff'])
            *(1+self.dfTable['Station'][self.nameStation]['boilOffEff'])
            *(1+self.dfTable['Storage'][self.nameStorage]['boilOff']*self.storageDays))
        
        #annual production
        self.annualProduction=self.demand*self.overCapacity*365
        #annuity calculation
        self.annuity=self.getAnnuity(self.WACC, self.investLifetime)                
        
        if self.name=='Compressor':
            self.electricityDemand, self.electricityPower=aFun.getCompressionEnergy(self.pressureIn, self.pressureOut, self.demand)            
            #CAPEX
            self.CAPEX=(self.annuity*self.investBase/self.annualProduction
                *(self.electricityPower*self.overCapacity/self.investCompare)**self.investScale)
            #OPEX
            self.fixOPEX=(self.investOM*self.investBase/self.annualProduction*
                (self.electricityPower*self.overCapacity/self.investCompare)**self.investScale)
        else:
            self.electricityDemand=self.electricityDemandBase*(self.demand/self.electricityDemandCompare)**self.electricityDemandScale
            
            #CAPEX
            self.CAPEX=(self.annuity*self.investBase/self.annualProduction
                *(self.demand*self.overCapacity*self.overDimension/self.investCompare)**self.investScale)
            #OPEX
            self.fixOPEX=(self.investOM*self.investBase/self.annualProduction*
                (self.demand*self.overCapacity/self.investCompare)**self.investScale)
            
            if self.system==31: #in case of dehydrogenation we need additional compression!!!
                kind=1
                self.getValues('Compressor', kind)
                self.electricityDemandComp, self.electricityPowerComp=aFun.getCompressionEnergy(self.pressureIn, self.pressureOut, self.demand)            
                self.electricityDemand+=self.electricityDemandComp
                #CAPEX
                self.CAPEXComp=(self.annuity*self.investBase/self.annualProduction
                    *(self.electricityPowerComp*self.overCapacity/self.investCompare)**self.investScale)
                self.CAPEX+=self.CAPEXComp
                #OPEX
                self.OPEXComp=(self.investOM*self.investBase/self.annualProduction*
                    (self.electricityPowerComp*self.overCapacity/self.investCompare)**self.investScale)
                self.fixOPEX+=self.OPEXComp
        
        
        #Variable electricity price industry - uncomment to work with predefined
        self.annualElectricityDemand=self.electricityDemand*365*self.demand
        self.electricityCost=self.getElectricityCost(self.annualElectricityDemand)
        
        
        
        self.varOPEX=self.storagePart*(self.electricityDemand*self.electricityCost
                        +self.heatDemand*self.NGCost
                        -self.heatSupply*self.heatGain
                        +self.boilOffEff*self.costH2In)
        
        #TOTEX
        self.TOTEX=self.CAPEX+self.fixOPEX+self.varOPEX
        self.costH2Out=self.costH2In+self.TOTEX
        
        return self.TOTEX
        
    def getTotalCost(self):
        self.costH2Out=self.costH2In+self.getTOTEX()
        return self.costH2Out        

    def getDemand(self):
        self.getTOTEX()
        #Output: Hydrogen, ElectricityRES, electricityGrid, Heat, Diesel
        return 1+self.boilOffEff,0, self.electricityDemand, self.heatDemand, 0

    def getExpenditures(self):
        return self.CAPEX, self.fixOPEX, self.varOPEX