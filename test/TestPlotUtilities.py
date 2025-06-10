# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 15:12:06 2025

@author: aantonak
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.ftplotlib import FTPlotlib as FTP

df = pd.read_csv('test/data/FlightData.csv',sep=';')

fig,ax = plt.subplots()
    

Ts = df['time (s)']

Fdr = FTP.FTPlot(fig,ax,slider=True,XLim=[0,int(np.max(Ts))+1],Ngridy = 16)

Fdr.AddAxis(Name='Euler',GridHeight=4,GridPos=12,Unit='deg')
Fdr.AddAxis(Name='AoA',GridHeight=4,GridPos=9,Unit='deg')
Fdr.AddAxis(Name='Sideslip',GridHeight=4,GridPos=9,Unit='deg',Position='Right')
Fdr.AddAxis(Name='Rates',GridPos=6,GridHeight=4,Unit='deg/s',Position='Right')
Fdr.AddAxis(Name='Controls',GridHeight=4,GridPos=3,Unit='deg')
Fdr.AddAxis(Name='Altitude',GridPos=0,Unit='m')
Fdr.AddAxis(Name='Speed',GridPos=0,Unit='m/s',Position = 'Right')


AxisCurves = {'Euler':['$\\phi$','$\\theta$','$\\psi$'],
              'AoA':['$\\alpha$'],
              'Sideslip' : ['$\\beta$'],
              'Rates' : ['P','Q','R'],
              'Controls' : ['$\\delta l$','$\\delta m$','$\\delta n$'],
              'Altitude' : ['$Z_p$','$Z_e$'],
              'Speed' : ['TAS']}



for axis in AxisCurves.keys():

    for Name in AxisCurves[axis]:
        
        Color = np.random.uniform(0,.8,size=3)
        Fdr.AddCurve(Name, axis, Ts, df[Name], color = Color)


plt.show()

