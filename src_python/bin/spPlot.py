#!/usr/bin/python3
# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################


import os, errno
import yaml
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from numpy import pi, real, imag
from PlotUtils import plotArray, plotArrayResponse, plotArrayResponseCuts







def main():

    
    plt.ion()
    
    savePlotsToFile = True
    outputDirPlot = './plots'
    imageFormat = 'png'
    imageDPI = 300
    YAML_FILE = 'Optimized_Info.yaml'
    
    half = True
    
    
    try:
        os.makedirs(outputDirPlot)
    except OSError as exception:
        if exception.errno != errno.EEXIST:    raise


    try:
        with open(YAML_FILE) as f:
            ArrayInfo = yaml.load(f, Loader=yaml.FullLoader)
            f.close()
    except yaml.scanner.ScannerError as e:
        print('There is a syntax problem with the YAML file {0}.'.format(YAML_FILE))
        print(type(e))
        print(e)
        return
    except Exception as e:
        print('There was a problem while loading the YAML file {0}.'.format(YAML_FILE))
        print(type(e))
        print(e)
        return
        
    try:
        Kmin = ArrayInfo.get('Kmin')
        Kmax = ArrayInfo.get('Kmax')        
        Nsensors = ArrayInfo.get('Nsensors') 
        normalized = ArrayInfo.get('normalized', False)        
    except Exception as e:
        print('There was a problem while parsing the configuration file {0}'.format(YAML_FILE))
        print(type(e))
        print(e)
        return

    try:
        pos = np.atleast_2d(np.loadtxt('PossiblePositions.csv', comments='#', delimiter='\t'))
        CostMask = np.atleast_2d(np.loadtxt('CostMask.csv', comments='#', delimiter='\t'))
        pos_sol = np.atleast_2d(np.loadtxt('Optimized_ArrayLayout.csv', comments='#', delimiter='\t'))
    except Exception as e:
        print('There was a problem while loading a CSV file.')
        print(type(e))
        print(e)
        return
    
    
        
    
    ### Plotting
    
    fig = plotArray(pos_sol)
    plt.plot(pos[:,0], pos[:,1], 'o', ms=4, lw=0, alpha=1.0, mfc='blue')
    ax = fig.gca()

    xMin = np.min(pos[:,0]); xMax = np.max(pos[:,0]);
    yMin = np.min(pos[:,1]); yMax = np.max(pos[:,1]);
    xAperture = np.max([xMax - xMin, 0.8*(yMax - yMin)]) # avoid really stretched aspect ratios
    yAperture = np.max([yMax - yMin, 0.8*(xMax - xMin)])
    ax.set_xlim(xMin-0.05*xAperture, xMax+0.05*xAperture)
    ax.set_ylim(yMin-0.05*yAperture, yMax+0.05*yAperture)
    if savePlotsToFile: plt.savefig('{0}/Optimized_ArrayLayout.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)
    
    
    plotArrayResponse(pos_sol, 1.5*Kmax, Knum=200, normalized=normalized)
    plt.plot(CostMask[:,0], CostMask[:,1], 'o', ms=3, lw=0, alpha=0.8, mfc='blue')
    Thetavec = np.linspace(0, 2*np.pi, 180)
    plt.plot(Kmin*np.cos(Thetavec), Kmin*np.sin(Thetavec),'b:',linewidth=2)
    plt.plot(Kmax*np.cos(Thetavec), Kmax*np.sin(Thetavec),'b:',linewidth=2)
    if savePlotsToFile: plt.savefig('{0}/Optimized_ArrayResponse.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)    
    
    plotArrayResponseCuts(pos_sol, 1.5*Kmax, Knum=300, Thetanum=90, normalized=normalized, half=half)
    h = 1 if normalized else Nsensors**2
    plt.plot([Kmin, Kmin], [0, h],'b:',linewidth=2)
    plt.plot([Kmax, Kmax], [0, h],'b:',linewidth=2)
    CostMask_K = np.unique(np.sqrt(CostMask[:,0]**2 + CostMask[:,1]**2))
    plt.plot(CostMask_K, (h/2)*np.ones((len(CostMask_K),)), 'o', ms=3, lw=0, alpha=0.8, mfc='blue')
    if not half:
        plt.plot([-Kmin, -Kmin], [0, h],'b:',linewidth=2)
        plt.plot([-Kmax, -Kmax], [0, h],'b:',linewidth=2)
    if savePlotsToFile: plt.savefig('{0}/Optimized_ArrayResponseCuts.{1}'.format(outputDirPlot, imageFormat), dpi=imageDPI)
    
    
    
    plt.waitforbuttonpress(timeout=15)
    
    return



if  __name__ =='__main__':
    main()

