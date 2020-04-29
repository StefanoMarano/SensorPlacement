# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com     #
##################################################
#
# These plotting routines are borrowed from WaveDec
#

import matplotlib.pyplot as plt
import numpy as np
from numpy.matlib import repmat


def plotArray(array_info, SourcePosition=None):
    """
    Plots array layout
    """
    
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    
    xMin = np.min(pos_x)
    xMax = np.max(pos_x)
    yMin = np.min(pos_y)
    yMax = np.max(pos_y)
    if SourcePosition is not None:
        xMin = np.min([xMin, SourcePosition[0]])
        xMax = np.max([xMax, SourcePosition[0]])
        yMin = np.min([yMin, SourcePosition[1]])
        yMax = np.max([yMax, SourcePosition[1]])
        
    xAperture = np.max([xMax - xMin, 0.8*(yMax - yMin)]) # avoid really stretched aspect ratios
    yAperture = np.max([yMax - yMin, 0.8*(xMax - xMin)])
   
    fig = plt.figure()
    ax = fig.gca()
    ax.set_aspect('equal')
    ax.grid()
    ax.plot(pos_x, pos_y, 'o', ms=10, lw=2, alpha=0.7, mfc='orange')
    if SourcePosition is not None:
        ax.plot(SourcePosition[0], SourcePosition[1], '*', ms=18, lw=2, alpha=0.8, mfc='red')
    ax.set_xlim(xMin-0.05*xAperture, xMax+0.05*xAperture)
    ax.set_ylim(yMin-0.05*yAperture, yMax+0.05*yAperture)
    ax.set_xlabel('Easting [m]')
    ax.set_ylabel('Northing [m]')
    plt.show()
    
    return fig
    
    
def plotArrayResponse(array_info, Kmax, Knum=100, normalized=True):
    """
        plot the squared modulo of the array response H^2
        
        pos the sensor positions
        Kmax the largest wavenumber of interest [1/m]
        Knum grid size
        normalized if True, sets the largest value to one
    """
    
    Nsensors = np.shape(array_info)[0]
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    K_vec = np.linspace(-Kmax, Kmax, num=Knum)
    K_plane = repmat(K_vec, Knum, 1)
    H = np.zeros((Knum,Knum))  # the squared modulos is real
    
    # loop solely on the elements above diagonal
    for n1 in range(0,Nsensors):
        for n2 in range(n1+1, Nsensors):
            H += np.cos(2*np.pi*(K_plane*(pos_x[n1]-pos_x[n2]) + np.transpose(K_plane*(pos_y[n1]-pos_y[n2]))))
    
    H *= 2        # elements below diagonal
    H += Nsensors # elements on the diagonal (all equal to 1)

    if normalized:
        H = H / Nsensors**2
        
        
       
    fig = plt.figure()
    plt.pcolor(K_vec, K_vec, H, cmap='hot_r')
    cbar = plt.colorbar()

    if normalized:
        plt.clim(0,1)
        cbar.ax.set_ylabel('Normalized array response')
    else:
        plt.clim(0,Nsensors**2)
        cbar.ax.set_ylabel("Array response")
    plt.xlim(-Kmax, Kmax)
    plt.ylim(-Kmax, Kmax)
    
    plt.axes().set_aspect('equal')
    plt.xlabel('Wavenumber x [1/m]')
    plt.ylabel('Wavenumber y [1/m]')
    plt.show()
    return fig
    
def plotArrayResponseCuts(array_info, Kmax, Knum=200, Thetanum=60, half=False, normalized=True):
    """
        plot sections the squared modulo of the array response H^2
        
        pos the sensor positions
        Kmax the largest wavenumber of interest [1/m]
        Knum grid size
        Thetanum number of sections to be displayed
        normalized if True, sets the largest value to one
    """
    Nsensors = np.shape(array_info)[0]
    pos_x = array_info[:,0]
    pos_y = array_info[:,1]
    if half:
        K_vec = np.linspace(0, Kmax, num=Knum)
    else:
        K_vec = np.linspace(-Kmax, Kmax, num=Knum)
    Theta_vec = np.linspace(0, np.pi, num=Thetanum, endpoint=False)
    H = np.zeros((Knum,Thetanum))  # the squared modulos is real
    
    for tt in range(0,Thetanum):
        Kx_vec = K_vec * np.cos(Theta_vec[tt])
        Ky_vec = K_vec * np.sin(Theta_vec[tt])
        
        # loop solely on the elements above diagonal
        for n1 in range(0,Nsensors):
            for n2 in range(n1+1, Nsensors):
                H[:,tt] += np.cos(2*np.pi*(Kx_vec*(pos_x[n1]-pos_x[n2]) + np.transpose(Ky_vec*(pos_y[n1]-pos_y[n2]))))
    
    H *= 2        # elements below diagonal
    H += Nsensors # elements on the diagonal (all equal to 1)

    if normalized:
        H = H / Nsensors**2
        
       
    fig = plt.figure()
    plt.plot(K_vec, H, color='red', linewidth=1.0, alpha=0.5)
    if normalized:
        plt.ylim(0,1)
        plt.ylabel('Normalized array response')
    else:
        plt.ylim(0,Nsensors**2)
        plt.ylabel('Array response')
    plt.xlim(np.min(K_vec), np.max(K_vec))
    
    plt.xlabel('Wavenumber [1/m]')
    plt.show()
    return fig
