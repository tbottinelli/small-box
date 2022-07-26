#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2011  Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#

import argparse
import h5py
import os
from numpy import *
#from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def main():
    # define and parse command line arguments
    parser = argparse.ArgumentParser(prog='plot_pot.py')
    parser.add_argument('--range', type=int, nargs=2, help='select range of data points')
    parser.add_argument('--dump', metavar='FILENAME', help='dump plot data to filename')
    parser.add_argument('--no-plot', action='store_true', help='do not produce plots, but do the analysis')
    parser.add_argument('--group', help='particle group (default: %(default)s)', default='all')
    parser.add_argument('input1', metavar='INPUT', help='H5MD input file with data for state variables')
    parser.add_argument('input2', metavar='INPUT', help='H5MD input file with data for state variables')
    parser.add_argument('input3', metavar='INPUT', help='H5MD input file with data for state variables')#
    parser.add_argument('input4', metavar='INPUT', help='H5MD input file with data for state variables')

	
    args = parser.parse_args()
    s0=0
    TR=30
    AT=25
    a=0
    Delta=7.5
    x=1

    H5 = [ h5py.File(args.input1, 'r'), h5py.File(args.input2, 'r'), h5py.File(args.input3, 'r'), h5py.File(args.input4, 'r')]
    Epot = []
    for j in range(len(H5)):
        H5obs = H5[j]['observables']
        Epot.append( [H5obs['region{}/potential_energy/value'.format(i)][x:] for i in range(0,40)] )
    Epot = np.array(Epot)
    mean_epot = np.mean(Epot, axis = 2)
    #print(mean_epot)
    #templist = mean_temp/ Temperature.shape[0]
    #print(templist)
    dx =  2.5
	

    plt.rc('font', **{ 'family':'serif', 'serif' : ['ptm'], 'size' :12})
    plt.rc('text', usetex=True)
    plt.rc('text.latex' , preamble=(
        r'\usepackage{textcomp}',
        r'\usepackage{amsmath}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage{times}',
        r'\usepackage{txfonts}',
        ))
    plt.rc('legend', frameon=False, numpoints=1, fontsize=8, labelspacing=0.2, handlelength=2, handletextpad=0.5, borderaxespad=0.5)
    plt.rc('figure',figsize=(4.7,2))
    plt.rc('xtick', direction='in',top=True)
    plt.rc('ytick', direction='in',right=True)
    plt.rc('xtick.minor',visible=True,top=True)
    plt.rc('ytick.minor',visible=True,right=True)
    plt.rc('axes', linewidth=0.7 )
    plt.rc('lines', linewidth=1, markersize = 2,markeredgewidth=0)
    plt.rc('savefig', bbox='tight',pad_inches=0.05,dpi=600,transparent=False)
    plt.rc('ps',usedistiller='xpdf')
     
    xgrid = dx * np.arange(int((2*TR + AT + 2.0*Delta)/dx))+1.25 
   # print(xgrid)
        
    rdf0 = interp1d(xgrid, mean_epot[0,:] ,bounds_error=False, kind = 'quadratic')
    rdf1 = interp1d(xgrid, mean_epot[1,:] ,bounds_error=False, kind = 'quadratic')
    rdf2 = interp1d(xgrid, mean_epot[2,:] ,bounds_error=False, kind = 'quadratic')
    rdf3 = interp1d(xgrid, mean_epot[3,:] ,bounds_error=False, kind = 'quadratic')

    grids_at = np.linspace(0, 70, num = 55, endpoint = False )
    grids_adr = np.linspace(0,100, num =1000 , endpoint=False)
    sym = -50

    plt.plot(grids_adr + sym, rdf0(grids_adr) , '-',color='deepskyblue',linewidth=1.2,fillstyle='full', label = 'D8')
    plt.plot(grids_adr + sym, rdf1(grids_adr) , '-',color='royalblue',linewidth=1.2,fillstyle='full', label = 'D10')
    plt.plot(grids_adr + sym, rdf2(grids_adr) , '-',color='mediumblue',linewidth=1.2,fillstyle='full', label = 'D15')
    plt.plot(grids_adr + sym, rdf3(grids_adr) , '-',color='midnightblue',linewidth=1.2,fillstyle='full', label = 'D20')

    plt.legend()

    plt.xlabel(r"$x / \sigma$")
    plt.ylabel(r"$k_{B}T(x)/\varepsilon$") 
    #print(np.max(time0)-0.8, np.max(time1)-1,np.max(time2)-1.2)
    hm = 30
    source = 10
    slab = 20
    plt.xlim([0+sym,100+sym])
    #plt.ylim([0,5])
    #plt.legend(loc = 'upper left')
    plt.axvline(x=source+sym, color='k', linestyle='--',linewidth=0.4)
 

    plt.axvspan(-slab/2,slab/2, alpha=0.5, color='grey')
    plt.axvspan(0+sym, source+sym, alpha=0.5, color='gold')

    plt.savefig('epot.pdf')
   



if __name__ == '__main__':
    main()
