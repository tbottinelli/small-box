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
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def main():
    # define and parse command line arguments
    parser = argparse.ArgumentParser(prog='plot_density.py')
    parser.add_argument('--range', type=int, nargs=2, help='select range of data points')
    parser.add_argument('--dump', metavar='FILENAME', help='dump plot data to filename')
    parser.add_argument('--no-plot', action='store_true', help='do not produce plots, but do the analysis')
    parser.add_argument('--group', help='particle group (default: %(default)s)', default='all')
    parser.add_argument('input1', metavar='INPUT', help='H5MD input file with data for state variables')
    #parser.add_argument('input2', metavar='INPUT', help='H5MD input file with data for state variables')
    #parser.add_argument('input3', metavar='INPUT', help='H5MD input file with data for state variables')
    #parser.add_argument('input4', metavar='INPUT', help='H5MD input file with data for state variables')


    args = parser.parse_args()

    # open and read data file
    H5 =  h5py.File(args.input1, 'r')#, h5py.File(args.input2, 'r'), h5py.File(args.input3, 'r'), h5py.File(args.input4, 'r') ]
   # density = []
   # for j in range(len(H5)):
   #     H5obs = H5[j]['observables']
   #     density.append([ H5obs['region{0}/particle_number/value'.format(i)] for i in range(0,40)] )

    H5obs = H5['observables']
    density = [ H5obs['region{0}/particle_number/value'.format(i)] for i in range(0,16)  ]
    density = np.array(density)/(10*25*2.5)
    density_mean = np.mean(density, axis=1)
    print(density_mean.shape) #should be 1x40 in the end

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

    dx=2.5
    xgrid = dx * np.arange(int((40)/dx))+1.25 

    # select data for plotting the potential energy as function of time
   # y0 = interp1d(xgrid, density_mean[0,:], bounds_error=False, kind='quadratic')
   # y1 = interp1d(xgrid, density_mean[1,:], bounds_error=False, kind='quadratic')
    y2 = interp1d(xgrid, density_mean[:], bounds_error=False, kind='quadratic')
   # y3 = interp1d(xgrid, density_mean[3,:], bounds_error=False, kind='quadratic')

    grids_at = np.linspace(0, 70, num = 55, endpoint = False )
    grids_adr = np.linspace(0,40, num =500 , endpoint=False)
    sym = -20



    # generate plot
    if not args.no_plot:
        # plot potential energy versus time
    #    plt.plot(grids_adr+sym, y0(grids_adr), '-',  color = 'deepskyblue' , label='D8')
     #   plt.plot(grids_adr+sym, y1(grids_adr),'-', color = 'royalblue' , label='D10')
        plt.plot(grids_adr+sym, y2(grids_adr),'-', color = 'mediumblue' , label='D15')
      
      #plt.plot(grids_adr+sym, y3(grids_adr),'-', color = 'midnightblue' , label='D20')


        # add axes labels and finalise plot
        axis('tight')
        xlabel(r'x')
        ylabel(r'density profile')
        plt.legend()
        sym = -20
        source = 3
        slab = 18.6
        plt.xlim([0+sym,40+sym])
        plt.axvline(x=source+sym, color='k', linestyle='--',linewidth=0.4)
        plt.axvspan(-slab/2,slab/2, alpha=0.5, color='grey')
        plt.axvspan(0+sym, source+sym, alpha=0.5, color='gold')

        plt.savefig('density_smallbox.pdf')
   


def compute_average(data, label, nblocks = 10):
    """ compute and print average of data set
        The error is estimated from grouping the data in shorter blocks """

    # group data in blocks, discard excess data
    data = reshape(data[:(nblocks * (data.shape[0] / nblocks))], (nblocks, -1))

    # compute mean and error
    avg = mean(data, axis=0)
    #err = std(mean(data, axis=1)) / sqrt(nblocks - 1)
    #print '{0:s}: {1:.4f} ± {2:.2g}'.format(label, avg, err)

    # return as tuple
    return avg

if __name__ == '__main__':
    main()
