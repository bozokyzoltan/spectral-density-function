#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2012.11.21.

Purpose:
========
TauC and Spectral density function calculation based on R1, R1rho and
heteronuclear NOE experiments.

According to:
=============
Spectral density function maplting using 15N relaxation data exclusively.
Farrow NA, Zhang O, Szabo A, Torchia DA, Kay LE.
J Biomol NMR. 1995 Sep;6(2):153-62. PMID: 8589604

Requires:
=========
* matplotlib module
* NMR experimental data: R1, R1rho and hetNOE

"""

import math
from pylab import plt


class SpectralDensityFunction(object):
    """
    """
    def __init__(self, R1resi, R1value, R2resi, R2value, NOEresi, NOEvalue):
        """ """
        # ----------------------
        # Load experimental data
        self.resi = []
        self.r1   = []
        self.r2   = []
        self.noe  = []
        for element in R1resi:
            if (element in R2resi) and (element in NOEresi):
                self.resi.append(element)
                self.r1.append(R1value[R1resi.index(element)])
                self.r2.append(R2value[R2resi.index(element)])
                self.noe.append(NOEvalue[NOEresi.index(element)])
        # ----------------------
        # Constants declarations
        # ----------------------
        # Giromagnetic ratios:
        GiromagneticRatio_H1  = 2.675222005E+8
        GiromagneticRatio_N15 = -27.116E+6
        GiromagneticRatio_C13 = 67.262E+6
        #----------------------------
        # PI constant
        PI = math.pi
        #----------------------------
        # the permeability of free space = 4*PI*1E-7
        mu0 = 1.2566370614E-6 / (4*PI)
        #----------------------------
        # Planck constant
        h = 6.6260695729E-34 / (2*PI)
        #----------------------------
        # internuclear distances
        rNH = 1.02E-10  # 1H-15N internuclear distance = 1.02 A
        #----------------------------
        # Larmour frequencies - 600 MHz - in rad
        LarmourFrequency_1H  =  2*PI * 598.817 * 1E+6
        LarmourFrequency_2H  =  2*PI * 91.922  * 1E+6
        LarmourFrequency_15N = -2*PI * 60.684  * 1E+6
        LarmourFrequency_13C =  2*PI * 150.589 * 1E+6
        ########################################
        # Note, that d = ( mu0 * h * gammaH * gammaN) / (8*PI**2 )* (1/ rNH**3), BUT 8*PI**2 can be eliminated if the Plack's constant is h/2PI and the permeability is mu/4PI units
        d  = (mu0 * h * GiromagneticRatio_N15 * GiromagneticRatio_H1 )/ (rNH**3)
        d2 = d**2
        #
        wN = LarmourFrequency_15N
        c  = (wN* 170E-6)
        c2 = (2/15)*c**2
        #
        ###################################
        self.result = {}
        for i in range(len(self.resi)):
            r1  = self.r1[i]
            r2  = self.r2[i]
            noe = self.noe[i]
            ###################################

            JwH_0870 = ((r1 * (noe - 1.0)  * 4 * GiromagneticRatio_N15 ) /
                        ( d2 * GiromagneticRatio_H1 * 5)
                       )

            JwH_0921 = (0.870 / 0.921)**2 * JwH_0870
            JwH_0955 = (0.870 / 0.955)**2 * JwH_0870
            # Method 1 = Assuming that J(w) is effectively constant at w ~ wH and therefore J(0.921 wH) and J(0.955 wH) may be replaced by the experimentally determined J(0.870 wH)
            JwN_method1 = (r1 - (JwH_0870* 7 * d2/4.0)) / (c2 + 3 * d2/4.0)
            J0_method1 = (r2 - 3*JwN_method1*(d2/8 + c2/6) - 13*JwH_0870*d2/8) / (4*(d2/8+c2/6))
            # Method 2 = Assume that J(w) ~ 1 / w**2 and estimate the values using J(e*wH) = (0.870/e)**2 * J(0.870 wH)
            JwN_method2 = (r1 - (JwH_0921* 7 * d2/4.0)) / (c2 + 3 * d2/4.0)
            J0_method2 = (r2 - 3*JwN_method2*(d2/8 + c2/6) - 13*JwH_0955*d2/8) / (4*(d2/8+c2/6))
            ###################################

            TauC = (1 + math.sqrt(1 - 4 * JwH_0870**2 * LarmourFrequency_1H**2)) / (2 * JwH_0870 * LarmourFrequency_1H**2)
            ###################################
            self.result[self.resi[i]] = {}

            self.result[self.resi[i]]['J(0.870 wH)'] = JwH_0870
            self.result[self.resi[i]]['J(0.921 wH)'] = JwH_0921
            self.result[self.resi[i]]['J(0.955 wH)'] = JwH_0955

            self.result[self.resi[i]]['J(wN) - lower limit'] = JwN_method1
            self.result[self.resi[i]]['J(wN) - uplter limit'] = JwN_method2

            self.result[self.resi[i]]['J(0) - lower limit'] = J0_method1
            self.result[self.resi[i]]['J(0) - uplter limit'] = J0_method2

            self.result[self.resi[i]]['TauC'] = TauC
        #
        return None
    ### ==================================================================== ###
    def Plot(self, filename = ''):
        """ """
        if filename == '':
            plt.figure(figsize=(10, 5), dpi=72)
        else:
            plt.figure(figsize=(10, 5), dpi=300)
        order = ['Get_J_0870', 'Get_JwN', 'Get_J0']
        names = ['npR_J_0870', 'npR_JwN', 'npR_J0']
        yscale = ['J(0.870 wH)', 'J(wN)', 'J(0)']
        yylim = [4E-11, 5E-10, 7E-9]

        for i,o in enumerate(order):
            plt.subplot(len(order), 1, i+1)
            x = self.residue
            y = eval('self.' + o + '()')

            plt.plot(x, y, marker = 'o')
            self.SaveResults(x, y, names[i] + '_data.txt')

            plt.xlabel('residue #')
            plt.xticks([a for a in range(650, 850, 10)])
            plt.xlim(650, 840)
            plt.ylabel(yscale[i])
            #plt.ylim(ymax = yylim[i])
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename)
        ######################################
        # TauC
        if filename == '':
            plt.figure(figsize=(12, 4), dpi=72)
        else:
            plt.figure(figsize=(12, 4), dpi=300)

        x = self.residue
        y = self.Get_TauC()
        self.SaveResults(x, y, 'npR_tauC_data.txt')

        plt.plot(x, y, marker = 'o')
        plt.xlabel('residue #')
        plt.xticks([a for a in range(650, 850, 10)])
        plt.xlim(650, 840)
        plt.ylabel('TauC')

        if filename == '':
            plt.show()
        else:
            plt.savefig(filename + 'tauC')
        return None
    ### ==================================================================== ###
    def SaveResults(self, x, y, filename):
        """ """
        with open(filename, 'w') as datafile:
            for a, b in zip(x, y):
                datafile.write(' '.join([str(a), str(b), '\n']))
        #
        return None
    ### ==================================================================== ###
    def get_resi(self):
        """ """
        return self.resi
    ### ==================================================================== ###
    def get_any(self, string):
        """ """
        data = []
        for i in range(len(self.resi)):
            data.append(self.result[self.resi[i]][string])
        return data
    ### ==================================================================== ###
    def get_average(self, string):
        """ """
        data = []
        for i in range(len(self.resi)):
            data.append((self.result[self.resi[i]][string+' - lower limit'] + self.result[self.resi[i]][string+' - uplter limit'])/ 2.0)
        return data
    ### ==================================================================== ###
    def Get_J_0870(self):
        """ """
        return self.get_any('J(0.870 wH)')
    ### ==================================================================== ###
    def Get_J_0921(self):
        """ """
        return self.get_any('J(0.921 wH)')
    ### ==================================================================== ###
    def Get_J_0955(self):
        """ """
        return self.get_any('J(0.955 wH)')
    ### ==================================================================== ###
    def Get_TauC(self):
        """ """
        return self.get_any('TauC')
    ### ==================================================================== ###
    def Get_JwN(self):
        """ """
        return self.get_average('J(wN)')
    ### ==================================================================== ###
    def Get_J0(self):
        """ """
        return self.get_average('J(0)')
    ### ==================================================================== ###
    residue = property(get_resi)
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###




class ReadDataFile():
    """
    It read datafiles and do simply modifications on data

    Usage:
    ======
        d = zbasic.ReadDataFile('data.txt',0,[2,5,3])
        d.delta = -0.5
        print d.x
        print d.y(2)
        for i in range(d.NumberOfY):
            print d.y(i)
    """
    ### ==================================================================== ###
    def __init__(self, filename, key, columns = [], skipline = 0):
        """ """
        self.__deltaX = 0.0
        # return variable
        self.__data = {}
        self.__keys = []
        # read the actual datafile
        try:
            datafile = open(filename,'r')
        except IOError:
            print 'Error opening',filename,'\n'
            exit()
        lines = datafile.readlines()
        datafile.close()
        # Extract the info from the file
        for i, line in enumerate(lines):
            # Do anything if it is not an empty line
            if (line != '') or (line != '\n') or ((i+1) < skipline):
                # Split the line into colomns
                col = line.split()
                # Only if the key colomn exists
                if len(col) >= key:
                    # initiate an empty array
                    self.__data[col[key]] = []
                    self.__keys.append(col[key])
                    # if everything is needed
                    if columns == []:
                        for c in col:
                            self.__data[col[key]].append(c)
                    # if limited colomns are needed
                    else:
                        for pos in columns:
                            # only if the position exists
                            if len(col) >= pos:
                                self.__data[col[key]].append(col[pos])
                            else:
                                self.__data[col[key]].append('-n/a-')
        return None
    ### ==================================================================== ###
    def y(self, number):
        """ """
        y = []
        for i in range(len(self.__keys)):
            if number <= len(self.__data[self.__keys[i]]):
                try:
                    y.append(float(self.__data[self.__keys[i]][number]))
                except ValueError:
                    y.append(self.__data[self.__keys[i]][number])
        return y
    ### ==================================================================== ###
    def yk(self, key, data = []):
        """ """
        y = []
        if key in self.__data.keys():
            if data == []:
                for i in range(len(self.__data[key])):
                    data.append(i)
            for d in data:
                try:
                    y.append(float(self.__data[key][d]))
                except ValueError:
                    y.append(self.__data[key][d])
        return y
    ### ==================================================================== ###
    def _getx(self):
        """ """
        xx = []
        for i in range(len(self.__keys)):
            try:
                xx.append(float(self.__keys[i]) + self.__deltaX)
            except ValueError:
                xx.append(self.__keys[i])
        return xx
    ### ==================================================================== ###
    def _getnumberofydata(self):
        """ """
        return len(self.__data[self.__keys[0]])
    ### ==================================================================== ###
    def _getdeltaX(self):
        """ """
        return self.__deltaX
    ### ==================================================================== ###
    def _setdeltaX(self, number):
        """ """
        self.__deltaX = float(number)
        return None
    ### ==================================================================== ###
    x         = property(_getx)
    NumberOfY = property(_getnumberofydata)
    delta     = property(fget = _getdeltaX, fset = _setdeltaX)
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

R1  = ReadDataFile('npR_R1_rates.txt', 0, [])
R2  = ReadDataFile('npR_R2_rates.txt', 0, [])
NOE = ReadDataFile('npR_hetNOE.txt', 0, [])


S = SpectralDensityFunction(R1.x, R1.y(1), R2.x, R2.y(1), NOE.x, NOE.y(1))
S.Plot('')
