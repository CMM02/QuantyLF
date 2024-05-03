#!/usr/bin/python
import numpy as np

from scipy.optimize import least_squares, nnls
from scipy.interpolate import interp1d
import scipy as sp
import scipy.signal as sig
import lmfit
import math
import sys
import os.path
import time
import subprocess
import platform

from multiprocessing import Process
import pickle

import matplotlib.pyplot as plt

QuantyFile = "/Users/phr542/Documents/Packages/QuantyLF/src/QuantyLF/quanty/LF_3d_Td.lua"


class QuantyLF:

    def __init__(self):
        self.edge_jump_interp = None
        self.quanty_command = {"default": "Quanty"}
        self.platform = platform.system()
        self.par_list = []
        self.par_file = 'ParVals.txt'
        self.file_par_dict = {}
        self.__read_par_file__()

    def __read_par_file__(self):
        # check if the file exists
        if not os.path.isfile(self.par_file):
            return

        with open(self.par_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.split(maxsplit=1)
                self.file_par_dict[line[0]] = line[1].strip()


    ####### Model the XAS edge jump and add it to the calculated output ###########

    def __edge_jump__(self, e, pos, jump, slope):
        return jump/np.pi * (np.arctan((e-pos)*slope) + np.pi/2)

    def __get_quanty_command__(self):
        if self.platform in self.quanty_command.keys():
            return self.quanty_command[self.platform]
        else:
            return self.quanty_command["default"]


    """
    Configure the edge jump for the XAS calculation. The edge jump is modelled as a arctan function.

    Parameters:
    edge_jumps: List of edge jumps. Each element in the list is a list of 3 elements: [position, jump, slope]
    x_range: The range of x values to model the edge jump over
    y_offset: The y offset of the edge jump
    display: If True, the edge jump will be displayed along with the experimental XASs
    """
    def config_edge_jump(self, edge_jumps, x_range, y_offset=0, display=False):
        x_range = np.linspace(x_range[0], x_range[1], math.ceil((max(x_range)-min(x_range))/0.1))
        edge_jump_y = [y_offset] * len(x_range)
        for i in range(len(edge_jumps)):
            cur_jump = edge_jumps[i]
            edge_jump_y += self.__edge_jump__(x_range, cur_jump[0], cur_jump[1], cur_jump[2])
        
        self.edge_jump_interp = interp1d(x_range, edge_jump_y, kind='cubic', fill_value="extrapolate")

        if display:
            if self.expXAS is None:
                print("To display edge jump with experimental XAS, load the experimental XAS first")
            else:
                plt.plot(self.expXAS[:,0], self.expXAS[:,1])
            plt.plot(x_range, edge_jump_y)
            plt.show()
        

    #A function to read in the edge jump in XAS as interpolated from experiment
    # def edge_jump(self, EdgeJumpObject):
    #     with open(str(EdgeJumpObject),'rb') as f:
    #         return pickle.load(f)

    def __add_edge_jump__(self, energy_scale,calcSpec):
        if self.edge_jump_interp is None:
            raise ValueError("Configure edge jump first, before running calculation with edge jump")
        edge_jump_calc = self.edge_jump_interp(energy_scale)
        calc_edge = edge_jump_calc+calcSpec
        
        return calc_edge

    ###############################################################################
    #A function which shifts a spectrum in energy by amount par[0],
    #then interpolates it to the energy spacing of the second spectrum,
    #then scales it by factor par[1], and then returns the difference
    #between the spectra

    #With edge jump fitting
    def __e_shift_res_edge_jump__(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp, add edge jump and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        calcSpecNewJump = self.__add_edge_jump__(expSpec[:,0],calcSpecNew)
        return calcSpecNewJump - expSpec[:,1]

    #Similar as above but now without edge jump fitting
    def __e_shift_res__(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        return calcSpecNew - expSpec[:,1]


    #Similar to above, but now returns error in derivatives - better for peak matching?
    def __e_shift_res_deriv__(self, pars,calcSpec,expSpec):

    #shift calc by par, then interp and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        atanerr = np.zeros_like(calcSpecNew)
        
        for i in range(len(atanerr)):
            if(i>0):
                atanerr[i] = np.arctan2(calcSpecNew[i]-calcSpecNew[i-1],expSpec[i,0]-expSpec[i-1,0]) - np.arctan2(expSpec[i,1]-expSpec[i-1,1],expSpec[i,0]-expSpec[i-1,0])        
        return atanerr

    
    #A function which runs a Quanty calculation for a given set of parameters,
    #then compares the resulting spectrum against an experimental spectrum
    #read from file, and returns the difference between the spectra
    def __quanty_res__(self, pars,allPars,type):

        global lsiter
        
        #create a list of the names of the pars being fitted and the other pars which are not
        parNames = [x[0] for x in allPars if x[4] <= 0]
        parVals = [x[1] for x in allPars if x[4] <= 0]

        for k,v in pars.items():
            parNames.append(v.name)
            parVals.append(v.value)

        #Write the current parameters to file, so Quanty can read
        #them and do the calculation. Note for this part we make
        #sure not to write any RIXS energy parameters, because the
        #RIXS is fitted later. Here just XAS calculation is done.
        f = open("ParVals.txt","w")
        for i in range(len(parNames)):
            if(parNames[i] != "RIXS"):
                f.write(parNames[i])
                f.write(" ")
                f.write(str(parVals[i]))
                f.write("\n")   
        f.write("XAS 0\n")
        f.close()
        

        #Print out current values of parameters, so they can be tracked
        print("=====================")
        for i in range(len(parNames)):
            print(parNames[i],parVals[i])
        print("=====================")


        #run Quanty - it will read the parameters from file, calculate 
        #the XAS, and write the new spectrum/spectra to file XAS_Calc.dat
        # subprocess.call([self.__get_quanty_command__(), QuantyFile],stdout=subprocess.DEVNULL)
        subprocess.call([self.__get_quanty_command__(), QuantyFile])#,stdout=None)

        
        #load spectra and experiment to compare
        calcSpec = np.loadtxt("XAS_Calc.dat")        
        
        #find the energy of largest peak in calculated and
        #experiment, to give a first rough shift for aligning energy
        calcPeak = calcSpec[calcSpec[:,1].argmax(),0]
        expPeak = self.expXAS[self.expXAS[:,1].argmax(),0]
        
        #parameters for the least squares function to fit the shift
        #and amplitude of the calculated XAS to compare to experiment
        dE = np.array([expPeak-calcPeak, 1]) #need to give it a close guess for energy shift
        amp = np.array([1])
        lowlim = np.array([-1e5,0])
        highlim = np.array([1e5,np.inf])
        
        #perform a non-linear least squares fit of the energy shift and amplitude of the 
        #calculated XAS, in order to compare it to experiment
        #Two versions are below - one using regular difference for least squares, and one using derivatives
        
        res_fn = self.__e_shift_res_edge_jump__ if self.edge_jump_interp is not None else self.__e_shift_res__

        res = least_squares(res_fn,dE,bounds=(lowlim,highlim),max_nfev=200,args=(calcSpec,self.expXAS),verbose=0)#

        #Get the difference between calculated and experiment
        diff = res.fun
        diff_XAS = np.true_divide(diff,len(self.expXAS))
        
        #Write the XAS to file (this is on the calculated energy grid)
        calcSpec[:,0] = calcSpec[:,0] + res.x[0]
        calcSpec[:,1] = calcSpec[:,1] * res.x[1]
        
        if self.edge_jump_interp is not None:
            calcSpec[:,1] = self.__add_edge_jump__(calcSpec[:,0],calcSpec[:,1])
        
        np.savetxt("XAS_Fit.dat",calcSpec)


        #If requested by the user, now do a fit of parameters based
        #on RIXS spectra
        if(type == "RIXS"):

            #rewrite the par file now with RIXS energies included, shifted appropriately
            f = open("ParVals.txt","w")
            RIXSEner = []
            for i in range(len(parNames)):
                f.write(parNames[i])
                f.write(" ")
                if(parNames[i]=="RIXS"):
                    print(parNames[i],parVals[i])
                    f.write(str(parVals[i]-res.x[0])) #shifted resonance energy according to XAS shift determined above
                    RIXSEner.append(parVals[i])
                else:
                    f.write(str(parVals[i])) #other parameter (non RIXS)      
                f.write("\n")
            f.close()    
            
            #call Quanty to do the RIXS calculation with the current set of parameters
            subprocess.call([self.__get_quanty_command__(), QuantyFile])#,stdout=subprocess.DEVNULL)

            #load the calculated RIXS spectra
            calcRIXS = np.loadtxt("RIXS_Calc.dat")
            
            
            #fit scaling of the RIXS spectra to get best agreement
            #with experiment (linear least squares fit, single scaling for all)
            #to do this, concatenate all the RIXS spectra
            calcRIXS2 = np.copy(self.expRIXS)
            for i in range(len(calcRIXS2[0,:])-1):
                calcRIXS2[:,i+1] = np.interp(self.expRIXS[:,0],calcRIXS[:,0],calcRIXS[:,i+1])
                calcRIXSCat = np.copy(calcRIXS2[:,1])
                expRIXSCat = np.copy(self.expRIXS[:,1])
            for i in range(len(calcRIXS[0,:])-2):
                calcRIXSCat = np.hstack((calcRIXSCat,calcRIXS2[:,i+2]))
                expRIXSCat = np.hstack((expRIXSCat,self.expRIXS[:,i+2]))
            
            #do the linear least squares fit
            amp,res = nnls(np.array(calcRIXSCat[:,None]),expRIXSCat)

            print(amp, res)
            
            #Apply the scaling factor to the calculated RIXS
            #(Both the interpolated and the original calculated)
            for i in range(len(calcRIXS2[0,:])-1):
                calcRIXS2[:,i+1] = amp[0] * calcRIXS2[:,i+1]
                calcRIXS[:,i+1] = amp[0] * calcRIXS[:,i+1]

            #save RIXS, return concatenated differences
            np.savetxt("RIXS_Fit.dat",calcRIXS)
            np.savetxt("RIXS_Exp_Trimmed.dat",self.expRIXS)
            diff = calcRIXS2[:,1] - self.expRIXS[:,1]
            for i in range(len(calcRIXS2[0,:])-2):
                diff = np.hstack((diff,calcRIXS2[:,i+2] - self.expRIXS[:,i+2]))
            print("Chi2: ",np.dot(diff,diff))
            sys.stdout.flush()
            
            diff_RIXS = np.true_divide(abs(calcRIXS2[:,1] - self.expRIXS[:,1]),len(self.expRIXS))
            counter = 1
            for i in range(len(calcRIXS2[0,:])-2):
                diff_RIXS = np.hstack((diff_RIXS,np.true_divide(abs(calcRIXS2[:,i+2] - self.expRIXS[:,i+2]),len(self.expRIXS))))
                counter += 1
            
            diff_RIXS2 = np.dot(diff_RIXS,diff_RIXS)/counter
            diff_XAS2 = np.dot(diff_XAS,diff_XAS)
            
            if lsiter == 1:
                global initialXAS
                global initialRIXS
                
                initialXAS = diff_XAS2
                initialRIXS = diff_RIXS2
                                
            diff_weighted = 0.5*(diff_XAS2/initialXAS) + 0.5*(diff_RIXS2/initialRIXS)
                                    
            print("Diff XAS",diff_XAS2)
            print("Diff RIXS",diff_RIXS2)
            print("Diff weighted",diff_weighted)
            lsiter +=1
                                    
            return diff_weighted

            
        #XAS fitting, so return XAS error
        print("Chi2: ",np.dot(diff,diff))
        sys.stdout.flush()
        
        return  diff #should be able to return something like res.error? 

    
    def __fit_pars__(self,type):
        params = lmfit.Parameters()
        for i in self.par_list:
            if i[4] == 1:
                params.add(i[0],value=i[1],min=i[2],max=i[3])

    #QuantyRes(pars,type)
        global lsiter
        lsiter = 1

        minimizer = lmfit.Minimizer(self.__quanty_res__,params,fcn_args=(self.par_list,type))
        res = minimizer.minimize(method='powell',params=params,options={'xtol': 1e-12})

        print("Final values: ")
        print(res.params)


    """
    Add a parameter to the model

    ----------------
    Parameters:
    name: Name of the parameter
    init_val: Initial value of the parameter
    interv: List of two values, lower and upper bounds of the parameter
    from_file: If True (default True), the parameter is read from a file (file value overrides init_val)
    """
    def add_par(self, name, init_val, interv = None, from_file = True):
        if interv and (interv[1] - interv[0]) < 0:
            raise ValueError("Upper bound of parameter should be greater than lower bound")
        
        low = 0
        high = 0
        if interv:
            low = interv[0]
            high = interv[1]
            if name in self.file_par_dict.keys():
                init_val = float(self.file_par_dict[name]) if from_file else init_val
        else:
            if name in self.file_par_dict.keys():
                init_val = self.file_par_dict[name] if from_file else init_val
        self.par_list.append([name, init_val, low, high, 1 if interv else 0])


    """
    Fit the parameters of the model to the experimental data

    ----------------
    Parameters:
    mode: "XAS" or "RIXS". If "XAS", only the XAS data is fitted. If "RIXS", both XAS and RIXS data is fitted
    """
    def fit(self, type):        
        self.__fit_pars__(type)

    """
    Load the experimental XAS data from a file (no header, two columns: energy, intensity)

    ----------------
    Parameters:
    path: Path to the file containing the experimental XAS data
    """
    def load_exp_xas(self, path):
        self.expXAS = np.loadtxt(path)


    """
    Load the experimental RIXS data from a file (no header, first row: resonance energies, rest of the rows: energy, intensity)

    ----------------
    Parameters:
    path: Path to the file containing the experimental RIXS data
    RIXS_energies: List of resonance energies for which the RIXS data is available
    """
    def load_exp_rixs(self, path, RIXS_energies):        
        for RIXS_energy in RIXS_energies:
            self.add_par("RIXS",RIXS_energy,from_file=False)
        
        #load exp RIXS. For experimental, the first row are the resonance energies
        expRIXS = np.loadtxt(path)
        
        #trim the exp RIXS just to have the columns with resonant energies we are calculating
        indices = [0] #0th column is energy, which we will keep
        for i in range(len(RIXS_energies)):
            for j in range(len(expRIXS[0,:])):
                if(abs(RIXS_energies[i]-expRIXS[0,j]) < 0.1):
                    indices.append(j)
            
        self.expRIXS = expRIXS[1:,indices] #this removes the first row as well, which had the energy values (which we no longer need)


    """
    Set the path to the Quanty executable

    ----------------
    Parameters:
    platform: Platform for which the path is being set
    command: Path to the Quanty executable
    """
    def set_quanty_command(self, command, platform='default'):
        self.quanty_command[platform] = command


    """
    Set custom path to the parameter file (default: ParVals.txt)

    ----------------
    Parameters:
    par_file: Path to the parameter file
    """
    def set_par_file(self, par_file):
        self.par_file = par_file
        