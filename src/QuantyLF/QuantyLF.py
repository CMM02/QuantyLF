#!/usr/bin/python
import numpy as np

from scipy.optimize import curve_fit, least_squares, nnls
from scipy.interpolate import interp1d
import scipy as sp
import scipy.signal as sig
from lmfit.models import GaussianModel,VoigtModel, ConstantModel, LinearModel, QuadraticModel, LorentzianModel
from lmfit import Model
import lmfit
import math
import sys
import os.path
import time
import subprocess

from multiprocessing import Process
import pickle

QuantyCommand = "/Users/phr542/Documents/Quanty/Quanty_macOS"
QuantyFile = "/Users/phr542/Documents/Packages/QuantyLF/src/QuantyLF/quanty/LF_3d_Td.lua"

EdgeJumpObject = "xxxx.pickle"

class QuantyLF:

    def __init__(self):
        pass

    ###############################################################################


    ####### Model the XAS edge jump and add it to the calculated output ###########

    def get_current_path(self):
        path = os.path.expanduser('~')
        return path

    #A function to read in the edge jump in XAS as interpolated from experiment
    def edge_jump(self, EdgeJumpObject):
        with open(str(EdgeJumpObject),'rb') as f:
            return pickle.load(f)

    def add_edge_jump(self, EdgeJumpObject,energy_scale,calcSpec):
        edge_jump_calc = self.edge_jump(EdgeJumpObject)(energy_scale)
        calc_edge = edge_jump_calc+calcSpec
        
        return calc_edge

    ###############################################################################
    #A function which shifts a spectrum in energy by amount par[0],
    #then interpolates it to the energy spacing of the second spectrum,
    #then scales it by factor par[1], and then returns the difference
    #between the spectra

    #With edge jump fitting
    def e_shift_res_edge_jump(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp, add edge jump and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        calcSpecNewJump = self.add_edge_jump(EdgeJumpObject,expSpec[:,0],calcSpecNew)
        return calcSpecNewJump - expSpec[:,1]

    #Similar as above but now without edge jump fitting
    def e_shift_res(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        return calcSpecNew - expSpec[:,1]


    #Similar to above, but now returns error in derivatives - better for peak matching?
    def e_shift_res_deriv(self, pars,calcSpec,expSpec):

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
    def quanty_res(self, pars,allPars,type):

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
        #subprocess.call([QuantyCommand, QuantyFile],stdout=subprocess.DEVNULL)
        subprocess.call([QuantyCommand, QuantyFile])#,stdout=None)

        
        #load spectra and experiment to compare
        calcSpec = np.loadtxt("XAS_Calc.dat")        
        
        #find the energy of largest peak in calculated and
        #experiment, to give a first rough shift for aligning energy
        calcPeak = calcSpec[calcSpec[:,1].argmax(),0]
        expPeak = self.expSpec[self.expSpec[:,1].argmax(),0]
        
        #parameters for the least squares function to fit the shift
        #and amplitude of the calculated XAS to compare to experiment
        dE = np.array([expPeak-calcPeak, 1]) #need to give it a close guess for energy shift
        amp = np.array([1])
        lowlim = np.array([-1e5,0])
        highlim = np.array([1e5,np.inf])
        
        #perform a non-linear least squares fit of the energy shift and amplitude of the 
        #calculated XAS, in order to compare it to experiment
        #Two versions are below - one using regular difference for least squares, and one using derivatives
        
        res = least_squares(self.e_shift_res,dE,bounds=(lowlim,highlim),max_nfev=200,args=(calcSpec,self.expSpec),verbose=0)#

        #Get the difference between calculated and experiment
        diff = res.fun
        diff_XAS = np.true_divide(diff,len(self.expSpec))
        
        #Write the XAS to file (this is on the calculated energy grid)
        calcSpec[:,0] = calcSpec[:,0] + res.x[0]
        calcSpec[:,1] = calcSpec[:,1] * res.x[1]
        
        #calcSpec[:,1] = addEdgeJump(EdgeJumpObject,calcSpec[:,0],calcSpec[:,1])
        
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
                    f.write(str(parVals[i]-res.x[0])) #shifted resonance energy according to XAS shift determined above
                    RIXSEner.append(parVals[i])
                else:
                    f.write(str(parVals[i])) #other parameter (non RIXS)
                f.write("\n")
            f.close()    
            
            #call Quanty to do the RIXS calculation with the current set of parameters
            subprocess.call([QuantyCommand, QuantyFile])#,stdout=subprocess.DEVNULL)

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

    
    def fit_pars(self, pars,type):
        params = lmfit.Parameters()
        for i in pars:
            if i[4] == 1:
                params.add(i[0],value=i[1],min=i[2],max=i[3])

    #QuantyRes(pars,type)
        global lsiter
        lsiter = 1

        minimizer = lmfit.Minimizer(self.quanty_res,params,fcn_args=(pars,type))
        res = minimizer.minimize(method='powell',params=params,options={'xtol': 1e-12})

        print("Final values: ")
        print(res.params)

    def fit(self):
        parList = []

        # Set up ion and oxidation state
        parList.append(["ion",25,25,25,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit
        parList.append(["oxy",3,3,3,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        parList.append(["Gamma1",0.4120250470353196,0.4,1,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit

        # Crystal field contribution in D4h symmetry
        parList.append(["tenDq",-1.4494920445039643,-2,0.1,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit
        parList.append(["tenDqF",0.6815489483432551,0.01,1.0,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit

        # Multiplet contribution
        # spin orbit coupling
        parList.append(["zeta_2p",1.0196625781428472,0.8,1.02,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit
        parList.append(["zeta_3d",0.8403012992370478,0.8,1.02,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit
        parList.append(["Xzeta_3d",1.9487449271777082,0.8,1.02,0]) #name, start val, low lim, high lim, 0/1  for nofit/fit

        # Slater integrals (Coulomb repulsion/exchange correlation)
        parList.append(["Fdd",0.9397729329705585,0.8,1.0,1]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        parList.append(["XFdd",0.8137253445941214,0.8,1.0,1]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        parList.append(["Fpd",0.8098173584848158,0.8,1.0,1]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        parList.append(["Gpd",0.8053014352519605,0.8,1.0,1]) #name, start val, low lim, high lim, 0/1 for nofit/fit


        # Ligand field contribution
        # on-site energies (usually drops out of the equation in crystal field theory)
        parList.append(["Udd",6.543685631427877,2.0,7.0,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        parList.append(["Upd_Udd",4.001467895225598,0.5,5.0,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit

        # Crystal field contribution of ligand site
        parList.append(["tenDqL",0.022975132006965073,0.01,1.0,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit

        # Charge transfer contribution
        parList.append(["Delta",4.040660314729548,1.0,5.0,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit

        # Hybridization
        parList.append(["VfScale",0.9775495515653867,0.8,1.0,0]) #name, start val, low lim, high lim, 0/1 for nofit/fit
        ## Define the incident RIXS energies 0 457.958 463.569
        for i in range(len(self.RIXS_energies)):
            parList.append(["RIXS",self.RIXS_energies[i],0,900,0])
        ## Want to fit XAS or RIXS (incl. XAS)0 457.96 463.569
        self.fit_pars(parList,"RIXS")

    def load_exp_xas(self, path):
        self.expSpec = np.loadtxt(path)

    def load_exp_rixs(self, path, RIXS_energies):
        self.RIXS_energies = RIXS_energies
        
        #load exp RIXS. For experimental, the first row are the resonance energies
        expRIXS = np.loadtxt(path)
        
        #trim the exp RIXS just to have the columns with resonant energies we are calculating
        indices = [0] #0th column is energy, which we will keep
        for i in range(len(RIXS_energies)):
            for j in range(len(expRIXS[0,:])):
                if(abs(RIXS_energies[i]-expRIXS[0,j]) < 0.1):
                    indices.append(j)
            
        self.expRIXS = expRIXS[1:,indices] #this removes the first row as well, which had the energy values (which we no longer need)
        