##! usr/bin/env python
# -*- coding: utf-8 -*-

# By Yixin Chen. Last edited: 08/01/2024

import numpy as np
from scipy.optimize import curve_fit
from concentration_gradient import ConcentrationGradient

class BindingModel:
    
    @staticmethod
    def f1(x, kd1, n1):
        return n1*x/(x+kd1)
    
    @staticmethod
    def f2(x, kd1, kd2, n1, n2):
        return n1*x/(x+kd1)+n2*x/(x+kd2)

    @staticmethod
    def f2a(x, kd1, kd2):
        return 2.0*x/(x+kd1)+x/(x+kd2)
        
    @staticmethod
    def f3(x, kd1, kd2, kd3, n1, n2, n3):
        return n1*x/(x+kd1)+n2*x/(x+kd2)+n3*x/(x+kd3)
        
    @staticmethod
    def f3a(x, kd1, kd2, kd3):
        return x/(x+kd1)+x/(x+kd2)+x/(x+kd3)
    



class TitrationCurve(ConcentrationGradient):

    def __init__(self, expDataDir = "./ExpDataLib"):
        self.expDataDir = expDataDir
    
    
    def __fit_titration_curve(self, conc, occ):
        
        func = BindingModel.f3a
        
        try:
            popt, pcov = curve_fit(func, conc, occ, bounds = (0.0, 1.0), p0 = [1e-12, 1e-12, 4.7e-7])
            y_fit = func(conc, *popt)
            y_aver = np.average(occ)
            r_sq = 1 - np.dot(y_fit - occ, y_fit - occ)/np.dot(y_aver - occ, y_aver - occ)
            
            Kd = [(popt[i], np.sqrt(np.diag(pcov))[i]) for i in range(3)]
            Kd.sort(key = lambda x: x[0])
            
            print("R^2=%.4f"%r_sq)
            
            for i in range(3):
                print("K_d%d=%.1e +- %.1e"%(i+1, Kd[i][0], Kd[i][1]))
            
            return [conc, func(conc, *popt)]
        
        except RuntimeError:
        
            print("Curve fit failed")
            return 0
    
    
    def __calculateTheoreticalOccupancy(self, log_conc, which, 
                                        Kd1 = {"T": 16e-9, "D": 41e-9},
                                        Kd2 = {"T": 1.5e-6, "D": 6e-6},
                                        Kd3 = {"T": 29e-6, "D": 42e-6}):
        
        if which in ["T", "D"]:
            conc = 10**np.array(log_conc)
            return conc/(conc+Kd1[which]) + conc/(conc+Kd2[which]) + conc/(conc+Kd3[which])
        
        else:
            return np.zeros(len(log_conc))-1.0
    
    
    def __readExpTitCurve(self, which):
        
        log_conc, occ = [], []
        
        if which == "mix":
            tag = which
        else:
            tag = "A%sP"%which

        with open("%s/F1-ATPase_%s.bind.csv"%(self.expDataDir, tag), "r") as file:
            lines = file.readlines()
        for line in lines[1:]:
            line = line.split(",")
            log_conc.append(float(line[0]))
            occ.append(float(line[1]))
        
        return np.array(log_conc), np.array(occ)
    
    
    def __readModelTitCurve(self, which, filePath):
        
        log_c, _, occs, _, _, _ = self.readFile(filePath)
        
        return log_c[which], occs
    
    
    def __plotTitrationCurve(self, ax, model_tit, exp_tit):
        
        scat_exp = ax.scatter(10.0**exp_tit[0], exp_tit[1], label = "experimental data", marker = "o", color = "grey", s = 120)
        scat_model = ax.scatter(10.0**model_tit[0], model_tit[1]["tot"], label = "model prediction", marker = "s", edgecolor = "darkred", linewidth = 2, color = "white", s = 60)
        
        ax.scatter(10.0**model_tit[0], model_tit[1]["ATP"], label = "occ. by ATP", marker = "^", edgecolor = "darkblue", color = "white", s = 60, linewidth = 2)
        ax.scatter(10.0**model_tit[0], model_tit[1]["ADP"], label = "occ. by ADP", marker = "D", edgecolor = "darkblue", color = "white", s = 60, linewidth = 2)
        
        ax.set_ylabel("Occupancy", fontsize = 12)
        ax.set_ylim(0, 3.2)
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels([0, 1, 2, 3], fontsize = 10)
        
        ax.grid(axis = "y", ls = "--")
        
        '''
        ax.set_xlabel("[nucleotide] ($\\mathrm{\\mu}$M)", fontsize = 20)
        ax.set_xscale("log")
        ax.set_xlim(1e-8, 10**(-2.2))
        ax.set_xticks(10.0**np.arange(-8, -2.2, 1))
        ax.set_xticklabels([0.01, 0.1, 1, 10, 100, 1000], fontsize = 18)
        '''
        
        return scat_exp, scat_model
    
    
    def __call__(self, ax, modelFilePath, which):
        
        model_tit = self.__readModelTitCurve(which, modelFilePath)
        exp_tit = self.__readExpTitCurve(which)
        
        scat_exp, scat_model = self.__plotTitrationCurve(ax, model_tit, exp_tit)
        
        curve_exp_fit = ax.plot(10.0**model_tit[0], self.__calculateTheoreticalOccupancy(model_tit[0], which), color = "grey", ls = ":", lw = 2)#, label = "theoretical")
        
        res = self.__fit_titration_curve(10.0**model_tit[0], model_tit[1])
        if res != 0:
            curve_model_fit = ax.plot(res[0], res[1], ls = ":", color = "darkred", lw = 2)

        return scat_exp, scat_model, curve_exp_fit, curve_model_fit
