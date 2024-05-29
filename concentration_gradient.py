#! usr/bin/env python
# -*- coding: utf-8 -*-

# By Yixin Chen. Last edited: 08/01/2024

import numpy as np
from scipy.optimize import curve_fit

def convertDict(myDict):
    return {each: np.array(myDict[each]) for each in myDict}


class ConcentrationGradient:

    def __init__(self):
        self.turnover_exp = [np.array([-6, -5, -4, -3]), np.array([1.951, 16.00, 57.14, 76.92])]
    
    
    def readFile(self, filePath):
        
        with open(filePath, "r") as file:
            lines = file.readlines()
        
        log_c = {"T": [], "D": []}
        rates = {"k_cat": [], "k_rot": []}
        occs = {"tot": [], "ATP": [], "ADP": []}
        dict_rotModes = dict()
        
        line = lines[0].strip().split(",")[7:]
        N_rotModes = len(line)//4
        rotModes = []
        
        for each in line[:N_rotModes]:
            rotModes.append(each)
            dict_rotModes[each] = {"pop": []}
        
        transModes = []
        for each in line[N_rotModes::N_rotModes]:
            trans = each.split(":")[0]
            transModes.append(trans)
            for mode in dict_rotModes:
                dict_rotModes[mode][trans] = []
        #print(rotModes, transModes)
        
        for line in lines[1:]:
            
            line = np.array((line.strip()).split(","), dtype = float)
            log_c["T"].append(line[0])
            log_c["D"].append(line[1])
            rates["k_cat"].append(line[2])
            rates["k_rot"].append(line[3])
            occs["tot"].append(line[4])
            occs["ATP"].append(line[5])
            occs["ADP"].append(line[6])
            
            for i in range(N_rotModes):
                mode = rotModes[i]
                dict_rotModes[mode]["pop"].append(line[7+i])
                
                for j in range(1, 4):
                    dict_rotModes[mode][transModes[j-1]].append(line[7+i+j*N_rotModes])
        
        for mode in dict_rotModes:
            for each in dict_rotModes[mode]:
                dict_rotModes[mode][each] = np.array(dict_rotModes[mode][each])

        return convertDict(log_c), convertDict(rates), convertDict(occs),\
               rotModes, transModes, dict_rotModes
    
    
    def plotRates(self, ax, filePaths, which):
        
        ax.scatter(10.0**self.turnover_exp[0], self.turnover_exp[1], color = "grey", s = 120, marker = "o")
        
        for fp in filePaths:
        
            log_c, rates, _, _, _, _ = self.readFile(fp)
            
            concs = 10.00**log_c[which]
            
            ax.scatter(concs, rates["k_cat"], edgecolor = "darkred", s = 60, marker = "s", color = "white", linewidth = 2)
            ax.scatter(concs, 3*rates["k_rot"], edgecolor = "darkblue", s = 60, marker = "s", color = "white", linewidth = 2)
            
            res = self.__fit_rates(concs, rates["k_cat"])
            if res != 0:
                curve_model_fit = ax.plot(res[0], res[1], ls = ":", color = "darkred", lw = 2)
        
            res = self.__fit_rates(concs, 3*rates["k_rot"])
            if res != 0:
                curve_model_fit = ax.plot(res[0], res[1], ls = ":", color = "darkblue", lw = 2)
        
            res = self.__fit_rates(10.0**self.turnover_exp[0], self.turnover_exp[1])
            if res != 0:
                curve_model_fit = ax.plot(res[0], res[1], ls = ":", color = "grey", lw = 2)
        
        
        ax.set_ylabel("Rate (s$^{-1}$)", fontsize = 12)
        ax.set_yticks(np.linspace(0, 100, 6))
        ax.set_yticklabels(["%.0f"%y for y in np.linspace(0, 100, 6)], fontsize = 10)
        ax.set_ylim(-10, 100)
        ax.grid(axis = "y", ls = "--")
        
        ax.set_xlabel("[ATP] ($\\mathrm{\\mu}$M)", fontsize = 12)
        ax.set_xscale("log")
        ax.set_xlim(1e-8, 10**(-2.2))
        ax.set_xticks(10.0**np.arange(-8, -2.2, 1))
        ax.set_xticklabels([0.01, 0.1, 1, 10, 100, 1000], fontsize = 10)
        
        return 0
    
    
    def plotRotModesPopulation(self, ax, filePath, which):
        
        log_c, _, _, _, _, dict_rotModes = self.readFile(filePath)
        concs = 10.0**log_c[which]

        pop_cat = dict_rotModes["80"]["pop"] + dict_rotModes["200"]["pop"] + dict_rotModes["320"]["pop"]
        pop_wait = dict_rotModes["80-120"]["pop"] + dict_rotModes["200-240"]["pop"] + dict_rotModes["320-360"]["pop"]
        
        pop_irregular = 1. - pop_cat - pop_wait
        
        #ax.scatter(concs, pop_wait, marker = "s", edgecolor = "darkred", linewidth = 2, color = "white", s = 60)
        #ax.scatter(concs, pop_cat, marker = "s", edgecolor = "darkblue", linewidth = 2, color = "white", s = 60)
        #ax.scatter(concs, pop_irregular, marker = "s", edgecolor = "darkgray", linewidth = 2, color = "white", s = 60)
        
        ax.plot(concs, pop_wait, color = "tab:blue", alpha=0.6)
        ax.plot(concs, pop_cat, color = "tab:orange", alpha=0.6)
        ax.plot(concs, pop_irregular, color = "darkgrey", alpha=0.6)
        
        ax.set_xscale("log")
        ax.set_xlim(1e-8, 10**(-2.2))
        ax.set_xticks(10.0**np.arange(-8, -2.2, 1))
        ax.set_xticklabels([0.01, 0.1, 1, 10, 100, 1000], fontsize = 10)
        ax.set_xlabel("[ATP] ($\\mathrm{\\mu}$M)", fontsize = 12)
        
        
        ax.set_ylim(-0.1, 1.1)
        ax.set_yticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticklabels(["%.1f"%x for x in np.linspace(0.0, 1.0, 6)], fontsize = 10)
        ax.set_ylabel("Population", fontsize = 12)
        
        return 0
    
    
    @staticmethod
    def __MichaelisMenten(x, v_max, K_M):
        return v_max*x/(x+K_M)
    
    
    def __fit_rates(self, conc, rate):
        
        func = self.__MichaelisMenten
        
        try:
            popt, pcov = curve_fit(func, conc, rate, p0 = [80.0, 1e-5], bounds = (0.0, 1e3))
            y_fit = func(conc, *popt)
            y_aver = np.average(rate)
            r_sq = 1 - np.dot(y_fit - rate, y_fit - rate)/np.dot(y_aver - rate, y_aver - rate)
            
            params_MM = [(popt[i], np.sqrt(np.diag(pcov))[i]) for i in range(2)]
            
            print("R^2=%.4f"%r_sq)
            
            print("v_max=%.2f +- %.2f"%(params_MM[0][0], params_MM[0][1]))
            print("K_M=%.1e +- %.1e"%(params_MM[1][0], params_MM[1][1]))
            
            return [np.logspace(-9, 0, 51), func(np.logspace(-9, 0, 51), *popt)]
        
        except RuntimeError:
        
            print("Curve fit failed")
            return 0
