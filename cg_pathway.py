#! usr/bin/env python
# -*- coding: utf-8 -*-

# By Yixin Chen. Last edited: 16/01/2024

import numpy as np
from concentration_gradient import ConcentrationGradient
from matplotlib import pyplot as plt

class CoarseGrainedPathway(ConcentrationGradient):
    
    def __init__(self, shadowColors=None):
        
        self.N_rotModes = 0
        self.N_paramSets = 0
        self.rot_modes = []
        self.trans_modes = ["bind", "hydr", "release"]
        self.probs = dict()
        self.probs_sorted = dict()
        self.paramIndex_sorted = []
        
        if shadowColors is None:
            cmap1 = plt.cm.get_cmap("Set3")
            cmap2 = plt.cm.get_cmap("Dark2")
            self.shadowColors = [cmap1(i) for i in range(cmap1.N)]
            self.shadowColors += [cmap2(i) for i in range(cmap2.N)]


    def __readProbs(self, filePath, cT):

        log_c, rates, _, rotModes, transModes, dictRotModes = self.readFile(filePath)
        diff = abs(log_c["T"] - np.log10(cT))
        i = np.argmin(diff)
        print(log_c["T"][i])

        dictProbs = {trans: [] for trans in transModes}
        
        self.rot_modes = []
        for mode in rotModes:
            if mode == "80-120-200-240-320-360":
                continue
            self.rot_modes.append(mode)
        self.rot_modes.append("80-120-200-240-320-360")

        for trans in dictProbs:
            for mode in self.rot_modes:
                dictProbs[trans].append(dictRotModes[mode][trans][i]/(rates["k_cat"][i]/3.))
        
        return transModes, dictProbs
        
    
    def __getProbs(self, filePaths, cT):
        
        #self.rot_modes = []
        self.probs.clear()
        
        for fp in filePaths:
            transModes, dictProbs = self.__readProbs(fp, cT)

            for trans in transModes:
                if trans not in self.probs:
                    self.probs[trans] = []

                self.probs[trans].append(dictProbs[trans])
        
        self.N_rotModes = len(self.rot_modes)
        self.N_paramSets = len(filePaths)
        for trans in transModes:
            self.probs[trans] = np.array(self.probs[trans])
        
        '''
        for i in range(self.N_rotModes):
            plt.plot([0, 1], [i, i], color=self.shadowColors[i], lw = 2, label=self.rot_modes[i])
        plt.legend()
        plt.show()
        plt.close()
        '''

        return
    

    def __sortParamSets(self, probs, sortModes=None):
    
        paramGroup = {r: [] for r in range(self.N_rotModes)}
    
        for n in range(self.N_paramSets):
            r_max = np.argmax(probs[n])
            paramGroup[r_max].append(n)

        if sortModes is None:
            sortModes = range(self.N_rotModes)
    
        accumProb = np.zeros(self.N_paramSets)
        for r in sortModes:
            accumProb += probs[:, r]
            maxProb_with_index = list(zip(paramGroup[r], [accumProb[n] for n in paramGroup[r]]))
            maxProb_with_index.sort(key = lambda x: x[1], reverse = True)
        
            self.paramIndex_sorted += [x[0] for x in maxProb_with_index]
    
        for r in range(self.N_rotModes):
            if r not in sortModes:
                maxProb_with_index = list(zip(paramGroup[r], [probs[n][r] for n in paramGroup[r]]))
                maxProb_with_index.sort(key = lambda x: x[1], reverse = True)
        
                self.paramIndex_sorted += [x[0] for x in maxProb_with_index]
    
        return


    def __plotProbs(self, ax, probs, plotModes=None, plotThreshold=0.0):
    
        if plotModes is None:
            if plotThreshold == 0.0:
                plotModes = range(self.N_rotModes)
            else:
                plotModes = []
                for r in range(self.N_rotModes):
                    if max(probs[:, r]) > plotThreshold:
                        plotModes.append(r)
        
        accumProb = np.zeros(self.N_paramSets)
    
        for r in plotModes:
        
            ax.bar(np.arange(self.N_paramSets), probs[:, r], bottom = accumProb, label = self.rot_modes[r], color = self.shadowColors[r], width = 1.)

            accumProb += probs[:, r]
    
        ax.bar(np.arange(self.N_paramSets), np.ones(self.N_paramSets)-accumProb, bottom = accumProb, label = "Other rotation mode(s)", color = "darkred", width=1.)
        #ax.legend()
    
        return


    def __find_most_likely_pathway(self, outFileName, most_likely_ns = 3):
        
        with open("%s.pathway.csv"%outFileName, "w") as file:
            file.write("PARAM_INDEX")
            for trans in self.trans_modes:
                for i in range(most_likely_ns):
                    file.write(",%s-#%d_MODE,PROB"%(trans.upper(), i+1))
            file.write("\n")
            
            for n, par_index in enumerate(self.paramIndex_sorted):
                file.write("%d"%par_index)

                for trans in self.trans_modes:
                    rot_mode_probs_sorted = sorted([(i, self.probs_sorted[trans][n][i]) for i in range(self.N_rotModes)], key=lambda x: x[1], reverse=True)
                    for i in range(most_likely_ns):
                        file.write(",%s,%.6f"%(self.rot_modes[rot_mode_probs_sorted[i][0]],
                                               rot_mode_probs_sorted[i][1]))
                
                file.write("\n")
                
        return



    def __call__(self, filePaths, outFileName):
        
        
        with open("%s_cT=1e-7_old.pathway.csv"%outFileName, "r") as file:
            lines = file.readlines()
        for line in lines[1:]:
            line = line.strip().split(",")
            self.paramIndex_sorted.append(int(line[0]))
        

        #fig, ax = plt.subplots(3, 2, figsize = (6, 6), sharex = True, sharey = True)
        
        self.__getProbs(filePaths, 1e-7)
        '''
        #self.__sortParamSets(self.probs["bind"], sortModes = [0, 1])
        for a, trans in enumerate(self.trans_modes):
            self.probs_sorted[trans] = np.array([self.probs[trans][n] for n in self.paramIndex_sorted])
            self.__plotProbs(ax[a][0], self.probs_sorted[trans], plotThreshold = 0.)
        
        self.__find_most_likely_pathway(outFileName+"_cT=1e-7")
        
        self.__getProbs(filePaths, 1e-3)
        for a, trans in enumerate(self.trans_modes):
            self.probs_sorted[trans] = np.array([self.probs[trans][n] for n in self.paramIndex_sorted])
            self.__plotProbs(ax[a][1], self.probs_sorted[trans], plotThreshold = 0.)
            
        self.__find_most_likely_pathway(outFileName+"_cT=1e-3")
        
        for i in range(3):
            ax[i][0].set_ylim(0., 1.)
            ax[i][0].set_ylabel("Accu. probability")
        for i in range(2):
            ax[2][i].set_xlim(-0.5, 75.5)
            ax[2][i].set_xlabel("#ParameterSet") 
        
        
        for i in range(3):
          for j in range(2):
              ax[i][j].axvline(2.5, color="black", lw=1)
              ax[i][j].axvline(20.5, color="black", lw=1)
              ax[i][j].axvline(61.5, color="black", lw=1)
        
              #ax[i][j].legend()
        
        plt.tight_layout()
        plt.savefig("%s.png"%outFileName, dpi = 500, bbox_inches = "tight", transparent = True)
        plt.show()
        plt.close()
        '''
        
        for i in range(self.N_rotModes):
            plt.plot([0, 1], [i, i], color=self.shadowColors[i], lw = 2, label=self.rot_modes[i])
        plt.legend()
        plt.savefig("legend.png")
        plt.show()
        plt.close()
        
