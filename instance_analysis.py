#! usr/bin/env python
# -*- coding: utf-8 -*-

# By Yixin Chen. Last edited: 04/01/2024

import numpy as np
from instance import Instance


class InstanceAnalysis:

    def __init__(self, model, vec_para, workDir = ""):
        
        self.model = model
        self.vec_para = vec_para
        self.workDir = workDir
        self.instance = Instance(model)
        
        self.__bindStateAbbrev = {0: "T", 1: "D", 2: "E"}
        self.__confStateAbbrev = {0: "o", 1: "h", 2: "c", 3: "c*"}


    def __analyzeSteadyState(self, cT, cD, summaryFile=None):
        
        self.instance.set_trans_rate_matrix(self.vec_para, cT, cD)
        
        self.instance.calculate_steady_state_properties()
        self.instance.write_st_distr("%s/SteadyState_Distr_(%.1e,%.1e).csv"%(self.workDir, cT, cD))
        
        #self.instance.analyzeSteadyState_betaConfig("%s/SteadyState_BetaConfig_(%.1e,%.1e).csv"%(self.workDir, cT, cD))
        #self.instance.analyzeSteadyState_NucExchange("%s/SteadyState_Flux_(%.1e,%.1e).csv"%(self.workDir, cT, cD),
        #                                       "%s/SteadyState_NucExchange_(%.1e,%.1e).csv"%(self.workDir, cT, cD))
        self.instance.analyzeSteadyState_rotMode()

        print("Steady state ([ATP] = %.2e, [ADP] = %.2e)"%(cT, cD))
        print("k_cat = %.2f, k_rot = %.2f"%(self.instance.k_cat, self.instance.k_rot))
        
        if summaryFile is not None:
            
            summaryFile.write("%.2f,%.2f,%.6e,%.6e,%.6e,%.6e,%.6e"%(\
                              np.log10(cT), np.log10(cD), 
                              self.instance.k_cat, self.instance.k_rot,
                              self.instance.occ_tot, self.instance.occ_ATP, self.instance.occ_ADP))
            
            for n in range(self.model.N_rot_modes):
                summaryFile.write(",%.6e"%(self.instance.pop_rotModes[n]))
            
            for trans in self.model.trans_modes:
                for n in range(self.model.N_rot_modes):
                    summaryFile.write(",%.6e"%(self.instance.dictFlux_rotModes[trans][n]))
                    
            summaryFile.write("\n")
            
        return
    

    
    def __writeTitle(self, file):

        file.write("LOG([ATP]),LOG([ADP]),K_CAT,K_ROT,OCC_TOT,OCC_ATP,OCC_ADP")
        for mode in self.model.rot_modes:
            file.write(",%s"%("-".join([str(angle) for angle in mode])))
            
        for trans in self.model.trans_modes:
            for mode in self.model.rot_modes:
                file.write(",%s:%s"%(trans, "-".join([str(angle) for angle in mode])))
        
        file.write("\n")
      
        return 0


    def scan_conc_gradient(self, gradient = "ATP", concRange = {"T": np.logspace(-9, 0, 15), "D": [1e-9, 1e-6]}):
        
        if gradient in ["T", "D"]:
        
            c1_label = gradient
            
            if c1_label == "T":
                c2_label = "D"
            
            elif c1_label == "D":
                c2_label = "T"
            
            conc = {"T": 0.0, "D": 0.0}
            
            for c2 in concRange[c2_label]:
                
                conc[c2_label] = c2
                
                with open("%s/SteadyStates_c%s=%.2e.csv"%(self.workDir, c2_label, c2),"w") as file:
                    self.__writeTitle(file)
                    
                    for c1 in concRange[c1_label]:
                        
                        conc[c1_label] = c1
                        
                        self.__analyzeSteadyState(cT = conc["T"], cD = conc["D"], summaryFile = file)
        
        else:
                
            with open("%s/SteadyStates_mix.csv"%self.workDir,"w") as file:
                self.__writeTitle(file)

                for c in concRange["T"]:
                    self.__analyzeSteadyState(cT = c, cD = c, summaryFile = file)
        
        return 0