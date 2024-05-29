#! usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright by Yixin Chen (yixin.chen@mpinat.mpg.de).
# Last edited: 04/01/2023

from itertools import product
from variant import *

def angularDiff(x, y):
    
    d = abs(x - y)

    if d > 180:
        return 360-d
    else:
        return d


class MarkovModel:
    
    def __init__(self, 
                 confs="ohc",
                 v_index=1,
                 activeConf=[2],#list of integers, representing the conformations that are catalytically active (allow ATP hydrolysis/synthesis)
                 name="F1-ATPase",
                 angle_cat=80,
                 jointConfTrans=False):
        # Gamma orientation: with fold of symmetry = 3, (80, 200, 320) & (120, 240, 360) are two symmetric groups.
        #                    Use (120, 240, 360) rather than (0, 120, 240) to make gamma at 80 and 120 point to the
        #                    same beta, and name this beta beta_A.
        #                    Increased angle corresponds to counterclockwise rotation.
        # Binding state: 2 for empty, 0 for ATP-bound, 1 for ADP-bound
        # Conformations: 0 for o (open), 1 for h (half-closed), 2 for c (closed), 3 for c* (alternative closed)
        
        self.name = name
        self.activeConf = activeConf
        self.angle_cat = angle_cat
        self.__jointConfTrans = jointConfTrans
        
        self.__variant = Variant(confs, v_index, activeConf)
        self.version = self.__variant.version
        self.N_conf = len(self.__variant.confStates)
        self.outFileName = "%s_%s"%(self.name, self.version)
        
        self.__dictConfName = {0: "o", 1: "h", 2: "c", 3: "a"}
        
        self.stepwiseAngle = [360, self.angle_cat, 120, self.angle_cat+120]
        self.betaConfig = {each: set() for each in [self.angle_cat, 120, self.angle_cat+120, 240, self.angle_cat+240, 360]}
        
        self.__bindStates = [0, 1, 2]
        self.__bindStatesTag = {0: "ATP-bound", 1: "ADP-bound", 2: "empty"}
        self.__stepwiseBindingState = {2: [(0, "+ATP"), (1, "+ADP")],
                                       0: [(2, "-ATP"), (1, "HYD")],
                                       1: [(2, "-ADP"), (0, "SYN")]}

        
        self.N = -9999
        self.asyState = []# each row: [gamma, beta_A conf, beta_B conf, beta_C, conf, beta_A bind, beta_B bind, beta_C bind]
        self.state = dict()
        
        self.accSubspace = dict()
        # Items in accSubspace: key as state index i, value as a list of tuples, including all states j that can be reached from state i
        # Each tuple is (state index j, transition type, optional info)
        
        #=============================================================================================
        # Attributes for studying rotation modes (you don't need to worry about this part now, Niklas)
        
        self.__rot_modes_asy     = [(angle_cat,), 
                                    (angle_cat, 120), 
                                    (angle_cat, 120, angle_cat+120,),
                                    (angle_cat, 120, angle_cat+120, 240),
                                    (angle_cat, 120, angle_cat+120, 240, angle_cat+240),
                                    (120,),
                                    (120, angle_cat+120,),
                                    (120, angle_cat+120, 240),
                                    (120, angle_cat+120, 240, angle_cat+240),
                                    (120, angle_cat+120, 240, angle_cat+240, 360),
                                    (angle_cat, 120, angle_cat+120, 240, angle_cat+240, 360)]
        self.__rot_modes_asy_set = [set(each) for each in self.__rot_modes_asy]
        
        self.__rot_modes_all = dict()
        self.__rot_modes_sym = []
        self.__rot_modes_sym_set = []
        for rotMode in self.__rot_modes_asy[:-1]:
            self.__rot_modes_all[rotMode] = []
            
            for s in [1, 2]:
                
                rotMode1 = []
                for angle in rotMode:
                    angle_sym = (angle + s*120) % 360
                    if angle_sym == 0:
                        angle_sym = 360
                    rotMode1.append(angle_sym)
                    
                self.__rot_modes_all[rotMode].append(tuple(rotMode1))
                self.__rot_modes_sym.append(tuple(rotMode1))
                self.__rot_modes_sym_set.append(set(rotMode1))
                
        self.__rot_modes_all[self.__rot_modes_asy[-1]] = []
        
        self.__configs = dict()
        self.config2rot_mode = dict()
        
        self.trans_modes = ["bind", "release", "hydr"]
        self.trans2rot_mode = {trans: dict() for trans in self.trans_modes}
        self.N_rot_modes = -9999
        self.rot_modes = []
        
        #====================================================================================================

    
    #Functions for generating Markov states====================================================
    
    def generateAsyMarkovStates(self):
        
        bindCombs = list(product(self.__bindStates, repeat = 3))
        
        confCombs_cat = list(product(self.__variant.gbCoupling_cat[0],self.__variant.gbCoupling_cat[1], self.__variant.gbCoupling_cat[2]))
        confCombs_wait = list(product(self.__variant.gbCoupling_wait[0],self.__variant.gbCoupling_wait[1], self.__variant.gbCoupling_wait[2]))
        
        for each in confCombs_cat:
            
            self.betaConfig[self.angle_cat].add(each)
            
            for each1 in bindCombs:
                
                state = [self.angle_cat] + list(each) + list(each1)
                self.asyState.append(state)
        
        
        for each in confCombs_wait:
            
            self.betaConfig[120].add(each)
            
            for each1 in bindCombs:
                
                state = [120] + list(each) + list(each1)
                self.asyState.append(state)
        
        self.N = len(self.asyState)
        for i in range(self.N):
            self.state[tuple(self.asyState[i])] = i
        
        return


    def addSymmetricMarkovStates(self):
        
        # No return, but the input dictMarkovState is updated
        
        file = open(self.outFileName+"_state.markov", "w")
        
        for i in range(self.N):
        
            state = self.asyState[i]
            symState1 = (state[0]+120, state[3], state[1], state[2], state[6], state[4], state[5])
            symState2 = (state[0]+240, state[2], state[3], state[1], state[5], state[6], state[4])
            
            self.betaConfig[state[0]+120].add(symState1[1:4])
            self.betaConfig[state[0]+240].add(symState2[1:4])
            
            #120 degree rotation
            self.state[symState1] = i + self.N
            #240 degree rotation
            self.state[symState2] = i + 2*self.N
            
            file.write("%4d%4d%15s%15s%15s%15s%15s%15s"%(i+1, state[0], 
                       self.__variant.cStatesTag[state[1]], self.__variant.cStatesTag[state[2]], self.__variant.cStatesTag[state[3]],
                       self.__bindStatesTag[state[4]], self.__bindStatesTag[state[5]], self.__bindStatesTag[state[6]]))
            
            file.write("%8d%4d%15s%15s%15s%15s%15s%15s"%(i+self.N+1, symState1[0], 
                            self.__variant.cStatesTag[symState1[1]], self.__variant.cStatesTag[symState1[2]], self.__variant.cStatesTag[symState1[3]],
                            self.__bindStatesTag[symState1[4]], self.__bindStatesTag[symState1[5]], self.__bindStatesTag[symState1[6]]))
            
            file.write("%8d%4d%15s%15s%15s%15s%15s%15s\n"%(i+2*self.N+1, symState2[0], 
                            self.__variant.cStatesTag[symState2[1]], self.__variant.cStatesTag[symState2[2]], self.__variant.cStatesTag[symState2[3]],
                            self.__bindStatesTag[symState2[4]], self.__bindStatesTag[symState2[5]], self.__bindStatesTag[symState2[6]]))
        
        file.close()
        
        return
    
    #==========================================================================================
    
    
    # Functions for generating accessible subspace==============================================
    # Will then be used to determine transition rate matrix
    
    def searchAccessibleSubspaceByBindingStates(self, i, state):
        
        for b in [4, 5, 6]:# change binding state
            
            for each in self.__stepwiseBindingState[state[b]]:
                
                if ((each[1] in ["HYD", "SYN"]) and (state[b-3] in self.__variant.inactiveConf)) == False:#change with model
                
                    newState = state[0:b] + [each[0]] +state[(b+1):]
                    
                    self.accSubspace[i].add((self.state[tuple(newState)], each[1], state[b-3]))
        
        return


    def searchAccessibleSubspaceByBetaConfTrans(self, i, state):
        
        if self.__jointConfTrans == True:

          for comb in list(product(self.__variant.stepwiseConfTrans[state[1]]+[state[1]],
                                           self.__variant.stepwiseConfTrans[state[2]]+[state[2]],
                                           self.__variant.stepwiseConfTrans[state[3]]+[state[3]])):
            
              if comb != (state[1], state[2], state[3]):
                
                  newState = (state[0], comb[0], comb[1], comb[2], state[4], state[5], state[6])
                
                  if newState in self.state:
                      self.accSubspace[i].add((self.state[newState], "CT"))
        
        else:

            for beta in range(1, 4):
          
                c = state[beta]

                for each in self.__variant.stepwiseConfTrans[c]:
                    newState = [] + state
                    newState[beta] = each
                    newState = tuple(newState)

                    if newState in self.state:
                        self.accSubspace[i].add((self.state[newState], "CT"))
      
        return


    def searchAccessibleSubspaceByGammaRot(self, i, state):
        
        for dir in [1, -1]:
            
            newGamma = self.stepwiseAngle[self.stepwiseAngle.index(state[0])+dir]
            onlyRotState = [newGamma] + state[1:]
            
            if tuple(onlyRotState) in self.state:
                
                self.accSubspace[i].add((self.state[tuple(onlyRotState)], "RT"))
                    
        return
    #===========================================================================================
    

    def construct(self):

        # Generate Markov states===============================================================
        self.generateAsyMarkovStates()
        self.addSymmetricMarkovStates()
        #======================================================================================
        
        # Generate accessible subspace for each (asymmetric) Markov state======================
        
        for i in range(self.N):
            self.accSubspace[i] = set()
            
            state = self.asyState[i]
            
            self.searchAccessibleSubspaceByBindingStates(i, state)# change binding state
            self.searchAccessibleSubspaceByBetaConfTrans(i, state)# change beta conformation (no gamma rotation)
            self.searchAccessibleSubspaceByGammaRot(i, state)     # change gamma (no beta conformation change)
        #======================================================================================
        
        self.betaConfig["all"] = list(self.betaConfig[self.angle_cat] | self.betaConfig[120])
        
        return
    

    #Functions for finding rotation modes======================================================

    def findConfig(self):

        for s in self.state:

            if s[1:4] not in self.__configs:
                if s[0] <= 120:
                    self.__configs[s[1:4]] = [s[0]]
                
            else:
                if s[0] not in self.__configs[s[1:4]]:
                    self.__configs[s[1:4]].append(s[0])
        
        return
        

    def groupConfigByRotMode(self):

        for config in self.__configs:
            angles = set(self.__configs[config])
            if angles not in self.__rot_modes_asy_set:
                #print(config, angles)
                if angles not in self.__rot_modes_sym_set:
                    angles = [self.__configs[config][0]]
                    if len(self.__configs[config]) > 1:
                        #print("yes")
                        for angle_in in angles:
                            for angle in self.__configs[config][1:]:
                                #print(angle_in, angle)
                                if angularDiff(angle, angle_in) < 120 and angle not in angles:
                                    angles.append(angle)
                    angles = set(angles)
                    if angles not in self.__rot_modes_asy_set:
                        continue
                else:
                    continue
            rotMode = self.__rot_modes_asy[self.__rot_modes_asy_set.index(angles)]
            
            if rotMode not in self.config2rot_mode:
                self.config2rot_mode[rotMode] = []
                
            self.config2rot_mode[rotMode].append(config)
        
        return


    def addSymRotMode(self):

        for rotMode in self.__rot_modes_all:
            if rotMode in self.config2rot_mode:
                
                for s, rotMode_sym in enumerate(self.__rot_modes_all[rotMode]):
                    self.config2rot_mode[rotMode_sym] = []
                    
                    for config in self.config2rot_mode[rotMode]:
                        self.config2rot_mode[rotMode_sym].append((config[(-1-s)%3], config[(-s)%3], config[(1-s)%3]))
        
        with open("%s.config.csv"%self.outFileName, "w") as file:
            for config in self.__configs:
                file.write("%s-%s-%s"%(self.__dictConfName[config[0]], self.__dictConfName[config[1]], self.__dictConfName[config[2]]))
            
                for angle in self.__configs[config]:
                    file.write(",%d"%angle)

                file.write("\n")
        
            for rotMode in self.config2rot_mode:
                file.write("-".join(list(map(lambda x: str(x), rotMode))))
                for config in self.config2rot_mode[rotMode]:
                    file.write(",%s-%s-%s"%(self.__dictConfName[config[0]], self.__dictConfName[config[1]], self.__dictConfName[config[2]]))
                file.write("\n")
            
        return


    def groupStateByRotMode(self):
        
        self.state2rot_mode = {mode: [] for mode in self.config2rot_mode}
        
        for s in self.state:
        
            for mode in self.config2rot_mode:
            
                if (s[0] in mode) and (s[1:4] in self.config2rot_mode[mode]):
                    self.state2rot_mode[mode].append(self.state[s])
                    break
        
        count = 0
        for mode in self.state2rot_mode:
            count += len(self.state2rot_mode[mode])
        if count != 3*self.N:
            print(3*self.N, count)
            print("Error in grouping states according to rotation mode!")

        return


    def groupTransitionByRotMode(self):
        
        allState = sorted(self.state.items(), key = lambda x: x[1])
        
        for mode in self.state2rot_mode:
            
            for trans in self.trans2rot_mode:
                self.trans2rot_mode[trans][mode] = []
            
            for i in self.state2rot_mode[mode]:
                
                if allState[i][0][4] == 2:#beta_1 is empty
                    
                    j = self.state[allState[i][0][:4]+(0,)+allState[i][0][5:]]
                    self.trans2rot_mode["bind"][mode].append((i, j))
                    
                elif allState[i][0][4] == 1:#beta_1 is ADP bound
                    
                    j = self.state[allState[i][0][:4]+(2,)+allState[i][0][5:]]
                    self.trans2rot_mode["release"][mode].append((i, j))
                
                elif allState[i][0][4] == 0:#beta_1 is ATP bound
                    
                    if allState[i][0][1] in self.activeConf:
                        j = self.state[allState[i][0][:4]+(1,)+allState[i][0][5:]]
                        self.trans2rot_mode["hydr"][mode].append((i, j))
        
        return
    
    
    def findRotModes(self):

        self.findConfig()
        self.groupConfigByRotMode()
        self.addSymRotMode()
        
        self.rot_modes = [mode for mode in self.__rot_modes_asy+self.__rot_modes_sym if mode in self.config2rot_mode]
        self.N_rot_modes = len(self.rot_modes)
        
        self.groupStateByRotMode()
        self.groupTransitionByRotMode()
        
        return
    #======================================================================================

