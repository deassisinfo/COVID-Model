# -*- coding: utf-8 -*-

import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
random.seed(0)

# Virus = 0
# ACE2 = 0
# PKC = 0
# ANG_2_T1R = 0
# ANG_2 = 0
# ADAM_17 = 0
# SIL6R = 0
# TLR4 = 0
# RIG1 = 0
# NFKB = 0
# TNF = 0
# IRF3 = 0
# STAT1 = 0
# ISG = 0
# C_FLIP = 0
# INF_A_B = 0
# NRLP3 = 0
# CASP1 = 0
# FOXO3A = 0
# IFNR = 0
# ProCASP8 = 0
# DISC = 0
# BCL_2 = 0
# tBID = 0
# Bax_Bak = 0
# CASP9 = 0
# ROS = 0
# TNFR = 0
# FADD = 0
# Gasdermin = 0
# Pyroptosis = 0
# IL1 = 0
# MLKL = 0
# Necroptosis = 0
# RIPK1_3 = 0
# CASP8 = 0
# ProCASP3_7 = 0
# CASP3_7 = 0
# Apoptosis = 0

# values represent indices, not boolean (x contain booleans at index positions)


# find frequencies of cell death, see what you think what should be happening


Virus,ACE2,PKC,ANG_2_T1R,ANG_2,ADAM_17,SIL6R,TLR4,RIG1,NFKB,TNF,IRF3,STAT1,ISG,C_FLIP,INF_A_B,NRLP3,CASP1,FOXO3A,IFNR,ProCASP8,DISC,BCL_2,tBID,Bax_Bak,CASP9,ROS,TNFR,FADD,Gasdermin,Pyroptosis,IL1,MLKL,Necroptosis,RIPK1_3,CASP8,ProCASP3_7,CASP3_7,Apoptosis = range(0,39)
x = [Virus,ACE2,PKC,ANG_2_T1R,ANG_2,ADAM_17,SIL6R,TLR4,RIG1,NFKB,TNF,IRF3,STAT1,ISG,C_FLIP,INF_A_B,NRLP3,CASP1,FOXO3A,IFNR,ProCASP8,DISC,BCL_2,tBID,Bax_Bak,CASP9,ROS,TNFR,FADD,Gasdermin,Pyroptosis,IL1,MLKL,Necroptosis,RIPK1_3,CASP8,ProCASP3_7,CASP3_7,Apoptosis]
#print(x)
#x = np.random.randint(0,2,len(x))
x = np.zeros(len(x))
x[Virus] = 1
#print(x)
x_begin = x.copy()

r = random.sample(range(len(x)), len(x))
n_iter = 100
mat = np.zeros((n_iter, len(x)))
for j in range(0, n_iter):
    for i in r:
        #print(x[r[i]]);
        if (i == Virus):
            x[Virus] = not (x[STAT1] or x[RIG1]) and x[Virus]
        if (i == ACE2):
            x[ACE2] = x[Virus]
        if (i == PKC):
            x[PKC] = x[ANG_2_T1R]
        if (i == ANG_2_T1R):
            x[ANG_2_T1R] = x[ANG_2]
        if (i == ANG_2):
            x[ANG_2] = not x[ACE2]
        if (i == ADAM_17):
            x[ADAM_17] = x[ANG_2_T1R] 
        if (i == SIL6R):
            x[SIL6R] = x[ADAM_17]
        if (i == TLR4): 
            x[TLR4] = x[Virus]
        if (i == RIG1): 
            x[RIG1] = x[Virus] ## check online again
        if (i == NFKB): 
            x[NFKB] = x[ANG_2_T1R] or x[PKC] or x[RIG1] or x[C_FLIP] or x[TNFR]
        if (i == TNF):
            x[TNF] = x[NFKB]
        if (i == IRF3):
            x[IRF3] = x[RIG1] and not x[Virus]
        if (i == STAT1):
            x[STAT1] = x[ISG] or x[IFNR]
        if (i == ISG):
            x[ISG] = x[Virus]
        if (i == C_FLIP):
            x[C_FLIP] = x[NFKB] and not x[FOXO3A]
        if (i == INF_A_B):
            x[INF_A_B] = x[IRF3]
        if (i == NRLP3):
            x[NRLP3] = x[NFKB]
        if (i == CASP1):
            x[CASP1] = x[NRLP3]
        if (i == FOXO3A):
            x[FOXO3A] = x[Virus]
        if (i == IFNR):
            x[IFNR] = x[INF_A_B]
        if (i == ProCASP8):
            x[ProCASP8] = x[DISC]
        if (i == DISC):
            x[DISC] = (x[ProCASP8] or x[FADD]) and not x[C_FLIP]
        if (i == BCL_2):
            x[BCL_2] = x[NFKB]
        if (i == tBID):
            x[tBID] = x[CASP8]
        if (i == Bax_Bak):
            x[Bax_Bak] = x[BCL_2] or x[tBID]
        if (i == CASP9):
            x[CASP9] = x[Bax_Bak]
        if (i == ROS):
            x[ROS] = 0 #nothing happens
        if (i == TNFR):
            x[TNFR] = x[TNF]
        if (i == FADD):
            x[FADD] = x[TNFR]
        if (i == Gasdermin):
            x[Gasdermin] = x[CASP1]
        if (i == Pyroptosis):
            x[Pyroptosis] = x[Gasdermin]
        if (i == IL1):
            x[IL1] = x[MLKL]
        if (i == MLKL):
            x[MLKL] = x[RIPK1_3] and not (x[NRLP3] or x[CASP1])
        if (i == Necroptosis):
            x[Necroptosis] = x[MLKL]
        if (i == RIPK1_3):
            x[RIPK1_3] = not x[CASP8] or x[RIG1]
        if (i == CASP8):
            x[CASP8] = x[Virus] or x[FADD] or x[ROS] or x[TNFR] or x[DISC]
        if (i == ProCASP3_7):
            x[ProCASP3_7] = x[CASP8]
        if (i == CASP3_7):
            x[CASP3_7] = x[ProCASP3_7]
        if (i == Apoptosis):
            x[Apoptosis] = x[CASP3_7]
    mat[j,:]=x
print("Complete")
sns.heatmap(mat)
plt.show()
print("Nothing changed") if sum(np.subtract(x,x_begin)) == 0 else print(np.subtract(x,x_begin))