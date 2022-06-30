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
#x = [Virus,ACE2,PKC,ANG_2_T1R,ANG_2,ADAM_17,SIL6R,TLR4,RIG1,NFKB,TNF,IRF3,STAT1,ISG,C_FLIP,INF_A_B,NRLP3,CASP1,FOXO3A,IFNR,ProCASP8,DISC,BCL_2,tBID,Bax_Bak,CASP9,ROS,TNFR,FADD,Gasdermin,Pyroptosis,IL1,MLKL,Necroptosis,RIPK1_3,CASP8,ProCASP3_7,CASP3_7,Apoptosis]
#print(x)
#x = np.random.randint(0,2,len(x))

x = np.zeros(39)
x[Virus] = 1

#print(x)
x_begin = x.copy()

r = random.sample(range(len(x)), len(x))
n_iter = 25
mat = np.zeros((n_iter + 1, len(x)))
# future: add k loop with average values, reset before each loop, consider a cumulative matrix (divide by k iter)
# consider simplification in the future, add negative feedback loops
# got email that hipergator account will be deactivated
for j in range(0, n_iter):
    for i in r:
        #print(x[r[i]]);
        if (i == Virus):
            x[Virus] = not (x[STAT1] or x[RIG1]) and x[Virus]
            # x[Virus] = 0 # possibly make virus = 0 
        if (i == ACE2):
            # PKC mediates ACE2 shedding from tubular cells
            x[ACE2] = x[PKC]
            # x[ACE2] = not x[Virus] 
            # not sure about if there is a relation since Virus just relies on ACE2 to enter cells, the presence of ACE2 promotes disease prog
        if (i == PKC):
            x[PKC] = x[ANG_2_T1R]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9000463/
        if (i == ANG_2_T1R):
            x[ANG_2_T1R] = x[ANG_2] 
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8193025/
            # Ang II acts on angiotensin type 1 (AT1) receptors and activates the NADPH-oxidase complex producing superoxide and promoting cell pro-oxidative and pro-inflammatory responses
        if (i == ANG_2):
            # ACE2 converts Ang II into Ang-(1–7)
            x[ANG_2] = not x[ACE2]
        if (i == ADAM_17):
            # A2_t1r activates ADAM 17, which promotes ACE2 shedding (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185693/)
            x[ADAM_17] = x[ANG_2_T1R]
        if (i == SIL6R):
            x[SIL6R] = x[ADAM_17]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7224649/
        if (i == TLR4): 
            x[TLR4] = x[Virus]
            # Spike glycoprotein, the major infective surface protein of SARS-CoV-2 has been found as a ligand for human TLR4
            # https://www.futuremedicine.com/doi/full/10.2217/fvl-2021-0249
        if (i == RIG1):
            x[RIG1] = not x[ACE2]
            # antiviral activity of RIG-1 may comprise inhibition of viral entry into the host cell by preventing the expression of its receptor, ACE2
            # https://www.news-medical.net/news/20210215/RIG-1-like-receptors-may-play-dominant-role-in-suppressing-SARS-CoV-2-infection.aspx
        if (i == NFKB): 
            x[NFKB] = x[ANG_2_T1R] or x[PKC] or x[RIG1] or x[C_FLIP] or x[TNFR]
            # should be good, may need to find what exactly inhibits NFKB
            # common drug therapy is inhibiting NFKB (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648206/)
        if (i == TNF):
            x[TNF] = x[NFKB]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7114322/
        if (i == IRF3):
            x[IRF3] = x[RIG1] and not x[Virus]
            # https://journals.asm.org/doi/10.1128/CMR.00299-20
            # SARS-CoV-2 membrane protein binds to importin karyopherin subunit alpha-6 (KPNA6) to inhibit interferon regulatory factor 3(IRF3) nuclear translocation
            # https://www.frontiersin.org/articles/10.3389/fcimb.2021.766922/full
        if (i == STAT1):
            x[STAT1] = x[ISG] or x[IFNR] and not x[Virus]
            # After the infection, STAT1 activity is inhibited by the SARS-CoV-2 proteins, NSP1, and ORF6 
            # https://www.nature.com/articles/s41418-020-00633-7
        if (i == ISG):
            x[ISG] = x[Virus]
            # https://www.nature.com/articles/s41586-021-03234-7
        if (i == C_FLIP):
            x[C_FLIP] = x[NFKB] and not x[FOXO3A]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == INF_A_B):
            x[INF_A_B] = x[IRF3]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7995242/
        if (i == NRLP3):
            x[NRLP3] = x[NFKB]
            # https://www.nature.com/articles/ni.3772
        if (i == CASP1):
            x[CASP1] = x[NRLP3]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6651423/
        if (i == FOXO3A):
            #x[FOXO3A] = x[Virus]
            pass
            # no direct relation (drug targetting in place - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8187014/)
        if (i == IFNR):
            x[IFNR] = x[INF_A_B]
            # need to confirm more
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.606456/full
        if (i == ProCASP8):
            x[ProCASP8] = x[DISC]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8899608/
        if (i == DISC):
            x[DISC] = (x[ProCASP8] or x[FADD]) and not x[C_FLIP]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == BCL_2):
            x[BCL_2] = x[NFKB]
            # https://www.nature.com/articles/1204926
        if (i == tBID):
            x[tBID] = x[CASP8]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882451/
        if (i == Bax_Bak):
            x[Bax_Bak] = x[BCL_2] or x[tBID]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == CASP9):
            x[CASP9] = x[Bax_Bak]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == ROS):
            x[ROS] = 0 #nothing happens - still need to figure out what exactly should be happening
        if (i == TNFR):
            x[TNFR] = x[TNF]
            # Do more research later
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.585880/full
        if (i == FADD):
            x[FADD] = x[TNFR]
            # Do more research on this later
            # https://pubmed.ncbi.nlm.nih.gov/9430227/
        if (i == Gasdermin):
            x[Gasdermin] = x[CASP1]
            # Do more research on this later
            # https://www.cell.com/immunity/fulltext/S1074-7613(20)30237-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761320302375%3Fshowall%3Dtrue
        if (i == Pyroptosis):
            x[Pyroptosis] = x[Gasdermin]
            # https://www.nature.com/articles/s41467-019-09753-2
        if (i == IL1):
            x[IL1] = x[MLKL]
            # Look into this more later
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == MLKL):
            x[MLKL] = x[RIPK1_3] and not (x[NRLP3] or x[CASP1])
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
            # read more of this later
        if (i == Necroptosis):
            x[Necroptosis] = x[MLKL]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == RIPK1_3):
            x[RIPK1_3] = not x[CASP8] or x[RIG1]
            # read more later
            # https://www.sciencedirect.com/science/article/pii/S1097276514008661
        if (i == CASP8):
            x[CASP8] = x[Virus] or x[FADD] or x[ROS] or x[TNFR] or x[DISC]
            # too many, will assume all is true, only certain drugs inhibit casp 8
        if (i == ProCASP3_7):
            x[ProCASP3_7] = x[CASP8]
            # https://pubmed.ncbi.nlm.nih.gov/9765224/
        if (i == CASP3_7):
            x[CASP3_7] = x[ProCASP3_7]
            # obvious relationship based on naming
        if (i == Apoptosis):
            x[Apoptosis] = x[CASP3_7]
            # https://pubmed.ncbi.nlm.nih.gov/10200555/
    mat[j,:]=x
mat[n_iter, :] = np.average(mat,axis=0)
print("Complete")
yticklabels = [str(x) for x in range(1,n_iter + 1)]
yticklabels.append("Average")
xticklabels = ["Virus","ACE2","PKC","ANG_2_T1R","ANG_2","ADAM_17","SIL6R","TLR4","RIG1","NFKB","TNF","IRF3","STAT1","ISG","C_FLIP","INF_A_B","NRLP3","CASP1","FOXO3A","IFNR","ProCASP8","DISC","BCL_2","tBID","Bax_Bak","CASP9","ROS","TNFR","FADD","Gasdermin","Pyroptosis","IL1","MLKL","Necroptosis","RIPK1_3","CASP8","ProCASP3_7","CASP3_7","Apoptosis"]
ax = sns.heatmap(mat, cmap="viridis", linewidths=.05, xticklabels=xticklabels, yticklabels=yticklabels)
ax.tick_params(axis='y', which='major', labelsize= 10)

plt.show()
print("Nothing changed") if sum(np.subtract(x,x_begin)) == 0 else print(np.subtract(x,x_begin))