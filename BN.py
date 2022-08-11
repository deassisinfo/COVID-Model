import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

random.seed(0)

components = ["Virus","ACE2","PKC","ANG_2_T1R","ANG_2","ADAM_17","SIL6R","TLR4","RIG1","NFKB","TNF","IRF3","STAT1","ISG","C_FLIP","INF_A_B","NRLP3","CASP1","FOXO3A","IFNR","BCL_2","tBID","Bax_Bak","CASP9","ROS","TNFR","FADD","Pyroptosis","IL1","IL1R","MLKL","Necroptosis","RIPK1_3","CASP8","Apoptosis"]
x = {k : 0 for k in components}

x["TNF"] = 1
x["ACE2"] = 1
#x[ROS] = 1

# find experiments that do this to confirm (grid search of genes array)

n_iter = 25
mat = np.zeros((n_iter + 1, len(x)))
# future: add k loop with average values, reset before each loop, consider a cumulative matrix (divide by k iter)
# consider simplification in the future, add negative feedback loops
for j in range(0, n_iter):
    r = random.sample(list(x), len(x))
    for i in r:
        if (i == "Virus"):
            x[i] = not x["ISG"] and x["Virus"] # or x["RIG1"]
            #x[i] = 1
            # x[i] = 0 # possibly make virus = 0
        if (i == "ACE2"):
            # PKC mediates ACE2 shedding from tubular cells
            x[i] = not x["Virus"] #x["PKC"] and not x["ADAM_17"]
            # x[i] = not x["Virus"]
            # not sure about if there is a relation since Virus just relies on ACE2 to enter cells, the presence of ACE2 promotes disease prog
        if (i == "PKC"):
            x[i] = x["ANG_2_T1R"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9000463/
        if (i == "ANG_2_T1R"):
            x[i] = x["ANG_2"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8193025/
            # Ang II acts on angiotensin type 1 (AT1") receptors and activates the NADPH-oxidase complex producing superoxide and promoting cell pro-oxidative and pro-inflammatory responses
        if (i == "ANG_2"):
            # ACE2 converts Ang II into Ang-(1–7")
            x[i] = not x["ACE2"]
        if (i == "ADAM_17"):
            # A2_t1r activates ADAM 17, which promotes ACE2 shedding (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185693/")
            x[i] = x["ANG_2_T1R"]
        if (i == "SIL6R"):
            x[i] = x["ADAM_17"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7224649/
        if (i == "TLR4"):
            x[i] = x["Virus"]
            # Spike glycoprotein, the major infective surface protein of SARS-CoV-2 has been found as a ligand for human TLR4
            # https://www.futuremedicine.com/doi/full/10.2217/fvl-2021-0249
        if (i == "RIG1"):
            x[i] = not x["ACE2"]
            # antiviral activity of RIG-1 may comprise inhibition of viral entry into the host cell by preventing the expression of its receptor, ACE2
            # https://www.news-medical.net/news/20210215/RIG-1-like-receptors-may-play-dominant-role-in-suppressing-SARS-CoV-2-infection.aspx
        if (i == "NFKB"):
            x[i] = x["TLR4"] or x["TNFR"] or x["IL1R"] #x["ANG_2_T1R"] or x["PKC"] or x["RIG1"] or
            # should be good, may need to find what exactly inhibits NFKB
            # common drug therapy is inhibiting NFKB (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648206/")
        if (i == "TNF"):
            x[i] = x["NFKB"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7114322/
        if (i == "IRF3"):
            x[i] = x["RIG1"] and not x["Virus"]
            # https://journals.asm.org/doi/10.1128/CMR.00299-20
            # SARS-CoV-2 membrane protein binds to importin karyopherin subunit alpha-6 (KPNA6") to inhibit interferon regulatory factor 3(IRF3") nuclear translocation
            # https://www.frontiersin.org/articles/10.3389/fcimb.2021.766922/full
        if (i == "STAT1"):
            x[i] = x["IFNR"] and not x["Virus"]
            # After the infection, STAT1 activity is inhibited by the SARS-CoV-2 proteins, NSP1, and ORF6
            # https://www.nature.com/articles/s41418-020-00633-7
        if (i == "ISG"):
            x[i] = x["STAT1"]
            # https://www.nature.com/articles/s41586-021-03234-7
        if (i == "C_FLIP"):
            x[i] = x["NFKB"] and not x["FOXO3A"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "INF_A_B"):
            x[i] = x["IRF3"] and not x["Virus"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7995242/
        if (i == "NRLP3"):
            x[i] = x["NFKB"]
            # https://www.nature.com/articles/ni.3772
        if (i == "CASP1"):
            x[i] = x["NRLP3"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6651423/
        if (i == "FOXO3A"):
            #x[i] = x["Virus"]
            x[i] = x["FOXO3A"]
            # no direct relation (drug targetting in place - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8187014/")
        if (i == "IFNR"):
            x[i] = x["INF_A_B"]
            # need to confirm more
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.606456/full
        if (i == "BCL_2"):
            x[i] = x["NFKB"]
            # https://www.nature.com/articles/1204926
        if (i == "tBID"):
            x[i] = x["CASP8"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882451/
        if (i == "Bax_Bak"):
            x[i] = x["BCL_2"] or x["tBID"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "CASP9"):
            x[i] = x["Bax_Bak"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "ROS"):
            x[i] = 0 #nothing happens - still need to figure out what exactly should be happening
        if (i == "TNFR"):
            x[i] = x["TNF"]
            # Do more research later
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.585880/full
        if (i == "FADD"):
            x[i] = x["TNFR"]
            # Do more research on this later
            # https://pubmed.ncbi.nlm.nih.gov/9430227/
        if (i == "Pyroptosis"):
            x[i] = x["CASP1"]
            # https://www.nature.com/articles/s41467-019-09753-2
        if (i == "IL1"):
            x[i] = x["MLKL"] or x["NFKB"]
            # Look into this more later
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == "IL1R"):
            x[i] = x["IL1"]
        if (i == "MLKL"):
            x[i] = x["RIPK1_3"] #and not (x["NRLP3"]")# or x["CASP1"]")
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
            # read more of this later
        if (i == "Necroptosis"):
            x[i] = x["MLKL"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == "RIPK1_3"):
            x[i] = (x["RIG1"] or x["TLR4"] or x["STAT1"]) and not x["CASP8"]
            # read more later
            # https://www.sciencedirect.com/science/article/pii/S1097276514008661
        if (i == "CASP8"):
            x[i] = (x["FADD"] or x["ROS"]) and not x["C_FLIP"] #or x["Virus"]
            # too many, will assume all is true, only certain drugs inhibit casp 8
        if (i == "Apoptosis"):
            x[i] = x["CASP8"]
            # https://pubmed.ncbi.nlm.nih.gov/10200555/
    mat[j,:]=list(x.values())
mat[n_iter, :] = np.average(mat[0:n_iter-1,:],axis=0)
print("Complete")
yticklabels = [str(x) for x in range(1,n_iter + 1)]
yticklabels.append("Average")
ax = sns.heatmap(mat, cmap="viridis", linewidths=.05, xticklabels=components, yticklabels=yticklabels)
ax.tick_params(axis='y', which='major', labelsize= 10)

print(mat[n_iter, :])

fig = plt.gcf()
fig.set_size_inches(10, 7, forward=True)
plt.tight_layout()
plt.show()