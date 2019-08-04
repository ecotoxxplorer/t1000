# Copyright (C) 2017, Xia Lab @ McGill,
# McGill University.
# All rights reserved.
# 
# (JUST TENTATIVE LICENSE POLICY !!!!!!!)
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the lab, McGill nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL the lab, McGill BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# Note: This program complies with any terms appears in the licenses of the used
# software programs.
#
# Name          : mainhelper.py
# Creation      : Nov 28, 2017
# Subject       : Holds functions needed to start running the main program.
# Author        : Othman Soufan othman.soufan@mcgill.ca
#
# usage           : import to main.py
# python_version  :2.7.13

from genemining.modelbuilders import *
from resvisualizer.statplots import *

from sklearn.neighbors import KNeighborsClassifier

from sklearn.model_selection import cross_val_predict

from sklearn import metrics
from sklearn.metrics import confusion_matrix

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

from sklearn.naive_bayes import GaussianNB

import os.path

def evaluatepreds(X, y, resfile, clf):
    
    if(clf == "QDA"):
        mdl = QuadraticDiscriminantAnalysis()
    elif(clf == "DT"):
        mdl = DecisionTreeClassifier()
    elif(clf == "KNN"):
        mdl = KNeighborsClassifier(n_neighbors=5)
    elif(clf == "NBC"):
        mdl = GaussianNB()
    elif(clf == "RF"):
        mdl = RandomForestClassifier(n_estimators=1, random_state=0, n_jobs=-1)
        
    np.random.seed(1)
    predicted = cross_val_predict(mdl.fit(X, y), X, y, cv=5)
     
    res = metrics.confusion_matrix(y, predicted)
    tp = res[0][0]
    tn = res[1][1]
    fp = res[1][0]
    fn = res[0][1]

    sens = tp*1.0/(tp+fn)
    spec = tn*1.0/(tn+fp)
    prec = tp*1.0/(tp+fp)

    scores_txt = str(metrics.recall_score(y, predicted))+"\t"+str(spec)+"\t"+str(metrics.precision_score(y, predicted))+"\t"+str(np.sqrt(sens*spec))
    scores_txt = scores_txt+"\t"+str(metrics.f1_score(y, predicted))+"\t"+str(metrics.fbeta_score(y, predicted, 0.5))+"\n"
    resfile.write(scores_txt)
    
def runRegressionAnal(dfile):
    # Read data
    print("Started reading data files for regression analysis.")
    path = dfile
    data = np.transpose(np.array(pd.read_csv(path+"fullmat", sep=" ", header=-1)))
    varnames = data[0,:]
    X = data[1:,:].astype(np.float)
    # Remove cases with assigned targets as NaN
    y = np.array(pd.read_csv(path+"y", header=-1)).astype(np.float)[:,0]
    dlvl = np.array(pd.read_csv(path+"dose_level", header=-1))
 
    #X = X[:, selfeaidx]
    X = X[~np.isnan(y), :]
    dlvl = dlvl[~np.isnan(y)]
    y = y[~np.isnan(y)]
    
    # Run experiments
    print("Started experiments for training models and activity prediction.")
     
    mdltypes = ["RT", "LR", "KNN", "BayesR"]
    preproctype = ["PCA", "Normal"]
     
    topk = 1000
    for m in mdltypes:
        for p in preproctype:
            print "Processing model: ", m, " with preprocessing type ",p
            runexps(X, y, dlvl, varnames, topk, m, p)
    
    print("Done regression analysis...")
    
def runGeneSelection(dfile, feasz, method):
    # Read data
    path = dfile
    data = np.transpose(np.array(pd.read_csv(path+"fullmat", sep=" ", header=-1)))
    varnames = data[0,:]
    X = data[1:,:].astype(np.float)
    # Remove cases with assigned targets as NaN
    y = np.array(pd.read_csv(path+"y", header=-1)).astype(np.float)[:,0]
    dlvl = np.array(pd.read_csv(path+"dose_level", header=-1))
    
    datafname = "/home/soufanom/tmp/X_tr"
    dosefname = "/home/soufanom/tmp/dose_level_tr"
    np.savetxt(datafname, np.transpose(np.round(X, 6)), fmt="%s")
    np.savetxt(dosefname, dlvl, fmt="%s")
      
    cmd = "Rscript /home/soufanom/git/gene-prioritization/geneprioritization/R/FeaSelMethods.R "+datafname+" "+dosefname+" "+str(feasz)+" "+method #" Limma"
    selidx = check_output(cmd.split())
    selidx = np.array(map(int, selidx.split("\n")[0].split(",")))-1
      
    ofile = open("/home/soufanom/Downloads/tmp_genes.txt", "w")
    for v in varnames[selidx]:
        ofile.write(v+"\n")
    ofile.close()
    print("Done") 
   
def plotSelFeaPCA(dfile, feafile, imgname, title, resfile, clf):
    # Read data
    path = dfile
    data = np.transpose(np.array(pd.read_csv(path+"fullmat", sep=" ", header=-1)))
    varnames = data[0,:]
    X = data[1:,:].astype(np.float)
    # Remove cases with assigned targets as NaN
    y = np.array(pd.read_csv(path+"y", header=-1)).astype(np.float)[:,0]
    dlvl = np.array(pd.read_csv(path+"dose_level", header=-1))
    
    #Map input indices
    #selfeaidx = np.array(map(int, [lines.strip() for lines in open("data/Gene_sets/biofeasel_idx.txt")]))
    if(feafile != ""):
        selfeaidx = featureidx(varnames, feafile)
        selfeaidx = np.unique(selfeaidx)
        if(len(selfeaidx) < 2000):
            selfeaidx = selfeaidx[0:999]
        print(len(selfeaidx))
        varnames = varnames[selfeaidx]
        X = X[:, selfeaidx]
    else:
        print(X.shape)
    
    X = X[~np.isnan(y), :]
    dlvl = dlvl[~np.isnan(y)]
    y = y[~np.isnan(y)]
   
    X_tmp, imp_vals, selfea =  impmissings(X)
    X_tmp = mapmissings(X, imp_vals, selfea)
    min_thresh = 96
    max_thresh = 104
    
    if("in_vivo" in dfile):
        y_tmp = np.zeros(len(y))
        y_tmp[y == 1] = 2
        y_tmp[y == 2] = 1
        y = y_tmp
    else:
        cond = np.array((y < min_thresh) | (y > max_thresh))
        y[cond] = 1
        y[~cond] = 2

    #pca = PCA(n_components=3)
    #X_tmp = pca.fit_transform(X_tmp)
    
    
    if(clf != ""):
        evaluatepreds(X_tmp, y, resfile, clf)
         
    #X_tmp, dum, y, dlvl = preprocess(X_tmp, X_tmp, y, dlvl, "RU")
 
 
    clfpcaplot(X_tmp, y, imgname, title)
    
def runClassificationRandomAnal(path):
    
    classifiers = ["RF"]

    directory = "Classification/Round2/Random/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    resfile = open(directory+"/prediction_scores_RF.txt", "w")

    #rat_genes = np.asarray([lines.strip() for lines in open('data/Gene_sets/rat_genes.txt')])

    #np.random.seed(1)
    for e in range(0, 1000):
        exp_name = "Random_"+str(e+1)
        #idx = np.random.randint(low=0, high=len(rat_genes), size=1000)	
        #genes = rat_genes[idx]
        feafile = "Classification/Random/"+"GeneLists/"+exp_name+".txt"
        #np.savetxt(feafile, genes, fmt="%s")
    
        for clf in classifiers:
	    exp_name2 = exp_name+"_"+clf       
            resfile.write(exp_name2+"\t") 
            imgname = 'Classification/Figures-pca/Rat_liver_random_ivivo_v3.tiff'
            title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (Random)"
            plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
    resfile.close()


def runClassificationAnal(path):
    
    classifiers = ["NBC", "KNN", "QDA", "DT", "RF"]
    
    for clf in classifiers:
        directory = "Classification/Round2/"+clf
        
        if not os.path.exists(directory):
            os.makedirs(directory)
    
        resfile = open(directory+"/prediction_scores_thousand.txt", "w")
        
        resfile.write("T1000\t")
        imgname = 'Classification/Figures-pca/Rat_liver_T1000_ivivo_v3.tiff'
        title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (T1000)"
        feafile = 'data/Gene_sets/ToxGPrio_1000chip_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        #resfile.write("All Genes\t") 
        #imgname = 'Classification/Figures-pca/Rat_liver_all_ivivo_v3.tiff'
        #title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (All Genes)"
        #feafile = 'data/Gene_sets/ToxGPrio_1000chip_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        #plotSelFeaPCA(path, "", imgname, title, resfile, clf)
        
        resfile.write("Limma\t") 
        imgname = 'Classification/Figures-pca/Rat_liver_limma_ivivo_v3.tiff'
        title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (Limma)"
        feafile = 'data/Gene_sets/limma_genes_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        resfile.write("CD\t")
        imgname = 'Classification/Figures-pca/Rat_liver_cd_ivivo_v3.tiff'
        title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (CD)"
        feafile = 'data/Gene_sets/CD_genes_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        #resfile.write("Random\t") 
        #imgname = 'Classification/Figures-pca/Rat_liver_random_ivivo_v3.tiff'
        #title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (Random)"
        #feafile = 'data/Gene_sets/random_genes_rat_r3.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        #plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        resfile.write("L1000\t")  
        imgname = 'Classification/Figures-pca/Rat_liver_l1000_ivivo_v3.tiff'
        title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (L1000)"
        feafile = 'data/Gene_sets/L1000_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        resfile.write("S1500\t")  
        imgname = 'Classification/Figures-pca/Rat_liver_s1500_ivivo_v3.tiff'
        title = "Rat $\it{in}$ $\it{vivo}$ Liver Expression Space (S1500)"
        feafile = 'data/Gene_sets/S1500_rat.txt'#_ecotoxchip/ecotoxchip_genes_ML_ranked.txt'
        plotSelFeaPCA(path, feafile, imgname, title, resfile, clf)
        
        resfile.close()


def runexps(X, y, dlvl, varnames, topk, mdltype="LR", preproctype="Normal", k=5):
    """ Run set of experiments based on variant inputs for evaluation of predictive ability of genes. 

    Parameters
    ----------
    X : double numpy matrix.
        Feature matrix where rows correspond to different perturbations and columns to the list of genes.
    y : double numpy array.
        Target variable containing assigned LDH% scores for regression or predication task.
    varnames    : numpy array.
        Full genes list in the data that corresponds to the assigned variable names for columns in X.
        
        
    Returns
    -------
        (Void) : No returns
    """
    all_results = {}
    all_scores = {}
    sel_features = {}
    
    feature_groups = ["All", "feasel4", "T1000", "L1000", "S1500", "Limma", "CD"]#["feasel2", "feasel1", "All", "S1500", "Limma", "EcoTox-M"]#,"feasel2", "feasel3" 
    input_files = {"L1000":"L1000_rat.txt", "T1000":"T1000_1000chip_rat.txt", "feasel1": "RFE", "Limma": "limma_genes_rat.txt", "CD": "CD_genes_rat.txt", "feasel4": "Random", "DeMAND":"DeMAND_genes_top1000.txt", "S1500": "S1500_rat.txt", "Manual": "EcoToxChip.txt"}#{"S1500": "S1500.txt","Limma":"limma_genes.txt", "EcoTox-M": "EcoToxChip.txt", "feasel1":"T1000", "feasel2":"Limma", "feasel3":"CD"} #
    
    for i in range(0, len(feature_groups)):
        all_sel_features = []
        f = feature_groups[i]
        if(f.find("feasel") != -1):# Case with feature selection
            resfname = "Regression/Results/"+input_files[f]+"_"+mdltype+"_topk"+str(topk)+"_"+preproctype+"_results.txt"
            scorefname = "Regression/Results/"+input_files[f]+"_"+mdltype+"_topk"+str(topk)+"_"+preproctype+"_scores.txt"
            genesfname = "Regression/Results/"+input_files[f]+"_topk"+str(topk)+"_"+preproctype+"_genes.txt"
            
            if(os.path.isfile(genesfname)):
                all_sel_features = np.array(pd.read_csv(genesfname, sep=" ", header=-1))
            
            if(not(os.path.isfile(resfname))):
                results, y_true, y_pred, all_sel_features  = evaluateperf(X, y, dlvl, mdltype, preproctype, k, topk, input_files[f], varnames, all_sel_features)                    
                np.savetxt(resfname, results, fmt="%s")
                np.savetxt(genesfname, all_sel_features, fmt="%s")
                np.savetxt(scorefname, np.array([y_true, y_pred]), fmt="%s")
                all_scores[input_files[f]] = [y_true, y_pred]
                all_results[input_files[f]] = results
            else:
                all_scores[input_files[f]] = np.array(pd.read_csv(scorefname, sep=" ", header=-1))
                all_results[input_files[f]] = np.array(pd.read_csv(resfname, sep=" ", header=-1))
        elif(f != "All"):# Case with a predefined set of features
            resfname = "Regression/Results/"+f+"_"+mdltype+"_topk"+str(topk)+"_"+preproctype+"_results.txt"
            scorefname = "Regression/Results/"+f+"_"+mdltype+"_topk"+str(topk)+"_"+preproctype+"_scores.txt"
            if(not(os.path.isfile(resfname))):
                if(f != "S1500"):
                    fea_idx = featureidx(varnames, 'data/Gene_sets/'+input_files[f])
                    fea_idx = np.unique(fea_idx)
                    fea_idx = fea_idx[0:topk]
                else:
                    fea_idx = featureidx(varnames, 'data/Gene_sets/'+input_files[f])
                results, y_true, y_pred, all_sel_features  = evaluateperf(X[:, fea_idx], y, dlvl, mdltype, preproctype, k)
                np.savetxt(resfname, results, fmt="%s")
                np.savetxt(scorefname, np.array([y_true, y_pred]), fmt="%s")
                all_scores[f] = [y_true, y_pred]
                all_results[f] = results
                
            else:
                all_scores[f] = np.array(pd.read_csv(scorefname, sep=" ", header=-1))
                all_results[f] = np.array(pd.read_csv(resfname, sep=" ", header=-1))
        else: # Case with all features
            resfname = "Regression/Results/allgenes"+"_"+mdltype+"_"+preproctype+"_results.txt"
            scorefname = "Regression/Results/allgenes"+"_"+mdltype+"_"+preproctype+"_scores.txt"
            if(not(os.path.isfile(resfname))):
                results, y_true, y_pred, all_sel_features = evaluateperf(X, y, dlvl, mdltype, preproctype, k)
                np.savetxt(resfname, results, fmt="%s")
                np.savetxt(scorefname, np.array([y_true, y_pred]), fmt="%s")
                all_scores[f] = [y_true, y_pred]
                all_results[f] = results
            else:
                all_scores[f] = np.array(pd.read_csv(scorefname, sep=" ", header=-1))
                all_results[f] = np.array(pd.read_csv(resfname, sep=" ", header=-1))
        
    
    summarybarplot(all_results, len(feature_groups), mdltype+"_topk"+str(topk)+"_"+preproctype) # 7 indicates total number of evaluation scores used
    summaryqqplot(all_scores, mdltype+"_topk"+str(topk)+"_"+preproctype)
    
    
    
    
    
