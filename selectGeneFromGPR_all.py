def selectGeneFromGPR_all(model, gene_names, gene_exp, parsedGPR, minSum):
    import numpy as np
    import pandas as pd
    
    expressionCol_expected=pd.read_csv('expressionRxns.csv')
    gene_used_expected=pd.read_csv('gene_used.csv')
    num_rxns = len(model['rxns'][0][0])
    numSamp = gene_exp.shape[1]
    
    # initialize dataframe gene_used with num_rxns rows and numSamp columns, all values are empty strings
    gene_used = pd.DataFrame(np.full((num_rxns, numSamp), np.nan, dtype=object))
    gene_used.columns = [f'gene_used{str(x+1)}' for x in range(numSamp)]
    SampStr = [str(zz) for zz in range(1, numSamp + 1)]
    # initialize dataframe expressionCol with num_rxns rows and numSamp columns, all values are -1
    expressionCol = pd.DataFrame(np.full((num_rxns, numSamp), -1.0))
    expressionCol.columns = [f'expressionRxns{str(x+1)}' for x in range(numSamp)]
    for i in range(len(model['rxns'][0][0])):
        curExprArr = parsedGPR[i][0][0]
        tmpStuc = {}
        for zz in range(numSamp):
            tmpStuc[f'curExpr{SampStr[zz]}'] = []
            tmpStuc[f'gene_potential{SampStr[zz]}'] = []
        for j in range(len(curExprArr)):
            if len(curExprArr[j]) >= 1:
                geneID = np.where(np.isin(gene_names, curExprArr[j]))[0]
                if len(geneID) > 0:
                    for zz in range(numSamp):
                        if minSum:# I don't think this ever gets called
                            tmpStuc[f'curExpr{SampStr[zz]}'].append(np.sum(gene_exp[geneID, zz]))
                            tmpStuc[f'gene_potential{SampStr[zz]}'].extend(gene_names[geneID, zz])
                        else:
                            minGenevalue = np.min(gene_exp[geneID, zz])
                            minID = np.argmin(gene_exp[geneID, zz])
                            tmpStuc[f'curExpr{SampStr[zz]}'].append(minGenevalue)
                            tmpStuc[f'gene_potential{SampStr[zz]}'].append(gene_names[geneID[minID]])
                           

        for samp in range(numSamp):
            if tmpStuc[f'curExpr{SampStr[samp]}']:
                if minSum:# I don't think this ever gets called
                    expressionCol.loc[i, f"expressionRxns{samp+1}"] = np.min(tmpStuc[f'curExpr{SampStr[samp]}'])
                    gene_used.loc[i, f"gene_used{samp+1}"] = tmpStuc[f'gene_potential{SampStr[samp]}'][
                        np.argmin(tmpStuc[f'curExpr{SampStr[samp]}'])][0]

                else:
                    expressionCol.loc[i, f"expressionRxns{samp+1}"] = np.max(tmpStuc[f'curExpr{SampStr[samp]}'])
                    gene_used.loc[i, f"gene_used{samp+1}"] = tmpStuc[f'gene_potential{SampStr[samp]}'][
                        np.argmax(tmpStuc[f'curExpr{SampStr[samp]}'])][0]

                   
    gene_used = gene_used.astype(float)
    return expressionCol, gene_used