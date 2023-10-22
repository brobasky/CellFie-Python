def selectGeneFromGPR_all(model, gene_names, gene_exp, parsedGPR, minSum):
    '''
    % Map gene expression to reaction expression using the GPR rules. An AND
    % will be replaced by MIN and an OR will be replaced by MAX.
    %
    % USAGE:
    %   expressionCol = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR, minMax)
    %
    % INPUTS:
    %   model:          COBRA model struct
    %   gene_names:     gene identifiers corresponding to gene_exp. Names must
    %                   be in the same format as model.genes (column vector)
    %                   (as returned by "findUsedGeneLevels.m")
    %   gene_exp:       gene FPKM/expression values, corresponding to names (column vector)
    %                   (as returned by "findUsedGeneLevels.m")
    %   parsedGPR:      GPR matrix as returned by "GPRparser.m"
    %
    % OPTIONAL INPUTS:
    %   minSum:         instead of using min and max, use min for AND and Sum
    %                   for OR
    %
    % OUTPUTS:
    %   expressionCol:  reaction expression, corresponding to model.rxns.
    %                   No gene-expression data and orphan reactions will
    %                   be given a value of -1.
    %
    % AUTHOR: Anne Richelle, May 2017
    '''

    #placeholder code

    import numpy as np
    import pandas as pd

    expressionCol=pd.read_csv('expressionRxns.csv')
    gene_used=pd.read_csv('gene_used.csv')



    '''
    if minSum is None:
        minSum = False

    numSamp = gene_exp.shape[1]
    gene_used = {}

    for i in range(len(model['rxns'])):
        for zz in range(numSamp):
            gene_used[(i, zz)] = ''

    SampStr = [str(zz) for zz in range(1, numSamp + 1)]

    expressionCol = np.full((len(model['rxns']), numSamp), -1.0)

    for i in range(len(model['rxns'])):
        curExprArr = parsedGPR[i]
        tmpStuc = {}

        for zz in range(numSamp):
            tmpStuc[f'curExpr{SampStr[zz]}'] = []
            tmpStuc[f'gene_potential{SampStr[zz]}'] = []

        for j in range(len(curExprArr)):
            if len(curExprArr[j]) >= 1:
                geneID = np.where(np.isin(gene_names, curExprArr[j]))[0]

                if len(geneID) > 0:
                    for zz in range(numSamp):
                        if minSum:
                            tmpStuc[f'curExpr{SampStr[zz]}'].append(np.sum(gene_exp[geneID, zz]))
                            tmpStuc[f'gene_potential{SampStr[zz]}'].extend(gene_names[geneID, zz])
                        else:
                            minGenevalue = np.min(gene_exp[geneID, zz])
                            #debug
                            print('geneID:')
                            print(type(geneID))
                            print(len(geneID))
                            print(geneID)


                            print('gene_names:')
                            print(type(gene_names))
                            print(len(gene_names))
                            print(gene_names)
                            print(gene_names.shape)

                            print(f'numSamp:{numSamp}')
                            print(zz)

                            minID = np.argmin(gene_exp[geneID, zz])
                            tmpStuc[f'curExpr{SampStr[zz]}'].append(minGenevalue)
                            #issue since this is a column vector and you're trying to use 2 indicies
                            #gene_names array (source of error) is gene_id var output by 1st function
                            #tmpStuc[f'gene_potential{SampStr[zz]}'].append(gene_names[geneID[minID], zz])
                            tmpStuc[f'gene_potential{SampStr[zz]}'].append(gene_names[geneID[minID]])

        for zz in range(numSamp):
            if tmpStuc[f'curExpr{SampStr[zz]}']:
                if minSum:
                    expressionCol[i, zz] = np.min(tmpStuc[f'curExpr{SampStr[zz]}'])
                    gene_used[(i, zz)] = tmpStuc[f'gene_potential{SampStr[zz]}'][
                        np.argmin(tmpStuc[f'curExpr{SampStr[zz]}'])]
                else:
                    expressionCol[i, zz] = np.max(tmpStuc[f'curExpr{SampStr[zz]}'])
                    gene_used[(i, zz)] = tmpStuc[f'gene_potential{SampStr[zz]}'][
                        np.argmax(tmpStuc[f'curExpr{SampStr[zz]}'])]
                        
        '''

    return expressionCol, gene_used