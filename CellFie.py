def CellFie(data,SampleNumber,ref,param, pathway, sampNames):
    import scipy.io
    import pandas as pd
    import numpy as np
    import math
    from findRxnIDs import findRxnIDs
    from findUsedGenesLevels_all import findUsedGenesLevels_all
    from selectGeneFromGPR_all import selectGeneFromGPR_all
    import warnings
    
    
    warnings.filterwarnings("ignore")
    
    #checking that input parameters are valid
    if data['value'].shape[1] != SampleNumber:
        raise ValueError('The number of samples defined is not the same as the size of the dataset')
    if data['value'].shape[0] != len(data['gene']):
        raise ValueError('data.value does not have the same number of rows as data.gene')
    if 'ref' not in locals() and 'ref' not in globals():
        raise ValueError('The reference model has not been defined - please choose a reference model')
    if 'param' not in locals() and 'param' not in globals():
        param = {}
        param['ThreshType'] = 'local'
        param['percentile_or_value'] = 'percentile'
        param['LocalThresholdType'] = 'minmaxmean'
        param['percentile_low'] = 25
        param['percentile_high'] = 75
    #handling what to do if 'param' not provided
    if 'param' not in locals() and 'param' not in globals():
        param = {}
        param['ThreshType'] = 'local'
        param['percentile_or_value'] = 'percentile'
        param['LocalThresholdType'] = 'minmaxmean'
        param['percentile_low'] = 25
        param['percentile_high'] = 75
        
    #loading taskstructure based on pathway arg (NEED TO CHANGE FILEPATHS)
    if pathway == "secretion":
        taskStructure = scipy.io.loadmat('taskStructure_sec.mat')
    elif pathway == "metabolism":
        taskStructure = scipy.io.loadmat('taskStructure.mat')
        taskInfos = taskStructure['taskStructure']
        taskInfos = [list(taskInfo) for taskInfo in taskInfos]
        taskInfos = [taskInfo[:5] for taskInfo in taskInfos]
    else:
        raise ValueError('Pathway not valid')
        
    #loading essential reactions based on pathway arg(CHANGE FILEPATHS)
    if pathway==('secretion'):
        essentialRxns=scipy.io.loadmat(f'essentialRxnsbyTask_{ref}_sec.mat')
    elif pathway==('metabolism'):
        essentialRxns=scipy.io.loadmat(f'essentialRxnsbyTask_{ref}.mat')
    else:
        raise ValueError('Pathway not valid')
    
    #check if gene in input data are present in model
    #count genes not included and remove them from data
    model=scipy.io.loadmat(ref)
    model=model['model']
    model['genes']=model['genes'][0]
    data_genes=data['gene']
    model_genes=model['genes']
    model_genes=model['genes'][0][0]
    
    
    ID_model = set()
    gene_notInModel = set() #indicies of data['gene'] where a gene is not in the model
    
    
    for i in range(len(data_genes)):
        if not np.isin(data_genes[i], model_genes):
            gene_notInModel.add(i)
        else:
            tmpid = np.where(np.array(data_genes[i]) == model_genes)[0]
            ID_model.add(tmpid[0])
    
    if not gene_notInModel:
        print('All genes provided in data are included in the reference model')
    else:
        print(str(len(gene_notInModel)) + ' genes provided are not included in the reference model:')
        #print(data_genes[gene_notInModel])
    # Remove genes and associated values that are not in the model    
    data['gene'] = [gene for i, gene in enumerate(data['gene']) if i not in gene_notInModel]
    if SampleNumber == 1:
        data['value'] = [value for i, value in enumerate(data['value']) if i not in gene_notInModel]
    else:
        data['value'] = np.delete(data['value'], list(gene_notInModel), axis=0)
    
    
    
    #getting threshold value for data
    #making histogram, print figure
    
    # Get the threshold value and the histogram for the complete dataset
    if SampleNumber > 1:
        linData = data['value'].reshape(-1, 1)
    else:
        linData = data['value']
    # Remove zeros from linData
    linData = linData[linData != 0]
    # Calculate the threshold value based on your conditions
    if param['ThreshType'] == 'global' and param['percentile_or_value'] == 'percentile':
        print('RUN - global: percentile')
        l_global = np.percentile(np.log10(linData), param['percentile'])
        data['ths'] = 10 ** l_global
    elif param['ThreshType'] == 'global' and param['percentile_or_value'] == 'value':
        print('RUN - global: value')
        data['ths'] = param['value']
    elif param['ThreshType'] == 'local' and param['LocalThresholdType'] == 'mean':
        print('RUN - local: mean')
    elif param['ThreshType'] == 'local' and param['LocalThresholdType'] == 'minmaxmean' and param[
        'percentile_or_value'] == 'percentile':
        print('RUN - local: minmaxmean: percentile')
        l_high = np.percentile(np.log10(linData), param['percentile_high'])
        data['ths_high'] = 10 ** l_high
        l_low = np.percentile(np.log10(linData), param['percentile_low'])
        data['ths_low'] = 10 ** l_low
    elif param['ThreshType'] == 'local' and param['LocalThresholdType'] == 'minmaxmean' and param[
        'percentile_or_value'] == 'value':
        print('RUN - local: minmaxmean: value')
        data['ths_high'] = param['value_high']
        data['ths_low'] = param['value_low']
    else:
        raise ValueError('No analysis triggered')
        
        
    
    #computing thresholds based on selected approach
    Gene_score = np.zeros_like(data['value'], dtype=float)  # Initialize Gene_score as an array of zeros
    if param['ThreshType'] == 'local':
        if param['LocalThresholdType'] == 'mean':
            # The threshold for each gene is equal to its mean value over all the samples
            threshold = np.mean(data['value'], axis=1)
            #handling future divide by zero errors to produce the same result as matlab code
            #does not produce the exact same scores but very close
            threshold[threshold==0.0] = 1e-100
            
        else:
            threshold = np.empty(len(data['gene']), dtype=float)
            for i in range(len(data['gene'])):
                expressionValue = data['value'][i, :]
                if np.mean(expressionValue) >= data['ths_high']:
                    threshold[i] = data['ths_high']
                else:
                    threshold[i] = max(np.mean(expressionValue), data['ths_low'])
    
        # Every single gene is associated with an expression score
        for i in range(SampleNumber):
            Gene_score[:, i] = 5.0 * np.log(1.0 + data['value'][:, i] / threshold)
    
    elif param['ThreshType'] == 'global':
        Gene_score = 5.0 * np.log(1.0 + data['value'] / data['ths'])
    
    # Table of gene scores to output
    geneScore = np.column_stack((data['gene'], Gene_score))
    geneScore_table = pd.DataFrame(geneScore, columns=sampNames)
    
    
    
    #mapping expression data to model, making expression dict, finding used gene levels
    
    # initializing a dictionary for expression data
    expression = {
        'gene': data['gene'],
        'Rxns': [],
        'gene_used': [],
        'count': []
    }
    
    minSum=False
    # Load parsedGPR for the model based on 'pathway'
    if pathway == "metabolism":
        parsedGPR = scipy.io.loadmat(f'parsedGPR_{ref}.mat')
    else:
        parsedGPR = scipy.io.loadmat(f'parsedGPR_{ref}_sec.mat')

    # Extract the parsedGPR data from the loaded data
    parsedGPR = parsedGPR['parsedGPR']
    expression['value']=Gene_score
    gene_id, gene_expr = findUsedGenesLevels_all(model, expression)
    
    
    #linking gene to model reactions
    #calling selectGeneFromGPR_all and manipulating data
    expressionRxns, gene_used = selectGeneFromGPR_all(model, gene_id, gene_expr, parsedGPR, minSum)


    #loop filling out count feild in expression dict
    #counting how many reactions each gene is associated with
    #loop also gets data for Rxns and gene_used feilds

    for zz in range(SampleNumber):
        gene_all=[] #used to store gene IDs
        for j in range(len(gene_used)):
            #check if cell in gene_used at j,zz is not empty
            if gene_used.iloc[j][zz]>=0:
                gene_all.append(gene_used.iloc[j][zz])
        gene_all_series=pd.Series(gene_all)
        countGene=gene_all_series.value_counts().reset_index()
        countGene.columns=['Value','Count']
        countGene['Percent']=(countGene['Count']/len(gene_all_series))*100
        count=[]
        #going thru the rows in gene_used again for column zz
        for k in range(len(gene_used)):
            #check if theres a gene associated with the reaction
            if gene_used.iloc[k][zz]>=0:
                #get count of association for specific gene
                gene_value = float(gene_used.iat[k, zz])
                tmp = countGene[countGene['Value']==gene_value]
                tmp=  tmp['Count'].item()
                count.append(tmp)
            else:
                #set count to 0 for no associated genes
                count.append(0)
        expression['count'].append(count)
        
    expression['Rxns']=expressionRxns
    expression['gene_used']=gene_used
    expression['count']=np.array(expression['count'])
    
    #cleaning essential reactions/ model reactions
    essentialRxns=essentialRxns['essentialRxns']
    essentialRxns=essentialRxns[0]
    model_rxns=model['rxns'][0][0]
    
    #initializing variables for score calculations
    expressionRxns=expression['Rxns']
    significance=1/expression['count']
    significance=significance.T
    #these next 2 lines handle dividing by 0 in matlab, maybe unnecessary in python 
    significance[np.isinf(significance)] = 0
    significance[np.isneginf(significance)] = 0
    ScorebyTask=[None]*len(taskInfos)
    ScorebyTask_binary=np.full((len(taskInfos),SampleNumber), None, dtype=object)
    
    #load np array of expected binary scores for testing
    binary_expected=np.genfromtxt('score_binary_expected.csv',delimiter=',')
    
    #computing scores/binary scores

    #looping through 195 indicies of taskInfos
    for i in range(len(taskInfos)):
        #checking if theres an associated rxn for index
        if len(essentialRxns[i])>0: #essentialrxns and taskinfos line up by index
            #finding each essential rxn's index in the model
            rxns=essentialRxns[i]
            rxnID=findRxnIDs(model_rxns,rxns)
            #checking if the rxn (essentialRxns[i]) is in the model
            if rxnID!=[]: 
                #if it is in the model, pull matching records from expressionrxns
                expValue=expressionRxns.loc[rxnID].reset_index(drop=True)
                if isinstance(expValue,pd.Series):
                    expValue=pd.DataFrame([expValue])
                #also pulling matching records from significance
                signValue=significance[rxnID]
                
                
                #if no gene is associated with a reaction, remove reation from count
                if (expValue.sum(axis=1)==-SampleNumber).any():
                    row_sums = np.sum(expValue, axis=1) 
                    indices = np.where(row_sums == -SampleNumber)[0]
                    signValue = np.delete(signValue, indices, axis=0) #deleting rows by index
                    expValue=expValue.drop(indices,axis=0)
                         
                #checking if expValue is empty
                if len(expValue)>0:
                    #updating score lists if expValue df has multiple rows
                    if len(expValue)>1:
                        #takes averages if expValue has more than 1 row
                        ScorebyTask[i]=((expValue.values * signValue).sum(axis=0))/expValue.shape[0]
                        Val=expValue.sum(axis=0)/expValue.shape[0] # .values?
                        #returning indicies where val is above a threshold
                        ID_up=np.where(Val>=(5*math.log(2)))[0]                       
                        #setting binary scores
                        ScorebyTask_binary[i]=[0]*SampleNumber
                        ScorebyTask_binary[i][ID_up]=1
                    
                    #updating score lists if expValue df is just 1 row
                    else: 
                        
                        ScorebyTask[i]=((expValue.values * signValue).sum(axis=0))/expValue.shape[0]
                        Val=expValue.sum(axis=0)/expValue.shape[0] 
                        ID_up=np.where(Val>=(5*math.log(2)))[0]
                        ScorebyTask_binary[i]=[0]*SampleNumber
                        ScorebyTask_binary[i][ID_up]=1
                        
                                            
                #updating score lists with -1's if expValue df is empty   
                else:
                    ScorebyTask[i]=[-1]*SampleNumber
                    ScorebyTask_binary[i]=[-1]*SampleNumber
            
            #updating score lists with -1's if rxnID is empty
            else:
                ScorebyTask[i]=[-1]*SampleNumber
                ScorebyTask_binary[i]=[-1]*SampleNumber
                
       
        #updating score lists with -1's if essentialRxns[i] is empty
        else:
            ScorebyTask[i]=[-1]*SampleNumber
            ScorebyTask_binary[i]=[-1]*SampleNumber
           
    ScorebyTask=np.vstack(ScorebyTask)
    taskInfos=np.vstack(taskInfos)
            
    
    return ScorebyTask, ScorebyTask_binary ,taskInfos
