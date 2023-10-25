def findUsedGenesLevels_all(model, exprData, print_level=0):
    '''
    % Returns vectors of gene identifiers and corresponding gene expression
    % levels for each gene present in the model ('model.genes').
    %
    % USAGE:
    %    [gene_id, gene_expr] = findUsedGenesLevels(model, exprData)
    %
    % INPUTS:
    %
    %   model:               input model (COBRA model structure)
    %
    %   exprData:            mRNA expression data structure
    %       .gene                cell array containing GeneIDs in the same
    %                            format as model.genes
    %       .value               Vector containing corresponding expression value (FPKM)
    %
    % OPTIONAL INPUTS:
    %    printLevel:         Printlevel for output (default 0);
    %
    % OUTPUTS:
    %
    %   gene_id:             vector of gene identifiers present in the model
    %                        that are associated with expression data
    %
    %   gene_expr:           vector of expression values associated to each
    %                        'gened_id'
    %
    %
    % Authors: - S. Opdam & A. Richelle May 2017
    '''


    import numpy as np

    gene_expr = []
    gene_id = model['genes'][0][0]
    tmpb = len(exprData['value'][0])

    for cur_ID in gene_id:
        dataID=[]
        if cur_ID in exprData['gene']:
            dataID.append(exprData['gene'].index(cur_ID))
        if len(dataID) == 0:
            gene_expr.append(-np.ones(tmpb))
        elif len(dataID) == 1:
            gene_expr.append(exprData['value'][dataID[0]])
        elif len(dataID) > 1:
            if print_level > 0:
                print(f'Double for {cur_ID}')
            gene_expr.append(np.mean(exprData['value'][dataID], axis=0))


    return gene_id, np.array(gene_expr)