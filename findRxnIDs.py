
def findRxnIDs(model,rxnList):
    import numpy as np
    

    
    
    # check if rxnList has multiple values
    # returns list of indicies as rxnID
    

    if len(rxnList)>1:
        tmp=np.isin(rxnList,model)
        rxnID_tmp=[]
        for i in range(len(tmp)):
            if tmp[i]==True:
                rxnStr=rxnList[i][0][0]
                rxnID_tmp.append(np.where(model==rxnList[i])) 
            else:
                rxnID_tmp.append([]) 
        
        #making rxn into a list of integers for indexing(theres a better way to do this probably)
        rxnID=[]
        for i in range(len(rxnID_tmp)):
            try:
                rxnID.append(rxnID_tmp[i][0][0])
            
            except IndexError:
                rxnID.append([]) 
                
            
            
                
    
    #excecutes next bit if rxnList has one or zero elements
    #gets list of indicies of where your string shows up in the array
    #if this is an empty list, set rxnID to empty list, return
    #if list is longer than one, return the first element as rxnID
    #if its just one thing, return the index
    else:
        rxnList=rxnList[0][0][0]
        rxnID_tmp=np.where(model==rxnList)[0] 
        if len(rxnID_tmp)==0:
            rxnID_tmp=[]
        if len(rxnID_tmp)>1:
            rxnID_tmp=rxnID_tmp[0]
            
        #getting index as integer (or empty list)
        try:
            rxnID=rxnID_tmp[0]
            rxnID=[rxnID]

            
            
        except IndexError:
            rxnID=[] 
        
        
    return rxnID #this is gonna be a list of indicies if the first bit excecutes or a single one