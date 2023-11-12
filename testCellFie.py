import scipy.io
import numpy as np
from CellFie import CellFie

#generation of all potential results


#loading/preparing data file
mat_data = scipy.io.loadmat('dataTest.mat')
data_void = mat_data['data'][0, 0]
data= {
    'gene': data_void['gene'],
    'Tissue': data_void['Tissue'],
    'value': data_void['value']
}

#setting params that remain constant
SampleNumber=3
ref='MT_recon_2_2_entrez'
pathway='metabolism'
sampNames=['samp1','samp2','samp3','samp4']



#dataRecon22_local_minmaxmean_value
param={
    'ThreshType':'local',
    'LocalThresholdType':'minmaxmean',
    'percentile_or_value':'value',
    'value_low': 25,
    'value_high': 75
}

scores,binaryScores,taskInfos=CellFie(data, SampleNumber, ref, param, pathway, sampNames)
np.savetxt('dataRecon22_local_minmaxmean_value_score.csv', scores, delimiter=',')
np.savetxt('dataRecon22_local_minmaxmean_value_score_binary.csv', binaryScores, delimiter=',')





#dataRecon22_local_mean
param={
    'ThreshType':'local',
    'LocalThresholdType':'mean'
}

scores,binaryScores,taskInfos=CellFie(data, SampleNumber, ref, param, pathway, sampNames)
np.savetxt('dataRecon22_local_mean_score.csv', scores, delimiter=',')
np.savetxt('dataRecon22_local_mean_score_binary.csv', binaryScores, delimiter=',')

#dataRecon22_global_value
param={
    'ThreshType':'global',
    'percentile_or_value':'value',
    'value':50
}

scores,binaryScores,taskInfos=CellFie(data, SampleNumber, ref, param, pathway, sampNames)
np.savetxt('dataRecon22_global_value_score.csv', scores, delimiter=',')
np.savetxt('dataRecon22_global_value_score_binary.csv', binaryScores, delimiter=',')

#dataRecon22_global_percentile
param={
    'ThreshType':'global',
    'percentile_or_value':'percentile',
    'percentile':50
}

scores,binaryScores,taskInfos=CellFie(data, SampleNumber, ref, param, pathway, sampNames)
np.savetxt('dataRecon22_global_percentile_score.csv', scores, delimiter=',')
np.savetxt('dataRecon22_global_percentile_score_binary.csv', binaryScores, delimiter=',')





