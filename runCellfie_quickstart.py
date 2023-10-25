#Function to run cellfie using the data and parameters from the quickstart on the wiki


import scipy.io
import numpy as np
from CellFie import CellFie

mat_data = scipy.io.loadmat('dataTest.mat')

gene_data = mat_data['data'][0, 0]['gene']

tissue_data = mat_data['data'][0, 0]['Tissue']

value_data = mat_data['data'][0, 0]['value']
data_void = mat_data['data'][0, 0]

data= {
    'gene': data_void['gene'],
    'Tissue': data_void['Tissue'],
    'value': data_void['value']
}


SampleNumber=3

ref='MT_recon_2_2_entrez'

param={
    'ThreshType':'local',
    'LocalThresholdType':'minmaxmean',
    'percentile_or_value':'percentile',
    'percentile_low': 25,
    'percentile_high': 75
}

pathway='metabolism'

sampNames=['samp1','samp2','samp3','samp4']

scores,binaryScores,taskInfos=CellFie(data, SampleNumber, ref, param, pathway, sampNames)

#write csv files for output data
np.savetxt('scores_output.csv', scores, delimiter=',')
np.savetxt('binary_scores_output.csv', binaryScores, delimiter=',')