import numpy as np
import pandas as pd 
import linecache



def las_to_df(file, depth_start, depth_stop, interval=True, nan=-999.25):
    '''read the las log file(exported from discovery or petrel) 
    to a DataFrame'''
    
    '''read the log name'''     
    lookup = '~C'
    with open(file, 'r') as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                lognamestart = num+1
    lookup = '~P'
    with open(file, 'r') as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                lognamestop = num
    lognames = []
    for num in np.arange(lognamestart, lognamestop, 1):
        line = linecache.getline(file, num)
        line = line.replace('.',' ')
        logname = line.split(' ')[0] #'.', 
        lognames.append(logname)

    '''read the data'''
    lookup = '~A'
    with open(file, 'r') as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                results = np.loadtxt(file, skiprows=num)
    logs = pd.DataFrame(results, columns=lognames)
    logs = logs.replace(nan, np.nan)
    if interval == True:
        if logs.iloc[:,0].min() > depth_start:
            print('Attention: depth_start should be larger')
        else: logs = logs[logs.iloc[:,0] >= depth_start]
    
        if logs.iloc[:,0].max() < depth_start:
            print('Attention: depth_stop should be smaller')
        else: logs = logs[logs.iloc[:,0] <= depth_stop]    
       
        logs = pd.DataFrame(np.array(logs), columns=lognames)
    
    return logs


