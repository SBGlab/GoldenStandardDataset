import numpy as np
import pandas as pd

def check_unique(query_list):
    if len(query_list) != len(set(query_list)):
        return True

# function to give the correct format to dataframes: wells in columns and timepoints for the indices. Also create a global time list
def correct_shape(dataframe, bouts):
    first = dataframe.columns[0]
    dataframe.index = list(dataframe[first])
    dataframe = dataframe.iloc[:, 1:]
    f = len(dataframe.columns) * bouts
    global time
    time = np.arange(0, f, bouts)
    dataframe.columns = time
    dataframe = dataframe.transpose()
    return dataframe

def len_list(list_data):
    l = len(list_data)
    if l%2 != 0:
        idx = int(l/2)
        return list_data[idx]
    else:
        idx1 = int(l/2)-1
        idx2 = int(l/2)
        return((list_data[idx1]+list_data[idx2])/2)

#round a number to nearest half (this assuming we measure each 0.5, could be optimized)
def round_half(number):
        return (round(number * 2) / 2)

#MAIN FUNCTION FOR KEREN ET AL. ALGORITHM
    #We apply the algorithm to individual samples instead of the whole plate
    #We minimize the variance in growthrate instead of dOD/dt
    #As there is the possibility that 2*doubling times exceeds the time of actual exponential phase, we can choose how many doubling times to use

def normalization(list_input, normalization_idx, value=1):
    norm_target = list_input[normalization_idx]
    list_output = []
    for data in list_input:
        if data == norm_target:
            norm_value = value
        else:
            norm_value = (data * value)/norm_target
        list_output.append(norm_value)
    return list_output

def y_fit_points(x_data, fit_data):
    deg = len(fit_data)
    y_data = 0
    for coef, exp in zip(range(0, deg), reversed(range(0, deg))):
        y_data = y_data + fit_data[coef] * x_data ** exp
    return y_data

def fit_equation(fit_data):
    deg = len(fit_data)
    eq = '0'
    for coef, exp in zip(range(0, deg), reversed(range(0, deg))):
        eq = eq  + ' + %f * x ** %d' %(fit_data[coef], exp)
    return eq

def exp_func(x, slope, intercept):
    A = np.exp(intercept)
    m = slope
    return A * np.exp(m * x)

def fit_func(x, coef1, coef2, coef3):
    return coef1 * x ** 2 + coef2 * x ** 1 + coef3 * x ** 0

def text_to_sample(string):
    samples = []    
    open_par = string.find('(') + 1
    close_par = string.find(')')    
    while True:
        out = string[open_par:close_par]
        if 'Samples' in out:
            break
        samples.append(out)        
        open_par = string.find('(', close_par) + 1
        close_par = string.find(')', open_par)
    return samples

def sample_to_wells(string):
    wells = []
    ap = string.find("'") + 1
    ap_next = string.find("'", ap)
    while True:
        out = string[ap:ap_next]
        if out == '':
            break
        elif ' ' in out:
            ap = string.find("'", ap_next) + 1
            ap_next = string.find("'", ap)
        else:
            ap = string.find("'", ap_next) + 1
            ap_next = string.find("'", ap)
            wells.append(out)
    return wells