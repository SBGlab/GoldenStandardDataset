import numpy as np
import pandas as pd

#function to check that elements do not repeat within a list
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

#function to get midpoint of a list
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

#round number to given bout, optimized version of round_half
def round_to_bout(number, bout):
    return (round(number * (1/bout)) * bout)

#MAIN FUNCTION FOR KEREN ET AL. ALGORITHM
    #We apply the algorithm to individual samples instead of the whole plate
    #We minimize the variance in growthrate instead of dOD/dt
    #As there is the possibility that 2*doubling times exceeds the time of actual exponential phase, we can choose how many doubling times to use

#normalize values of a list to a value of that list
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

#functions to generate plots from fits
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

#function to get sample info from the output string info
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

#match sample info to wells
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

#function to convert algorithm output info to actual values when the algorithm
#computes only considering some cropped timespan of the whole data
def index_crop_to_index(initial_crop_val, w_idcs_from_algorithm, time):    
    n = list(time).index(initial_crop_val)
    i = w_idcs_from_algorithm[0]
    r = w_idcs_from_algorithm[-1]
    k = list(time)[n+i]
    p = list(time).index(k)
    s = list(time)[n+r]
    j = list(time).index(s)
    w_idcs_no_crop = np.arange(p, j+1, 1)
    return w_idcs_no_crop

#list indx to X,2 shape coordinates
def convert_idx_to_coord(idx):
    row = int(idx/2)
    col = idx-2*row
    return [row, col]

wells_in_plate_96 = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
                  'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12',
                  'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12',
                  'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12',
                  'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12',
                  'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12',
                  'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12',
                  'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12']

wells_in_plate_384 = [	'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19', 'A20', 'A21', 'A22', 'A23', 'A24', 
						'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15', 'B16', 'B17', 'B18', 'B19', 'B20', 'B21', 'B22', 'B23', 'B24', 
						'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24',
						'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20', 'D21', 'D22', 'D23', 'D24', 
						'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E14', 'E15', 'E16', 'E17', 'E18', 'E19', 'E20', 'E21', 'E22', 'E23', 'E24', 
						'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23', 'F24',
						'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 
						'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'H17', 'H18', 'H19', 'H20', 'H21', 'H22', 'H23', 'H24', 
						'I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', 'I9', 'I10', 'I11', 'I12', 'I13', 'I14', 'I15', 'I16', 'I17', 'I18', 'I19', 'I20', 'I21', 'I22', 'I23', 'I24', 
						'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J19', 'J20', 'J21', 'J22', 'J23', 'J24', 
						'K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'K7', 'K8', 'K9', 'K10', 'K11', 'K12', 'K13', 'K14', 'K15', 'K16', 'K17', 'K18', 'K19', 'K20', 'K21', 'K22', 'K23', 'K24', 
						'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10', 'L11', 'L12', 'L13', 'L14', 'L15', 'L16', 'L17', 'L18', 'L19', 'L20', 'L21', 'L22', 'L23', 'L24', 
						'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21', 'M22', 'M23', 'M24', 
						'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15', 'N16', 'N17', 'N18', 'N19', 'N20', 'N21', 'N22', 'N23', 'N24', 
						'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10', 'O11', 'O12', 'O13', 'O14', 'O15', 'O16', 'O17', 'O18', 'O19', 'O20', 'O21', 'O22', 'O23', 'O24',
						'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24']

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:pink', 'tab:brown',
          'tab:olive', 'lightcoral', 'seagreen', 'mediumpurple', 'gold', 'sandybrown', 'deepskyblue', 'palegreen']

#function to get wells that have no samples in the experiment
def spot_removed_wells(trimmed_df_dict):
    df_key = list(trimmed_df_dict.keys())[0]
    wells_in_df = list(trimmed_df_dict[df_key].keys())
    wells_not_in_df = []
    for well in wells_in_plate:
        if well not in wells_in_df:
            wells_not_in_df.append(well)
    return wells_not_in_df

#output name for export file
def create_file(fname, format):
    return fname + format

#check correct format of input files
def check_file(file, f):
    file_name = list(file.keys())[0]
    if file_name.split('.')[1] == f:
        return True
    else:
        return False
