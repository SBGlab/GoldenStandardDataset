import numpy as np

#round a number to nearest half (this assuming we measure each 0.5, could be optimized)
def round_half(number):
        return (round(number * 2) / 2)

def algorithm_iteration(data, iterations, doub_time, bouts, minimal_points=3, **kwargs):
    
    #INITAL VARIABLES-------------------------------------------------------------------
    
    time = np.arange(0, len(data)*bouts, bouts)
    #convert dataframe to proper format for np calculations
    data = data.astype(float)
    
    #calculations within the input dataframe (OD, log2 for doubling time, ln for growthrate)
    sample_od = list(data)
    sample_log2 = list(np.log2(data))
    sample_ln = list(np.log(data))
    
    #create the embedded growthrate data with np.gradient function (differentiates with finite differences)
    gr_grad = np.gradient(sample_log2, bouts)
    gr_list = list(gr_grad)
    gr_ln_grad = np.gradient(sample_ln, bouts)
    gr_ln_list = list(gr_ln_grad)

    #for index in range(1, len(data)):
    #    growthrate = (sample_log2[index] - sample_log2[index-1])/bouts
    #    gr_list.append(growthrate)

    #find the maximum growthrate
    max_gr = max(gr_list)
    max_gr_idx = gr_list.index(max_gr)

    #FUNCTION-----------------------------------------------------------------------
    two_doubling_time = (1/max_gr) * doub_time
    #print('FIRST DT:', two_doubling_time)  

    #checking dt value. If 0 raise error, minimal points defaults to 3 for defining regression        
    if two_doubling_time < 0:
        for ar in kwargs:
            sample_name = kwargs[ar]
        #print('Could not find a proper window for %s. Doubling time < 0' %sample_name)
        return None

    #WINDOW CALCULATION BASED ON COMPUTED 2 * DOUBLING TIME-----------------------------------------------

    #compute window of points spanning 2 * doubling time (just a range of points)
    near_half = round_half(two_doubling_time)
    window = int(near_half/bouts) #how many points define the window
    #print(window)
    #if 2 * doubling times is less than set minimal points, use minimal points
    if window < minimal_points:
        window = minimal_points

    #print('WINDOW: ', window)
    #compute the actual points within that range that contain the max growthrate point
    windows_with_max = [] #stores the windows of dOD/dt containing the max gr point
    for i in range(len(gr_list)-(window-1)): #loops through every n-point window, being n defined previously on 'window' variable
        window_temp = gr_list[i:i+window]
        if gr_list[max_gr_idx] in window_temp: #append only the dOD/dt windows that contain the max gr point
            windows_with_max.append((window_temp, 
                                     list(range(i,i+window)), 
                                     #'iteration '+ str(x)))
                                    ))
    #inputs for the output result. Warning: if regression results are computed, notice that are calculated with list indices,
    #not with the actual values. Numeric results for slope and intercept will be different. Better to compute them correctly where necessary

    OD_window = [] #output container
    for w in windows_with_max:
        window_idx = w[1] #list of indices (now we no longer have a list shifted by one)
        #window_add1 = [x+1 for x in window_idx] #indices for raw data-related list with original indices

        OD_window.append((np.var(w[0]),                            #variance of gr (computed as log2)
                          window_idx,                              #window indices
                          [time[d] for d in window_idx],           #time values of window
                          [sample_od[d] for d in window_idx],      #OD values
                          [sample_log2[d] for d in window_idx],    #log2OD values
                          [sample_ln[d] for d in window_idx],      #lnOD values
                          [gr_list[d] for d in window_idx],        #growth rate values (computed as log2OD)
                          [gr_ln_list[d] for d in window_idx],     #growth rate values (computed as lnOD)
                          near_half,
                          #stats.linregress(window_add1, ([sample_ln[d] for d in window_add1]))[0],
                          #stats.linregress(window_add1, ([sample_ln[d] for d in window_add1]))[2],
                          #stats.linregress(window_add1, ([sample_ln[d] for d in window_add1]))[1],
                          #'iteration '+ str(x)))
                         ))
    #filter all the window points looking for the minimum variance in log2(OD)

    try:
        min_var = min(OD_window)
        idx = OD_window.index(min_var)
        output = OD_window[idx]
    except:
        for ar in kwargs:
            sample_name = kwargs[ar]
        print('Could not find a proper window for %s. Doubling time > timespan of experiment' %sample_name)
        return None

    return output