# General functions applicable to all techniques
import numpy as np
from scipy import interpolate

def yforx(x, xdata, ydata): #Used by other functions
    """Sorts xdata into ascending order (required by splrep) and solves y (as fa) for a value of x. Also returns spl to allow for further calculations to take place quicker. Used by other functions.

    Arguments
    ---------
    x : float
        x value you want to obtain the y for
    xdata : numpy.array
        numpy array containing the x data
    ydata : numpy.array
        numpy array containing the y data
    """

    #Sorts data into ascending order (required by splrep)
    data = xdata
    indices = np.argsort(data)
    data_index = data[indices]
    percent_index = ydata[indices]
    #interpolates data and solves y for a specific x
    spl = interpolate.splrep(data_index, percent_index, s=0) 
    #If errors, try s=1 smoothing to eradicate duplicates
    fa = interpolate.splev(x,spl,der=0)   # f(a)
    return fa, spl
    
def max_in_range(dataframe, dataframe_index, low, high):
    """Returns the maximum value of y as well as the corresponding x value in a given range of x (low to high)
    
    """
    signal = dataframe.values
    energy = dataframe_index.values
    # data = np.vstack((energy, signal))
    y_values = signal[np.logical_and(low < energy, energy < high)]
    x_values = energy[np.logical_and(low < energy, energy < high)]
    index_max_y = y_values.argmax()
    max_y = y_values[index_max_y]
    max_x = x_values[index_max_y]
    return max_x, max_y
