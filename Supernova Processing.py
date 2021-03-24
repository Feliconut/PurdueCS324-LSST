# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:30:13 2021

@author: blgnm
"""

import pandas as pd
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import numpy as np
import lmfit

#%%
def Data(Objectid):
    """
    Parameters
    ----------
    Objectid : String, ZTF ID
    
    Takes in ZTF object id and makes a dataframe of it's lightcurve data

    Returns
    -------
    DataFrame with lightcurve data

    """
    Supernovadf = pd.DataFrame(columns = ["mjd", "band", "magnitude", "error", "event", "class"])

    for i in Objectid:
        
        locus = get_by_ztf_object_id(i)
        Data = locus.lightcurve
    
        Data_frame1 = pd.DataFrame.from_dict(Data[Data['ant_passband']=='g']) 
        Data_frame2 = pd.DataFrame.from_dict(Data[Data['ant_passband']=='R']) 
        
        Data_frame1['ant_mag'] = Data_frame1['ant_mag'].replace(np.nan, 0)
        Data_frame2['ant_mag'] = Data_frame2['ant_mag'].replace(np.nan, 0)
        Data_frame1 = Data_frame1[Data_frame1.ant_mag > 0] 
        Data_frame2 = Data_frame2[Data_frame2.ant_mag > 0] 
    
        
        MJD1 = Data_frame1['ant_mjd']
        MJD2 = Data_frame2['ant_mjd']
        MagnitudeG = Data_frame1['ant_mag'] 
        MagnitudeR = Data_frame2['ant_mag']
    
        MJD1 = MJD1 - (MJD1.min() - 1)
        MJD2 = MJD2 - (MJD2.min() - 1)
        
        GBand = pd.DataFrame(columns = ["mjd", "band", "magnitude", "error", "event"])
        GBand["mjd"] = Data_frame1["ant_mjd"]
        GBand["band"] = pd.Series(np.zeros([len(MagnitudeG)]))
        GBand["magnitude"] = MagnitudeG    
        GBand['band'] = GBand['band'].replace(np.nan, 0)
        GBand['error'] = Data_frame1["ant_magerr"]
    
        RBand = pd.DataFrame(columns = ["mjd", "band", "magnitude", "error", "event"])
        RBand["mjd"] = Data_frame2["ant_mjd"]
        RBand["band"] = pd.Series(np.ones([len(MagnitudeR)]))
        RBand["magnitude"] = MagnitudeR  
        RBand['band'] = RBand['band'].replace(np.nan, 1)
        RBand['error'] = Data_frame2['ant_magerr']
        num = np.zeros(len(RBand))
        num1 = np.zeros(len(GBand))
                
        GBand['event'] = num1
        RBand['event'] = num
        GBand['event'] = GBand['event'].replace([0], [i])
        RBand['event'] = RBand['event'].replace([0], [i])
              
        Band = pd.concat([GBand, RBand], axis = 0, ).reset_index(drop=True)
        Supernovadf = pd.concat([Supernovadf, Band], axis = 0).reset_index(drop=True)
        
    return(Supernovadf)

#%%
import george
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.optimize import minimize
from astropy.stats import biweight_location
from scipy.optimize import curve_fit
import random
from astropy.stats import biweight_location



def Multi_Band_GP(x_range, x, y, y_err, dim, n_samples = False, sampling = False):
    """ Considers cross corrolations between multiple bands as dims, prone to holding the order of the bands too rigidly """
    """ Will optimize for 'best' parameters when given no parameters """
    """ x = mjd, y and y_err = measurment, dim and dim_err = wavelength in nm """
    length_scale = 20
    signal_to_noises = (np.abs(y) / np.sqrt(np.power(y_err,2) + (1e-2 * np.max(y))**2))
    scale = np.abs(y[signal_to_noises.argmax()])
    kernel = ((0.5 * scale)**2 * george.kernels.Matern32Kernel([length_scale**2, 6000**2], ndim=2))
    kernel.freeze_parameter('k2:metric:log_M_1_1')
    kernel.freeze_parameter('k1:log_constant') #Fixed Scale
    x_data = np.vstack([x, dim]).T
    gp = george.GP(kernel, mean = biweight_location(y))
    guess_parameters = gp.get_parameter_vector()
    gp.compute(x_data, y_err)
    x_pred = np.linspace(x.min(), x.max(), n_samples)
    x_pred = np.vstack([x, dim]).T
    pred, pred_var = gp.predict(y, x_pred, return_var=True)
    # bounds = [(0, np.log(1000**2))]
    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(y)
    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(y)
    result = minimize(
            neg_ln_like,
            gp.get_parameter_vector(),
            jac=grad_neg_ln_like,
            # bounds=bounds
            )
    if result.success:
            gp.set_parameter_vector(result.x)
    else:
            gp.set_parameter_vector(guess_parameters)    
    gp.set_parameter_vector(result.x)
    # print(kernel.get_parameter_vector(True))
    #print("\nFinal ln-likelihood: {0:.2f}".format(gp.log_likelihood(y)))
    if n_samples != False:
        x_pred = np.vstack([np.array(list(np.linspace(x_range.min(), x_range.max(), n_samples))*np.unique(dim).size),
                            np.array(np.sort(list(np.unique(dim))*n_samples))]).T 
        # x_pred = np.vstack([np.array(list(np.linspace(x_range.min(), x_range.max(), n_samples))*6),
        #                     np.array(np.sort([357, 476, 621, 754, 870, 1004]*n_samples))]).T 
        pred, pred_var = gp.predict(y, x_pred, return_var=True)
        output = [x_pred[:,0], pred, np.sqrt(pred_var), x_pred[:,1], []]
        return output
    elif sampling != False:
        x_pred = np.vstack([np.array(sampling[0]),
                            np.array(sampling[1])]).T  
        pred, pred_var = gp.predict(y, x_pred, return_var=True)
        output = [x_pred[:,0], pred, np.sqrt(pred_var), x_pred[:,1], []]
        return output

def band_to_color(inp):
    labels = [357, 476, 621, 754, 870, 1004]
    # labels = [0,1,2,3,4,5]
    labels_2=['green', 'red', 'goldenrod', 'blue', 'pink', 'grey']
    outlst = []
    for x in inp:
        out = labels.index(int(x))
        out = labels_2[out]
        outlst.append(out)
    return outlst

def band_to_wvlg(inp):
    labels = [0,1,2,3,4,5]
    labels_2=[357.0, 476.7, 621.5, 754.5, 870.75, 1004.0]
    outlst = []
    for x in inp:
        out = labels.index(int(x))
        out = labels_2[out]
        outlst.append(out)
    return outlst

def expfun(x, a, b):
    return np.multiply(np.exp(np.multiply(x, b)), a)

def randomoff(inp, off = 0.25):
    outlist = []
    for i in inp:
        value = random.random()
        outlist += [i+value*off]
    return outlist

def Spectra_Model():
    return 0

#%%
def GaussianRegression(DataFrame, Graph=False, classification = False):
    """
    

    Parameters
    ----------
    DataFrame : Pandas DataFrame 
        DataFrame of light curves that has been processed by "Data()"
    Graph : Boolean, optional
        Choose whether or not to plot the data. The default is False.
    Classification: Boolean, optional
        Choose wether or not to include classification data from DataFrame

    Returns
    -------
    GaussianFitted : Pandas DataFrame
        DataFrame of light curve data processed by Gaussian Regression 

    """
    sf = DataFrame
    SN_uqlst = sf.event.unique()
    for i in SN_uqlst:
        SNdf = sf[sf['event']==i]
        SNdf['mjd'] = SNdf['mjd'] - (SNdf['mjd'].min() -1)
        mjdrange = np.asarray([min(SNdf['mjd'].tolist()),max(SNdf['mjd'].tolist())])
        D = Multi_Band_GP(x_range = mjdrange, x=SNdf['mjd'].to_numpy(),
                          y=SNdf['magnitude'].to_numpy(), y_err=SNdf['error'].to_numpy(),
                          dim=band_to_wvlg(SNdf['band'].to_numpy()),
                          n_samples= 500)
        if Graph ==True:
            plt.errorbar(D[0],D[1],D[2], c=band_to_color(D[3]), alpha = 0.2, ls = 'None')
            plt.errorbar(x = SNdf['mjd'].to_numpy(), y = SNdf['magnitude'].to_numpy(), yerr = SNdf['error'].to_numpy(), ls = 'None')
            plt.xlabel("MJD")
            plt.ylabel("Magnitude")
            plt.title(SNdf['event'].unique())
            plt.gca().invert_yaxis()
            plt.show()
        GaussianFitted = pd.DataFrame()
        GaussianFitted['mjd'] = D[0]
        GaussianFitted['magnitude'] = D[1]
        GaussianFitted['error'] = D[2]
        GaussianFitted['band'] = D[3]
        y = pd.Series(data=np.zeros(1000)).astype(int)
        y = y.replace(0,i)
        GaussianFitted['event'] = y
        if classification == True:
            x = pd.Series(data = np.zeros(1000)).astype(int)
            x = x.replace(0, SNdf['class'].unique()[0])
            GaussianFitted['class'] = x
            
    return GaussianFitted

#%%
def Fit(Data, Bounds = None, Graph = False, TFI = 10, TRI = -2, classification = False):
    """
    

    Parameters
    ----------
    Data : Data Frame
        Data Frame of light curve data
    Bounds : list, optional
        Set custom bounds for each parameter for fitting. The default is None.
        Bound Syntax Example:
            Bounds = [[min=1, min = 0, min=x.min()-100, min=0, min = -10],[max=y.max()*10, max = np.inf, max=x.max(), max = 100, max=0]]

    Graph: True or Flase, optional
        Choose whether or not to display a graph of the fitting. The default is False
    TFI : float, optional
        Initial guess for fall time. The default is 10.
    TRI : float, optional
        Initial guess for rise time. The default is -2.
    
    Returns
    -------
    TYPE
        Data frame of function parameters

    """

    if Bounds == None:
        Function_Parameters = pd.DataFrame(columns = ["band", 'A', 'B', 't_int', 't_fall', 't_rise', 'class'])
        Data_uqlst = Data['event'].unique()
        for i in Data_uqlst:
            GaussianFitted = GaussianRegression(Data[Data['event']==i])
        
        
        
    
            def Lightcurve(x, A, B, t_int, t_fall, t_rise):
                """1-d Lightcurve: Lightcurve(x, A, B, t_int, t_fall, t_rise)"""
                return A*((np.exp(-1*(x - t_int)/t_fall))/(1 + np.exp((x - t_int)/t_rise))) + B
            Uqe_Bands = GaussianFitted['band'].unique()  #test all unique bands
            xt = []  #Temporary lists for creating accurate estimants for fitting
            yt = []
            FilOut = 0 #Filter value (increaces with more bands stimated to be unfitable)
            PurityFill = 0 #another filter value
            
            for Band in Uqe_Bands:
                if classification == True:
                    Class = GaussianFitted[GaussianFitted['band']==Band]['class']
                
                x = GaussianFitted[GaussianFitted['band']==Band]['mjd'].astype(float)  #Creates x for fitting estimants
                y = GaussianFitted[GaussianFitted['band']==Band]['magnitude'].astype(float) #Creates y for fitting
                yer = GaussianFitted[GaussianFitted['band']==Band]['error'].astype(float) #Creates y_err for fitting estimants
                yl = list(y.values.flatten()) #Creates all band length y list for estimants
                if yl.index(y.max())==0 or yl.index(y.max())==len(y):  #Checks that the highest magnitude in each band indicates an event has occoured (magnitude is hagher tan 50 at a single or multiple points)
                    FilOut+=1
                xt.append(x[y.idxmax()]) #appends fitting estimant list
                yt.append(y.max()) #appends fitting estimant list
            MX = xt[yt.index(max(yt))] #finds the highest of the fitting estimants time
            MMX = max(yt) #Sets filter value
    
                
            for UBand in Uqe_Bands: #Begin fitting one supernova
                #if GaussianFitted[GaussianFitted['band']==UBand].empty==False and len(GaussianFitted[GaussianFitted['band']==UBand]['mjd'])>=5: #Begin fitting one band
                    
                if classification == True:
                    Class = GaussianFitted[GaussianFitted['band']==UBand]['class']
                
                x = GaussianFitted[GaussianFitted['band']==UBand]['mjd'].astype(float) #Sets x for fitting
                y = GaussianFitted[GaussianFitted['band']==UBand]['magnitude'].astype(float) #Sets y for fitting
                y_err = GaussianFitted[GaussianFitted['band']==UBand]['error'].astype(float) #Sets y_err for fitting
                df2 = GaussianFitted[GaussianFitted['band']==UBand].copy() #Copys band for Chi^2 calculation
                lcmodel = lmfit.Model(Lightcurve, nan_policy='propagate') #Sets the fitting to the lightcurve shape and sets failed fits to return NAN rather than crash
                lcmodel.set_param_hint('A', min=-(y.max()*10), max=0, value=-abs(y.min()-y.max())) #Sets guess for amplitude
                lcmodel.set_param_hint('B', min=y.max(), value=y.max()+.5) #Sets guess for y-axis shift or baseline
                lcmodel.set_param_hint('t_int', min = x.min()-100, max=x.max(), value=MX-10) #Sets the time initial value which falls just before the maximum
                lcmodel.set_param_hint('t_fall', min= 0, value=TFI) #Sets the fall time
                lcmodel.set_param_hint('t_rise', max= 0, value=TRI) #Sets the rise time
                result = lcmodel.fit(y, x=x) #Outputs restults of the fit
                CST = result.best_values
                AA = CST['A'] #Sets the initial fit values to the weighted fit guesses
                BB = CST['B']
                TI = CST['t_int']
                TF = CST['t_fall']
                TR = CST['t_rise']
                        
                if Graph == True:
                    time=np.linspace(x.min(), x.max(), 100)
                    LC = (AA*((np.exp(-1*(time-TI)/TF))/(1+np.exp((time-TI)/TR)))+BB)
                    plt.plot(time, LC, 'r', label = '$\chi^{2}$ Fit')
                    plt.errorbar(GaussianFitted[GaussianFitted['band']==UBand]['mjd'].astype(float),
                                 GaussianFitted[GaussianFitted['band']==UBand]['magnitude'].astype(float),
                                 yerr= (GaussianFitted[GaussianFitted['band']==UBand]['error'].astype(float)),
                                 ls='none',
                                 ecolor='b',
                                 label = 'Magnitude Data',
                                 alpha = .2)
                    plt.axvline(x=TI, color =  'c', label = '$t_0$')
                    plt.axvline(x=TI+TF, color = 'r', label = '$t_{fall}$ (Relative to $t_0$)')
                    plt.axvline(x=TI+TR, color = 'b', label = '$t_{rise}$ (Relative to $t_0$)')
                    plt.title(i)
                    plt.gca().invert_yaxis()
                    plt.xlabel('Days')
                    plt.ylabel('Magnitude')
                    plt.legend()
                    plt.show()
                    plt.close()
                
                Features = pd.DataFrame(columns = ['band', 'A', 'B', 't_int', 't_fall', 't_rise'])
                
                Feat1 = pd.Series([AA])
                Feat2 = pd.Series([BB])
                Feat3 = pd.Series([TI])
                Feat4 = pd.Series([TF])
                Feat5 = pd.Series([TR])
                Feat6 = pd.Series([UBand])
                
                Features['A'] = Feat1
                Features['B'] = Feat2
                Features['t_int'] = Feat3
                Features['t_fall'] = Feat4
                Features['t_rise'] = Feat5
                Features['band'] = Feat6
                if classification == True:
                    Features['class'] = Class.unique()[0]
                
                Function_Parameters = pd.concat([Function_Parameters, Features], axis =0 )
                #if Function_Parameters['class'].unique()[0] == np.nan:
                    #Function_Parameters =Function_Parameters.drop(['class'],axis=1)
                
    else:
        Function_Parameters = pd.DataFrame(columns = ["band", 'A', 'B', 't_int', 't_fall', 't_rise', 'class'])
        Data_uqlst = Data['event'].unique()
        for i in Data_uqlst:
            GaussianFitted = GaussianRegression(Data[Data['event']==i])
        
            def Lightcurve(x, A, B, t_int, t_fall, t_rise):
                """1-d Lightcurve: Lightcurve(x, A, B, t_int, t_fall, t_rise)"""
                return A*((np.exp(-1*(x - t_int)/t_fall))/(1 + np.exp((x - t_int)/t_rise))) + B
            Uqe_Bands = GaussianFitted['band'].unique()  #test all unique bands
            xt = []  #Temporary lists for creating accurate estimants for fitting
            yt = []
            FilOut = 0 #Filter value (increaces with more bands stimated to be unfitable)
            PurityFill = 0 #another filter value
            
            for Band in Uqe_Bands:
                if classification == True:
                    Class = GaussianFitted[GaussianFitted['band']==Band]['class']
                x = GaussianFitted[GaussianFitted['band']==Band]['mjd'].astype(float)  #Creates x for fitting estimants
                y = GaussianFitted[GaussianFitted['band']==Band]['magnitude'].astype(float) #Creates y for fitting
                yer = GaussianFitted[GaussianFitted['band']==Band]['error'].astype(float) #Creates y_err for fitting estimants
                yl = list(y.values.flatten()) #Creates all band length y list for estimants
                if yl.index(y.max())==0 or yl.index(y.max())==len(y):  #Checks that the highest magnitude in each band indicates an event has occoured (magnitude is hagher tan 50 at a single or multiple points)
                    FilOut+=1
                xt.append(x[y.idxmax()]) #appends fitting estimant list
                yt.append(y.max()) #appends fitting estimant list
            MX = xt[yt.index(max(yt))] #finds the highest of the fitting estimants time
            MMX = max(yt) #Sets filter value
            for UBand in Uqe_Bands: #Begin fitting one supernova
                if GaussianFitted[GaussianFitted['band']==UBand].empty==False and len(GaussianFitted[GaussianFitted['band']==UBand]['mjd'])>=5: #Begin fitting one band
                    if classification == True:
                        Class = GaussianFitted[GaussianFitted['band']==UBand]['class']
                    
                    
                    x = GaussianFitted[GaussianFitted['band']==UBand]['mjd'].astype(float) #Sets x for fitting
                    y = GaussianFitted[GaussianFitted['band']==UBand]['magnitude'].astype(float) #Sets y for fitting
                    y_err = GaussianFitted[GaussianFitted['band']==UBand]['error'].astype(float) #Sets y_err for fitting
                    df2 = GaussianFitted[GaussianFitted['band']==UBand].copy() #Copys band for Chi^2 calculation
                    lcmodel = lmfit.Model(Lightcurve, nan_policy='propagate') #Sets the fitting to the lightcurve shape and sets failed fits to return NAN rather than crash
                    lcmodel.set_param_hint('A', min=Bounds[0][0], max=Bounds[0][1], value=abs(y.min()-y.max())) #Sets guess for amplitude
                    lcmodel.set_param_hint('B', min=Bounds[1][0], max =Bounds[1][1], value=0) #Sets guess for y-axis shift or baseline
                    lcmodel.set_param_hint('t_int', min = Bounds[2][0], max=Bounds[2][1], value=MX-10) #Sets the time initial value which falls just before the maximum
                    lcmodel.set_param_hint('t_fall', min= Bounds[3][0], max = Bounds[3][1], value=TFI) #Sets the fall time
                    lcmodel.set_param_hint('t_rise', min = Bounds[4][0] , max= Bounds[4][1], value=TRI) #Sets the rise time
                    result = lcmodel.fit(y, x=x) #Outputs restults of the fit
                    CST = result.best_values
                    
                    AA = CST['A'] #Sets the initial fit values to the weighted fit guesses
                    BB = CST['B']
                    TI = CST['t_int']
                    TF = CST['t_fall']
                    TR = CST['t_rise']
                if Graph == True:
                    time=np.linspace(x.min(), x.max(), 100)
                    LC = (AA*((np.exp(-1*(time-TI)/TF))/(1+np.exp((time-TI)/TR)))+BB)
                    plt.plot(time, LC, 'r', label = '$\chi^{2}$ Fit')
                    plt.errorbar(GaussianFitted[GaussianFitted['band']==UBand]['mjd'].astype(float),
                                 GaussianFitted[GaussianFitted['band']==UBand]['magnitude'].astype(float),
                                 yerr= (GaussianFitted[GaussianFitted['band']==UBand]['error'].astype(float)),
                                 ls='none',
                                 ecolor='b',
                                 label = 'Magnitude Data',
                                 alpha = .2)
                    plt.axvline(x=TI, color =  'c', label = '$t_0$')
                    plt.axvline(x=TI+TF, color = 'r', label = '$t_{fall}$ (Relative to $t_0$)')
                    plt.axvline(x=TI+TR, color = 'b', label = '$t_{rise}$ (Relative to $t_0$)')
                    plt.title(i)
                    plt.gca().invert_yaxis()
                    plt.xlabel('Days')
                    plt.ylabel('Magnitude')
                    plt.legend()
                    plt.show()
                    plt.close()
                
                Features = pd.DataFrame(columns = ['band', 'A', 'B', 't_int', 't_fall', 't_rise'])
                
                Feat1 = pd.Series([AA])
                Feat2 = pd.Series([BB])
                Feat3 = pd.Series([TI])
                Feat4 = pd.Series([TF])
                Feat5 = pd.Series([TR])
                Feat6 = pd.Series([UBand])
                
                Features['A'] = Feat1
                Features['B'] = Feat2
                Features['t_int'] = Feat3
                Features['t_fall'] = Feat4
                Features['t_rise'] = Feat5
                Features['band'] = Feat6
                if classification == True:
                    Features['class'] = Class.unique()[0]
                
                Function_Parameters = pd.concat([Function_Parameters, Features], axis =0 )
                #if Function_Parameters['class'].unique()[0] == np.nan:
                    #Function_Parameters =Function_Parameters.drop(['class'],axis=1)
                   
    return Function_Parameters.reset_index(drop=True)
#%%
def CorrPlot(Data):
    """
    

    Parameters
    ----------
    Data : Pandas DataFrame
        Fit Parameters from output of "Fit()"
        
    Generates a correlation plot of fitted light curve parameters.
    
    Returns
    -------
    None.

    """
    #Data['class'] = Data['class'].replace(['Type 1a', 'Type 1b', 'Type 1c', 'Type 2', 'Type 2b', 'Type 2n', 'Type 2p'], [0, 1, 2, 3, 4, 5, 6])
    import matplotlib as mpl
    from matplotlib import cm

    GBand = Data[Data['band']==357.0]
    RBand = Data[Data['band']==476.7]
    
    AvsAG = (GBand['A'])
    AvsAR = RBand['A']
    AvstfallG =(GBand['A'])/(GBand['t_fall'])
    AvstfallR = (RBand['A'])/(RBand['t_fall'])
    AvstriseG = (GBand['A'])/(GBand['t_rise'])
    AvstriseR = (RBand['A'])/(RBand['t_rise'])
    tfallvstriseG = (GBand['t_fall'])/(GBand['t_rise'])
    tfallvstriseR = (RBand['t_fall'])/(RBand['t_rise'])
    tfallG = GBand['t_fall']
    tfallR = RBand['t_fall']
    triseG = GBand['t_rise']
    triseR = RBand['t_rise']
    
    fig, axs = plt.subplots(2, 3, figsize=(25, 15))
       
    axs[0, 0].scatter(np.abs(AvsAR), np.abs(AvsAG), alpha=.1, c=GBand['labels'])
    axs[0, 1].scatter(np.abs(AvstfallR), np.abs(AvstfallG), alpha=.1,c=GBand['labels'])
    axs[0, 2].scatter(np.abs(AvstriseR), np.abs(AvstriseG), alpha=.1,c=GBand['labels'])
    axs[1, 0].scatter(np.abs(tfallvstriseR), np.abs(tfallvstriseG), alpha=.1,c=GBand['labels'])
    axs[1, 1].scatter(tfallR, tfallG, alpha=.1,c=GBand['labels'])
    axs[1, 2].scatter(np.abs(triseR), np.abs(triseG), alpha=.1,c=GBand['labels'])
    
    axs[0, 0].set_title('G Band A vs R Band A')
    axs[0, 1].set_title('G Band A vs R Band t fall')
    axs[0, 2].set_title('G Band A vs R Band t rise')
    axs[1, 0].set_title('G Band t fall vs R Band t rise')
    axs[1, 1].set_title('G Band tfall vs R Band tfall')
    axs[1, 2].set_title('G Band t rise vs R Band t rise')
    
    axs[0,0].set_xscale('log')
    axs[0,0].set_yscale('log')
    
    axs[0,1].set_xscale('log')
    axs[0,1].set_yscale('log')
    
    axs[0,2].set_xscale('log')
    axs[0,2].set_yscale('log')
    
    axs[1,0].set_xscale('log')
    axs[1,0].set_yscale('log')
    
    axs[1,1].set_xscale('log')
    axs[1,1].set_yscale('log')
    
    axs[1,2].set_xscale('log')
    axs[1,2].set_yscale('log')
    
    axs[0,0].set_xlim([10**-15,10**10])
    axs[0,1].set_xlim([10**-20,10**10])
    axs[0,2].set_xlim([10**-15,10**10])
    
    
    #handles, labels =  axs[0, 0].scatter(AvsAR*(-1), AvsAG*(-1), c = GBand['class']).legend_elements()
    
    #legend = fig.legend(handles , ['Type 1a: 21', 'Type 1b: 12', 'Type 1c: 12', 'Type 2: 21', 'Type 2b', 'Type 2n', 'Type 2p', 'Bogus Events'], loc="upper right", title="Classification", prop={'size': 15})
    
    
    fig.suptitle("Correlation Between Fitted Light Curve Parameters", fontsize = 25)
    
    plt.show()
#%%
#%%
#Some spectrally classified supernova to run tests on
Type1a = [
"ZTF20aahptds",
"ZTF19abimzvh",
"ZTF20acezpem",
"ZTF19abitbcj",
"ZTF20abccojn",
"ZTF20accsfkt",
"ZTF20acrdemq",
"ZTF19aavnwxf",
"ZTF20acgnrqu",
"ZTF19acyfakk",
"ZTF20aaljgcp",
"ZTF20acjhhqx", 
"ZTF20achzugy",
"ZTF20acikuon",
"ZTF20achvmoj",
"ZTF20acgnzhj",
"ZTF20achawym",
"ZTF20achdwmq",
"ZTF20achksbc",
"ZTF20achbekh",
"ZTF20acgnrqu"] #21 #type 1a

Type1b = [
"ZTF20acvebcu", #1b - 12
"ZTF20abvquuo",
"ZTF20abjpvce",
"ZTF18aakkrjm",
"ZTF20aalcyih",
"ZTF20aajcdad",
"ZTF19acmelor",
"ZTF19abztknu",
"ZTF19abqxibq",
"ZTF19abqqrgy",
"ZTF19aamgghn",
"ZTF18abktmfz"]

Type1c = [
"ZTF19acyogrm", #1c - 12
"ZTF19abupned", 
"ZTF21aaaubig",
"ZTF20acwofhd",
"ZTF20acueziy",
"ZTF20actpqgc",
"ZTF20acfqngt",
"ZTF20abxpoxd",
"ZTF20abvvnqh",
"ZTF19abfiqjg",
"ZTF19abdoior",
"ZTF19abcegvm"
]
Type2 = [
"ZTF21aadkhte", #type 2 -21
"ZTF21aabygea", 
"ZTF20acytfkf",
"ZTF20acyqzeu",
"ZTF20acwgxhk",
"ZTF20acvjlev",
"ZTF20acvezdt",
"ZTF20acuhgar",
"ZTF20acpkses",
"ZTF20acnzkxb",
"ZTF20acocohy",
"ZTF20acnezne",
"ZTF20acmaaan",
"ZTF19abgpgyp",
"ZTF19abdkgwo",
"ZTF19abbxykm",
"ZTF19aatqzim",
"ZTF18aaqkdwu",
"ZTF19aatlqdf",
"ZTF19aathllr",
"ZTF19aapafqd"
]
Type2b = [
"ZTF21aaabrpu", #type 2b - 12
"ZTF19aanbpus",
"ZTF19abxtcio",
"ZTF19abxjrge",
"ZTF19abqmsnk",
"ZTF19aazwnwy",
"ZTF19aapadxs",
"ZTF19aalzvnu",
"ZTF19aaknate",
"ZTF18acrcyqw",
"ZTF18acnmifq",
"ZTF18acbzvpg"
]
Type2n = [
"ZTF21aaabxbd", #type 2n - 12
"ZTF20aczgmai",
"ZTF20abonvte",
"ZTF20abpmqnr",
"ZTF20abjbgjj",
"ZTF20aayvmyh",
"ZTF20aadtarr",
"ZTF20aaaweke",
"ZTF20aacbyec",
"ZTF19acyjviz",
"ZTF19ablpasc",
"ZTF19ablojrw"
]
Type2p = [
"ZTF20abxmwwd", #2p - 8
"ZTF19acbhvgi",
"ZTF19aaycrgf",
"ZTF19aavjzrz",
"ZTF19aaduufr",
"ZTF18acyybvg",
"ZTF18abckutn",
"ZTF21aaagypx",
]

bogus = [
"ZTF18aaadcmc",
"ZTF19aauhwqv",
"ZTF18aahajkt",
"ZTF18aaadqsa",
"ZTF19aasceup",
"ZTF18aabdajj",
"ZTF17aabtrfz",
"ZTF20aajbvnb",
"ZTF19aanrsug",
"ZTF18acerlqm",
"ZTF20aakrfek",
"ZTF17aabueix",
"ZTF19aascegf",
"ZTF17aabxife",
"ZTF19aanwzip",
"ZTF19aaapejy",
"ZTF18aajjddu",
"ZTF19abhjvps",
"ZTF18acepfiu",
"ZTF18aarzkhb",
"ZTF19aaimcud",
"ZTF18aawqcsw",
"ZTF18aayiynz",
"ZTF18aavhyxv",
"ZTF18accwpyq",
    ]
