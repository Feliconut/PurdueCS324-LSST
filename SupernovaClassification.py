# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 20:15:51 2021

@author: blgnm
"""

import george
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import random
from astropy.stats import biweight_location
import numpy as np
import pandas as pd
from imblearn.over_sampling import SMOTE
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import sklearn
import matplotlib.pyplot as plt
import pywt
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from antares_client.search import search
from antares_client._api.models import Locus
from antares_client.search import get_by_ztf_object_id
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from imblearn.pipeline import Pipeline as imbpipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import BorderlineSMOTE 
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn import datasets, metrics, model_selection, svm
from sklearn.model_selection import RepeatedStratifiedKFold
from imblearn.ensemble import BalancedRandomForestClassifier



class ZTFData:
    def __init__(self,ZTF_id = None):
        self.ZTF_id = ZTF_id
        
    
    def get_id(self):
        return self.ZTF_id
    
    def get_raw_data(self):
        for i in self.ZTF_id:
            locus = get_by_ztf_object_id(i)
            return locus.lightcurve
    
    def get_lightcurve(self):
        """
        

        Returns
        -------
        Light Curve of all ztf objects stored within ZTFData()

        """
        for i in self.ZTF_id:
            plt.figure(num=None, figsize=(12, 5), dpi=100)
            locus = get_by_ztf_object_id(i)
            lc = locus.lightcurve
            lc['alert_type'] = lc['alert_id'].apply(lambda x: x[:x.index(':')])
            for (pb, at), df in lc.groupby(['ant_passband', 'alert_type']):
                is_candidate = at == 'ztf_candidate'
                plt.errorbar(
                    x=df['ant_mjd'],
                    y=df['ant_mag'] if is_candidate else df['ant_maglim'],
                    yerr=df['ant_magerr'],
                    #  uplims=(at!='ztf_candidate'),
                    label=pb + '  ' + at[4:],
                    color=pb,
                    fmt=('o' if is_candidate else '^') + pb.lower(),
                    alpha=1 if is_candidate else 0.3)
                plt.plot(df['ant_mjd'], df['ant_mag'])
            plt.title(str(i))
            plt.xlabel('MJD')
            plt.ylabel('Magnitude')
            plt.legend()
            plt.gca().invert_yaxis()
            plt.show()
          

    def Data(self, labels = None):
        """
        

        Parameters
        ----------
        labels : Pandas Series, optional
            Classification labels, if applicable (I don't remember if that part works tbh). The default is None.

        Returns
        -------
        Pandas data frame of processed event information for all provided ZTF ID's.

        """
       
        Supernovadf = pd.DataFrame(columns = ["mjd", "band", "magnitude", "error", "event","class"])
            #plt.title(locus.properties['ztf_object_id'])
        loop = tqdm(total = len(self.ZTF_id), position =0, leave = False)
        for i in self.ZTF_id:
            
            locus = get_by_ztf_object_id(i)
            try:
                Data = locus.lightcurve
            except Exception:
                pass
            
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
            
            GBand['class'] = num1
            RBand['class'] = num
            
            
            
            GBand['event'] = GBand['event'].replace([0], [str(i)])
            RBand['event'] = RBand['event'].replace([0], [str(i)])
            
        
            Band = pd.concat([GBand, RBand], axis = 0, ).reset_index(drop=True)
            Supernovadf = pd.concat([Supernovadf, Band], axis = 0).reset_index(drop=True)
            loop.set_description("Fetching Data...".format(len(i)))
            loop.update(1)
        
        loop.close()
            
        return(Supernovadf)

    def GaussianRegression(self, data = None, classification = True, DateRange = None, n_samples = 100):
        """
        

        Parameters
        ----------
        data : Pandas Data Frame, optional
            Pandas Data Frame with info from Data(); if None, it will pull data based off stored ZTF ID's. 
            The default is None.
        classification : Boolean, optional
            If you are making a training set to True (I always keep it True personally, not sure if it works otherwise). 
            The default is True.
        DateRange : Integer, optional
            How many days you want the classifier to look at. The default is None.
        n_samples : Integer, optional
            The number of samples GP Regression takes from the data. The default is 100.

        Returns
        -------
        Pandas Data Frame
            Pandas Data Frame of GP Regression Data.

        """
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
        
        if data is None:
            sf = self.Data()
        else:
            sf = data
        Gaus = pd.DataFrame()
        pd.options.mode.chained_assignment = None  # default='warn'
        SN_uqlst = sf.event.unique()
        loop = tqdm(total = len(SN_uqlst), position =0, leave = False)
        
        
        for i in SN_uqlst:
            SNdf = sf[sf['event']==i]
            SNdf['mjd'] = SNdf['mjd'] - (SNdf['mjd'].min() -1)
            if DateRange is not None:
                SNdf = SNdf[SNdf['mjd'] < DateRange]
            b = SNdf['band'].unique() == np.array([0.0, 1.0])
            if b[0] != True or b[1] != True:
                continue
            mjdrange = np.asarray([min(SNdf['mjd'].tolist()),max(SNdf['mjd'].tolist())])
            D = Multi_Band_GP(x_range = mjdrange, x=SNdf['mjd'].to_numpy(),
                              y=SNdf['magnitude'].to_numpy(), y_err=SNdf['error'].to_numpy(),
                              dim=band_to_wvlg(SNdf['band'].to_numpy()),
                              n_samples= n_samples)
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
            Gaus = pd.concat([Gaus, GaussianFitted])
            loop.set_description("Computing GPR...".format(len(i)))
            loop.update(1)
        loop.close()
        return Gaus
    
    def Gpr_Graph(self, DateRange = None ,n_samples = 100):
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
            Gaussian = self.GaussianRegression(DateRange = DateRange,n_samples = n_samples)
            for i in Gaussian['event'].unique():
                plt.errorbar(Gaussian[Gaussian['event']==i]['mjd'],Gaussian[Gaussian['event']==i]['magnitude'],Gaussian[Gaussian['event']==i]['error'], c=band_to_color(Gaussian[Gaussian['event']==i]['band']), alpha = 0.2, ls = 'None')
                
                #plt.errorbar(x = SNdf['mjd'].to_numpy(), y = SNdf['magnitude'].to_numpy(), yerr = SNdf['error'].to_numpy(), ls = 'None')
                plt.xlabel("MJD")
                plt.ylabel("Magnitude")
                plt.title(i)
                plt.gca().invert_yaxis()
                plt.show()
                plt.close()
                
    def Wavelet(self, Data = None, WaveletType = 'sym2', classification = True, Date = None, length = 150):
        """
        

        Parameters
        ----------
        Note: This version Processes both bands seperately, see Wavelet2() for multiband processing
        
        Data : Pandas Data Frame, optional
            Pandas DataFrame processed by Data(). The default is None.
        WaveletType : Str, optional
            Type of Wavelet transformation to be used. The default is 'sym2'.
        classification : Boolean, optional
            If you are making a training set to True (I always keep it True personally, not sure if it works otherwise). 
            The default is True.
        Date : Integer, optional
            How many days you want the classifier to look at. The default is None. The default is None.
        length : Integer, optional
            Set maximum event length; all events longer than set length are filtered out. The default is 150.

        Returns
        -------
        Function_Parameters : Pandas DataFrame
            Event Information such as ZTF ID and classification.
        Coefficients : Numpy Array
            List of Wavelet transformation Coefficients.

        """
       
        from tqdm import tqdm
        Function_Parameters = pd.DataFrame()
        Coefficients = list()
        
        if Data is None:    
            Data = self.Data()
        Gaussian = self.GaussianRegression(data = Data, DateRange = Date)
        
        Data_uqlst = Data['event'].unique()
        loop = tqdm(total = len(Data['event'].unique()), position =0, leave = False)
        

        for i in Data_uqlst:
            
            b = Data[(Data['event']==i)]['band'].unique() == np.array([0.0, 1.0])
            if b[0] != True or b[1] != True:
                continue
        
            if max(Data[Data['event']==i]['mjd'])-min(Data[Data['event']==i]['mjd']) > length:
                #print(len(Data[Data['event']==i]['mjd']))
                continue
            else:
                GaussianFitted = Gaussian[Gaussian['event']==i]
                Uqe_Bands = GaussianFitted['band'].unique()  
                
                for UBand in Uqe_Bands:
                    if classification == True:
                        Class = GaussianFitted[GaussianFitted['band']==UBand]['class']
                    
                    x = GaussianFitted[GaussianFitted['band']==UBand]['mjd'].astype(float)
                    y = GaussianFitted[GaussianFitted['band']==UBand]['magnitude'].astype(float)                    
                    y_err = GaussianFitted[GaussianFitted['band']==UBand]['error'].astype(float)                   
                    signal = y.values.squeeze()
    
                    ca = np.array(pywt.swt(np.array(signal), WaveletType, level = 2, axis = 0))
    
                    npoints=len(ca[0, 0, :])
                    coefficients =ca.reshape(2*2, npoints)
                        
                    Features = pd.DataFrame(data = {'band': [UBand], 'event': str(i), 
                                                    'delta':y.values.max()-y.values.min(), 'variance':y.var(), 
                                                    'duration': max(Data[Data['event']==i]['mjd'])-min(Data[Data['event']==i]['mjd'])})
                    if classification == True:
                        Features['class'] = Class.unique()[0]
                        
                    Coefficients.append(coefficients.flatten())
                    
                    Function_Parameters = pd.concat([Function_Parameters, Features], axis =0 )
                    Function_Parameters = Function_Parameters.reset_index(drop=True)
            loop.set_description("Computing Wavelet Transform...".format(len(i)))
            loop.update(1)
        loop.close()
        Function_Parameters['class'] = Function_Parameters['class'].replace(['SN Ia', 'SN II', 'SN Ib/c', 'SLSN'], [0,1,2,3])
        #Function_Parameters['class'] = Function_Parameters['class'].replace(['SN Ia', 'SN II', 'SN IIn', 'SN IIP', 'SN Ia-91T', 'SLSN-I', 'SLSN-II', 'SN Ic', 'SN Ib', 'SN Ic-BL', 'SN IIb', 'SN Ia-pec', 'SN Ibn', 'SN Ia-91bg'], [0,1,2,3,4,5,6,7,8,9, 10,11,12,13])

        return Function_Parameters, Coefficients

    
    def Wavelet2(self, Data = None, WaveletType = 'sym2', classification = True, Date = None, length = 150):
        """
        

        Parameters
        ----------
        Note: This version Processes both bands together, see Wavelet() for seperate band processing
        
        Data : Pandas Data Frame, optional
            Pandas DataFrame processed by Data(). The default is None.
        WaveletType : Str, optional
            Type of Wavelet transformation to be used. The default is 'sym2'.
        classification : Boolean, optional
            If you are making a training set to True (I always keep it True personally, not sure if it works otherwise). 
            The default is True.
        Date : Integer, optional
            How many days you want the classifier to look at. The default is None. The default is None.
        length : Integer, optional
            Set maximum event length; all events longer than set length are filtered out. The default is 150.

        Returns
        -------
        Function_Parameters : Pandas DataFrame
            Event Information such as ZTF ID and classification.
        Coefficients : Numpy Array
            List of Wavelet transformation Coefficients.

        """
       
        from tqdm import tqdm
        Function_Parameters = pd.DataFrame()
        Coefficients = list()
        
        if Data is None:    
            Data = self.Data()
        Gaussian = self.GaussianRegression(data = Data, DateRange = Date)
        
        Data_uqlst = Data['event'].unique()
        loop = tqdm(total = len(Data['event'].unique()), position =0, leave = False)
        

        for i in Data_uqlst:
            
            b = Data[(Data['event']==i)]['band'].unique() == np.array([0.0, 1.0])
            if b[0] != True or b[1] != True:
                continue
        
            if max(Data[Data['event']==i]['mjd'])-min(Data[Data['event']==i]['mjd']) > length:
                #print(len(Data[Data['event']==i]['mjd']))
                continue
            
            GaussianFitted = Gaussian[Gaussian['event']==i]
            
            
                
            if classification == True:
                Class = GaussianFitted['class']

            x = GaussianFitted['mjd'].astype(float)
            y = GaussianFitted['magnitude'].astype(float)                    
            y_err = GaussianFitted['error'].astype(float)                   
            signal = y.values.squeeze()
            if len(signal) == 0:
                continue
            from scipy import integrate
            Area = integrate.simpson(y, x)
            
            ca = np.array(pywt.swt(np.array(signal), WaveletType, level = 2, axis = 0))

            npoints=len(ca[0, 0, :])
            coefficients =ca.reshape(2*2, npoints)
            
            
            Features = pd.DataFrame(data = {'event': str(i), 
                                            'delta':y.values.max()-y.values.min(), 'variance':y.var(), 
                                            'duration': max(Data[Data['event']==i]['mjd'])-min(Data[Data['event']==i]['mjd']),
                                            'area':Area}, index=[0])
            if classification == True:
                Features['class'] = Class.unique()[0]

            Coefficients.append(coefficients.flatten())

            Function_Parameters = pd.concat([Function_Parameters, Features], axis =0 )
            Function_Parameters = Function_Parameters.reset_index(drop=True)
            loop.set_description("Computing Wavelet Transform...".format(len(i)))
            loop.update(1)
        loop.close()
        Function_Parameters['class'] = Function_Parameters['class'].replace(['SN Ia', 'SN II', 'SN Ib/c', 'SLSN'], [0,1,2,3])
        #Function_Parameters['class'] = Function_Parameters['class'].replace(['SN Ia', 'SN II', 'SN IIn', 'SN IIP', 'SN Ia-91T', 'SLSN-I', 'SLSN-II', 'SN Ic', 'SN Ib', 'SN Ic-BL', 'SN IIb', 'SN Ia-pec', 'SN Ibn', 'SN Ia-91bg'], [0,1,2,3,4,5,6,7,8,9, 10,11,12,13])

        return Function_Parameters, Coefficients
    
    def DimensionalityReduction2(self, Coefficients =None, labels=None, smot = False, n = 20, Trainset = True):
        """
        

        Parameters
        ----------
        Use this version if you used Wavelet2() (Multiband processing)
        
        Coefficients : Pandas Data Frame, optional
            Provide your own wavelet coefficients. The default is None.
        labels : Pandas Data Frame, optional
            Provide your own labels. The default is None.
        smot : Boolean, optional
            Choose whether or not to use SMOTE. The default is False.
        n : Integer, optional
            Output Dimension. The default is 20.
        Trainset : Boolean, optional
            Specify if this is the training set or if its unlabeled data. The default is True.

        Returns
        -------
        Pandas Data Frame
            Pandas Data Frame of PCA reduced Wavelet Coefficients.
        Function
            If Trainset = True, returns PCA() to use on unlabeled data.


        """
        if Coefficients is not None:    
            labels = labels
            Coefficients = Coefficients
        else:
            labels, Coefficients = self.Wavelet()
        
        
        Coefficients = pd.concat([pd.DataFrame(data=labels),pd.DataFrame(data=Coefficients)],axis=1)

        Coeff = Coefficients.iloc[:,6:]
        
        pca = PCA(n_components = n, whiten = True)
        if smot == True:
            sm = SMOTE()
            Coeff, labels= sm.fit_resample(Coeff, Coefficients['class'].ravel())
        
        print(Coeff)
        final = pca.fit_transform((Coeff))
        #RBand2, GBand2 = pd.DataFrame(data = {'Rdelta': RBand['delta'], 'Rvariance': RBand['variance']}), pd.DataFrame(data = {'Gdelta':GBand['delta'], 'Gvariance': GBand['variance']})
        if smot == True:
            events =pd.DataFrame(data = {'class': labels}).reset_index(drop=True)
        if smot == False:
            events =pd.DataFrame(data = {'event': Coefficients['event'], 'class': Coefficients['class']}).reset_index(drop=True)
        if Trainset == True:
            return pd.concat([events, pd.DataFrame(final)],axis=1), pca
        if Trainset == False:
            return pd.concat([events, pd.DataFrame(data = Coeff).reset_index(drop=True)],axis=1)

    
    def DimensionalityReduction(self, Coefficients =None, labels=None, smot = False, n = 20, Trainset = True):
        """
        

        Parameters
        ----------
        Use this version if you used Wavelet() (One band at a time processing)
        
        Coefficients : Pandas Data Frame, optional
            Provide your own wavelet coefficients. The default is None.
        labels : Pandas Data Frame, optional
            Provide your own labels. The default is None.
        smot : Boolean, optional
            Choose whether or not to use SMOTE. The default is False.
        n : Integer, optional
            Output Dimension. The default is 20.
        Trainset : Boolean, optional
            Specify if this is the training set or if its unlabeled data. The default is True.

        Returns
        -------
        Pandas Data Frame
            Pandas Data Frame of PCA reduced Wavelet Coefficients.
        Function
            If Trainset = True, returns PCA() to use on unlabeled data.

        """
        if Coefficients is not None:    
            labels = labels
            Coefficients = Coefficients
        else:
            labels, Coefficients = self.Wavelet()
        
        
        Coefficients = pd.concat([pd.DataFrame(data=labels),pd.DataFrame(data=Coefficients)],axis=1)

        GBand, RBand = Coefficients[Coefficients['band']==357.0].reset_index(drop=True), Coefficients[Coefficients['band']==476.7].reset_index(drop=True)
        print(RBand)
        pca = PCA(n_components = n, whiten = True)
        RBand = pd.concat([RBand.iloc[:,6:].reset_index(drop=True),GBand.iloc[:,6:].reset_index(drop=True)],axis=1, ignore_index=True)
        if smot == True:
            sm = SMOTE()
            RBand, labels= sm.fit_resample(RBand, GBand['class'].ravel())
        
        
        final = pca.fit_transform((RBand))
        #RBand2, GBand2 = pd.DataFrame(data = {'Rdelta': RBand['delta'], 'Rvariance': RBand['variance']}), pd.DataFrame(data = {'Gdelta':GBand['delta'], 'Gvariance': GBand['variance']})
        if smot == True:
            events =pd.DataFrame(data = {'class': labels}).reset_index(drop=True)
        if smot == False:
            events =pd.DataFrame(data = {'event': GBand['event'], 'class': GBand['class']}).reset_index(drop=True)
        if Trainset == True:
            return pd.concat([events, pd.DataFrame(final)],axis=1), pca
        if Trainset == False:
            return pd.concat([events, pd.DataFrame(data = RBand).reset_index(drop=True)],axis=1)



    def SupernovaTrainer(self, Train = None, y = None, **kwargs):
        """
        

        Parameters
        ----------
        Trains for Supernova vs Bogus Classification
        
        Train : Pandas Data Frame, optional
            Provide your own wavelet coefficients. The default is None.
        y : Pandas Data Frame, optional
            Provide data labels. The default is None.
        **kwargs : Classifier arguments
            Input arguments for classifier.

        Returns
        -------
        Function
            Trained Classifier.

        """
        if Train is None:
            Data = self.DimensionalityReduction()
            Train = Data.iloc[:,2:].reset_index(drop=True)
            y = Data['class'].reset_index(drop=True)
        svc = RandomForestClassifier(random_state=iterations, n_estimators = 30, min_samples_split = 6)
        
        if kwargs:
            classifier=AdaBoostClassifier(**kwargs)
        else:
            classifier=AdaBoostClassifier(base_estimator=svc,n_estimators=30, learning_rate =4)
        
        if evaluate == True:
            pipeline = imbpipeline(steps = [['classifier', classifier]])
            stratified_kfold = StratifiedKFold(n_splits=3, shuffle=True)
            print(cross_validate(pipeline, np.array(Train),np.array(y).ravel(),scoring = 'accuracy'))
            y_pred = cross_val_predict(pipeline, np.array(Train),np.array(y), cv = stratified_kfold)
            conf_mat = confusion_matrix(y, y_pred)
            plot_confusion_matrix1(conf_mat, ['SN', 'Bogus'], cmap = 'Blues')
            
        return classifier.fit(Train, y)
        
    def SupernovaTypeClassifierTrainer(self, Train, y, evaluate = True , smot = True,
                                       Ada = True, KNN = False,roc = True, Rand = False,
                                       grid = False, n = 1, fold = 3,n_components = 20, 
                                       metric = 'accuracy', param_grid = None,**kwargs):
        """
        

        Parameters
        ----------
        Trains for Supernova Type Classification
        
        Train : Pandas Data Frame, optional
            Provide your own wavelet coefficients. The default is None.
        y : Pandas Data Frame, optional
            Provide data labels. The default is None.
        evaluate : Boolean, optional
            Choose whether or not to show model performance. The default is True.
        Ada : Boolean, optional
            Choose to use Ada Boosted Random Forrest. The default is True.
        KNN : Boolean, optional
            Choose to use K-nearest neighbors. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            If you set both KNN and Ada to false, raises an error.

        Returns
        -------
        Function
            Trained Classifier.

        """
        if Train is not None:
            TrainingData, u = pd.concat([pd.DataFrame(data=y),pd.DataFrame(data=Train)],axis=1).reset_index(drop=True), y
        #else *** Remember to make this load in default training data
        svc = RandomForestClassifier(n_estimators = 30, min_samples_split = 6)
        TrainingData = TrainingData.sample(frac = 1).reset_index(drop=True)
        if kwargs:
            if Ada ==True:
                classifier = AdaBoostClassifier(**kwargs)
            if KNN == True:
                classifier = KNeighborsClassifier(**kwargs)
            if Rand == True:
                #classifier = RandomForestClassifier(**kwargs)
                classifier = BalancedRandomForestClassifier(**kwargs)
              
        else:
            classifier=AdaBoostClassifier(base_estimator=svc,n_estimators=30, learning_rate =2)
        #classifier = KNeighborsClassifier(n_neighbors=1500)

        if evaluate == True:
            
            if smot == True:
                pipeline = imbpipeline(steps = [['scale',MinMaxScaler()],['smote', SMOTE()],['classifier', classifier]])
            if smot == False:
                pipeline = imbpipeline(steps = [['scale',MinMaxScaler()],['classifier', classifier]])
            stratified_kfold = StratifiedKFold(n_splits=fold, shuffle=True)
            repeatstratified_kfold = RepeatedStratifiedKFold(n_splits=fold, n_repeats=n)
            cross = cross_validate(pipeline, np.array(TrainingData.iloc[:,1:]),np.array(TrainingData.iloc[:,0]),scoring = metric, cv = repeatstratified_kfold, n_jobs = -1)
            print(f'The mean {metric} over {fold} fold stratified crossvalidation repeated {n} times is {np.mean(cross["test_score"])}, with a standard deviation of {np.std(cross["test_score"])}')
            y_pred = cross_val_predict(pipeline, np.array(TrainingData.iloc[:,1:]),np.array(TrainingData.iloc[:,0]), cv = stratified_kfold, n_jobs = -1)
                        
            #conf_mat = confusion_matrix(y, y_pred)
            conf_mat = confusion_matrix(TrainingData.iloc[:,0], y_pred)
            disp = ConfusionMatrixDisplay(confusion_matrix=conf_mat)
            disp.plot(cmap = 'Blues')
            
            if grid == True:
                clf = GridSearchCV(pipeline, param_grid, n_jobs = -1, cv = stratified_kfold, scoring = 'f1_micro', verbose = 1)
                clf.fit(TrainingData.iloc[:,1:], TrainingData.iloc[:,0])
                
            
               
            
            plot_confusion_matrix1(conf_mat, ['Type 1a','Type 2', 'Type 1b/c', 'SLSN'], cmap = 'Blues')
  
        Classifier = pipeline.fit(TrainingData.iloc[:,1:], TrainingData.iloc[:,0])
        
        if grid == False:
            return Classifier
        if grid == True:
            return clf
    

def plot_confusion_matrix1(cm,
                          target_names,
                          title='Confusion matrix',
                          cmap=None,
                          normalize=True):
    """
    given a sklearn confusion matrix (cm), make a nice plot

    Arguments
    ---------
    cm:           confusion matrix from sklearn.metrics.confusion_matrix

    target_names: given classification classes such as [0, 1, 2]
                  the class names, for example: ['high', 'medium', 'low']

    title:        the text to display at the top of the matrix

    cmap:         the gradient of the values displayed from matplotlib.pyplot.cm
                  see http://matplotlib.org/examples/color/colormaps_reference.html
                  plt.get_cmap('jet') or plt.cm.Blues

    normalize:    If False, plot the raw numbers
                  If True, plot the proportions

    Usage
    -----
    plot_confusion_matrix(cm           = cm,                  # confusion matrix created by
                                                              # sklearn.metrics.confusion_matrix
                          normalize    = True,                # show proportions
                          target_names = y_labels_vals,       # list of names of the classes
                          title        = best_estimator_name) # title of graph

    Citiation
    ---------
    http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html

    """
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools

    accuracy = np.trace(cm) / np.sum(cm).astype('float')
    misclass = 1 - accuracy

    if cmap is None:
        cmap = plt.get_cmap('Blues')

    plt.figure(figsize=(8, 6))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()

    if target_names is not None:
        tick_marks = np.arange(len(target_names))
        plt.xticks(tick_marks, target_names, rotation=45)
        plt.yticks(tick_marks, target_names)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]


    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        if normalize:
            plt.text(j, i, "{:0.4f}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="black" if cm[i, j] > thresh else "black")
        else:
            plt.text(j, i, "{:,}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="black" if cm[i, j] > thresh else "black")


    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass))
    plt.show()
