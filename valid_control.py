
import pylab as plt
import numpy as np
import matplotlib.pylab as plt
import scipy.stats

from scipy.optimize import curve_fit
#from ..kdtree import KDTree
import pickle
import random
import math


sourceFluxField='base_PsfFlux'

color = {'all': 'grey', 'bright': 'blue',
         'iqr': 'green', 'rms': 'red'}

srdSpec = {
    'levels':("design", "minimum", "stretch"),
    'PA1':{"design": 5, "minimum": 8, "stretch": 3}, 'pa1Units':'mmag',
    'PF1':{"design": 10, "minimum": 20, "stretch": 5}, 'pf1Units':'%',
    'PA2':{"design": 15, "minimum": 15, "stretch": 10}, 'pa2Units':'mmag',
    'D1':5, 'd1Units':'arcmin',
    'AM1':{"design": 10, "minimum": 20, "stretch": 5}, 'am1Units':'mas',
    'AF1':{"design": 10, "minimum": 20, "stretch": 5}, 'af1Units':'%',
    'AD1':{"design": 20, "minimum": 40, "stretch": 10},' ad1Units':'mas',
    'D2':20, 'd2Units':'arcmin',
    'AM2':{"design": 10, "minimum": 20, "stretch": 5}, 'am2Units':'mas',
    'AF2':{"design": 10, "minimum": 20, "stretch": 5}, 'af2Units':'%',
    'AD2':{"design": 20, "minimum": 40, "stretch": 10}, 'ad2Units':'mas',
    'D3':200, 'd3Units':'arcmin',
    'AM3':{"design": 15, "minimum": 30, "stretch": 10}, 'am3Units':'mas',
    'AF3':{"design": 10, "minimum": 20, "stretch": 5}, 'af3Units':'%',
    'AD3':{"design": 30, "minimum": 50, "stretch": 20}, 'ad3Units':'mas',
}



radToDeg = 180./np.pi
degToArcs = 3600.
radToArcs = radToDeg * degToArcs
radToMas = radToArcs *1000

def save_pkl(obj, name):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pkl(name):
    with open( name, 'rb') as f:
        return pickle.load(f)


class Comparaison_multiple:
#    """ '/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation_phVisits_peudevisites/sources_ndarray_grp_4visits.pkl'

    def __init__(self, liste=[4,5,6,20], dico=False):#liste=[4,5,6,10,20,30,50,60]):
        self.liste=liste
        if dico:
            self.SRCDICT={}
        for i in liste:
            num=str(i)
            print num
            exec("sourcesArray"+num+ " = load_pkl('/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation_phVisits_peudevisites/sources_ndarray_grp_"+num+"visits.pkl')")
       
            exec("self.RA"+num+" = sourcesArray"+num+ "['coord_ra']")
            exec("self.Dec"+num+" = sourcesArray"+num+ "['coord_dec']")
            if dico:
                self.SRCDICT[i]={}
                self.SRCDICT[i]['RA']=eval('np.array(self.RA'+num+')*radToDeg')
                self.SRCDICT[i]['Dec']=eval('np.array(self.Dec'+num+')*radToDeg')


    def plotPos(self):
        colors = ['m','c','y','r','w','b','g', 'k']
        taille_points = []
        print "self.liste",self.liste
        for i in range(len(self.liste)):
            taille_points.append(2**(len(self.liste)+2-i))

        compt = 0
        for i in self.liste:
            num=str(i)
            ra = eval('np.array(self.RA'+num+')*radToDeg')
            dec = eval('np.array(self.Dec'+num+')*radToDeg')
            plt.scatter(ra, dec, c=colors[compt%len(colors)], s=taille_points[compt], label= num+ ' visites')
            plt.xlabel('RA (deg)')
            plt.ylabel('Dec (deg)')
            plt.title('position des sources en fonction du nombre de visites (couleurs)')
            plt.legend()
            plt.xlim(min(ra), max(ra))
            plt.ylim(min(dec), max(dec))
            compt+=1






# fonctions prealables aux plots :

def plotOutlinedLinesHorizontal(ax, *args, **kwargs):
    """Plot horizontal lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    plotOutlinedLines(ax.axhline, *args, **kwargs)


def plotOutlinedLinesVertical(ax, *args, **kwargs):
    """Plot vertical lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    plotOutlinedLines(ax.axvline, *args, **kwargs)


def plotOutlinedLines(ax_plot, x1, x2, x1_color=color['all'], x2_color=color['bright']):
    """Plot horizontal lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    ax_plot(x1, color='white', linewidth=4)
    ax_plot(x2, color='white', linewidth=4)
    ax_plot(x1, color=x1_color, linewidth=3)
    ax_plot(x2, color=x2_color, linewidth=3)


def plotAstromErrModelFit(snr, dist, fit_params=None,
                          color='red', ax=None, verbose=True):
    """Plot model of photometric error from LSST Overview paper
    http://arxiv.org/abs/0805.2366v4

    Astrometric Errors
    error = C * theta / SNR

    Parameters
    ----------
    snr : list or numpy.array
        S/N of photometric measurements
    dist : list or numpy.array
        Separation from reference [mas]
    """
    if ax is None:
        ax = plt.figure()
        xlim = [10, 30]
    else:
        xlim = ax.get_xlim()

    if fit_params is None:
        fit_params = fitAstromErrModel(snr, dist)

    x_model = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), num=100)
    fit_model_mas_err = astromErrModel(x_model, **fit_params)
    label = r'$C, \theta, \sigma_{\rm sys}$ =' + '\n' + \
            '{C:.2g}, {theta:.4g}, {sigmaSys:.4g} [mas]'.format(**fit_params)

    if verbose:
        print(fit_params)
        print(label)

    ax.plot(x_model, fit_model_mas_err,
            color=color, linewidth=2,
            label=label)
    # Set the x limits back to their original values.
    ax.set_xlim(xlim)


def fitAstromErrModel(snr, dist): # from check.py
    """Fit model of astrometric error from LSST Overview paper

    Parameters
    ----------
    snr : list or numpy.array
        Signal-to-noise ratio of photometric observations
    dist : list or numpy.array
        Scatter in measured positions [mas]

    Returns
    -------
    dict
        The fit results for C, theta, sigmaSys along with their Units.
    """
    fit_params, fit_param_covariance = \
        curve_fit(astromErrModel, snr, dist, p0=[1, 0.01])

    params = {'C': 1, 'theta': fit_params[0], 'sigmaSys': fit_params[1],
              'cUnits': '', 'thetaUnits': 'mas', 'sigmaSysUnits': 'mas'}
    return params

def astromErrModel(snr, theta=1000, sigmaSys=10, C=1, **kwargs): # from check.py
    """Calculate expected astrometric uncertainty based on SNR.

    mas = C*theta/SNR + sigmaSys

    Parameters
    ----------
    snr : list or numpy.array
        S/N of photometric measurements
    theta : float or numpy.array, optional
        Seeing
    sigmaSys : float
        Systematic error floor
    C : float
        Scaling factor

    theta and sigmaSys must be given in the same units.
    Typically choices might be any of arcsec, milli-arcsec, or radians
    The default values are reasonable astronominal values in milliarcsec.
    But the only thing that matters is that they're the same.

    Returns
    -------
    np.array
        Expected astrometric uncertainty.
        Units will be those of theta + sigmaSys.
    """
    return C*theta/snr + sigmaSys

#plots supplementaires  necessaires a plotPhotometry
def plotPhotErrModelFit(mag, mmag_err, fit_params=None, color='red', ax=None, verbose=True):
    """Plot model of photometric error from LSST Overview paper (Eq. 4 & 5)

    Parameters
    ----------
    mag : list or numpy.array
        Magnitude
    mmag_err : list or numpy.array
        Magnitude uncertainty or variation in *mmag*.
    fit_params : list or numpy.array
        Fit parameters to display
    ax : matplotlib.Axis, optional
        The Axis object to plot to.
    verbose : bool, optional
        Produce extra output to STDOUT
    """

    if ax is None:
        ax = plt.figure()
        xlim = [10, 30]
    else:
        xlim = ax.get_xlim()

    if fit_params is None:
        fit_params = fitPhotErrModel(mag, mmag_err)

    x_model = np.linspace(*xlim, num=100)
    fit_model_mag_err = photErrModel(x_model, **fit_params)
    fit_model_mmag_err = 1000*fit_model_mag_err
    labelFormatStr = r'$\sigma_{{\rm sys}} {{\rm [mmag]}}$, $\gamma$, $m_5 {{\rm [mag]}}$=' + '\n' + \
                     r'{sigmaSysMmag:6.4f}, {gamma:6.4f}, {m5:6.3f}'
    label = labelFormatStr.format(sigmaSysMmag=1000*fit_params['sigmaSys'],
                                  **fit_params)

    if verbose:
        print(fit_params)
        print(label)

    ax.plot(x_model, fit_model_mmag_err,
            color=color, linewidth=2,
            label=label)

    return fit_params


def fitPhotErrModel(mag, mmag_err): # from check.py
    """Fit model of photometric error from LSST Overview paper

    Parameters
    ----------
    mag : list or numpy.array
        Magnitude
    mmag_err : list or numpy.array
        Magnitude uncertainty or variation in *mmag*.

    Returns
    -------
    dict
        The fit results for sigmaSys, gamma, and m5 along with their Units.
    """
    mag_err = mmag_err / 1000
    fit_params, fit_param_covariance = \
        curve_fit(photErrModel, mag, mag_err, p0=[0.01, 0.039, 24.35])

    params = {'sigmaSys': fit_params[0], 'gamma': fit_params[1], 'm5': fit_params[2],
              'sigmaSysUnits': 'mmag', 'gammaUnits': '', 'm5Units': 'mag'}
    return params


def photErrModel(mag, sigmaSys, gamma, m5, **kwargs):# from check.py
    """Fit model of photometric error from LSST Overview paper
    http://arxiv.org/abs/0805.2366v4

    Photometric errors described by
    Eq. 4
    sigma_1^2 = sigma_sys^2 + sigma_rand^2

    Eq. 5
    sigma_rand^2 = (0.04 - gamma) * x + gamma * x^2 [mag^2]
    where x = 10**(0.4*(m-m_5))

    Parameters
    ----------
    mag : list or numpy.array
        Magnitude
    sigmaSq : float
        Limiting systematics floor [mag]
    gamma : float
        proxy for sky brightness and readout noise
    m5 : float
        5-sigma depth [mag]

    Returns
    -------
    numpy.array
        Result of noise estimation function
    """
    x = 10**(0.4*(mag - m5))
    sigmaRandSq = (0.04 - gamma) * x + gamma * x**2
    sigmaSq = sigmaSys**2 + sigmaRandSq
    return np.sqrt(sigmaSq)




class LoadDataValidation:
    def __init__(self, file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/sources_ndarray_grp_16visits.pkl'):
        self.sources = load_pkl(file_path)


    def createVariablesForPlots(self, additional=False, Extended=False ):
        RA_grp = []
        Dec_grp = []
        mag_grp = []
        snr_grp = []
        magerr_grp = []
        psf_fwhm_grp = []
        if Extended:
            extendedness_grp = []

        self.RA_mean= []
        self.Dec_mean = []
        self.posRMS = []
        self.mag_mean = []
        self.snr_med = []
        self.mag_rms = []
        self.magerr_med = []
        if additional:
            self.deltaRAcosdecs = []
            self.deltaDecs  = []
            self.deltamags = []
            self.RAcosDec_mean = []
            self.RAcosDec_RMS = []
            self.Dec_RMS = []
            self.Psf_fwhm_mean = []
            self.medsnrlong = []
       
        if Extended:
            self.extendedness = []
            self.maxextendedness = []

          #  extended = np.max(cat.get(extendedKey))
 
         #   psfSnr >= safeSnr and extended < safeMaxExtended
            

        #test = []
        #testi = []
        #for i, j in enumerate(self.sources['Nb_group']):
        for i in range(len(self.sources['Nb_group'])):
           # print 'i, j', i #,j
            if i == 0:
                val = 0
            else:
                val = self.sources['Nb_group'][i-1]

            if self.sources['Nb_group'][i] == val:
                RA_grp.append(self.sources['coord_ra'][i])
                Dec_grp.append(self.sources['coord_dec'][i])
                #test.append(self.sources['Nb_group'][i])
                mag_grp.append(self.sources[sourceFluxField+'_mag'][i])
                magerr_grp.append(self.sources[sourceFluxField+'_magerr'][i])
                snr_grp.append(self.sources[sourceFluxField+'_snr'][i])
                psf_fwhm_grp.append(self.sources['PSF-FWHM'][i])
                if Extended:
                    extendedness_grp.append(self.sources['base_ClassificationExtendedness_value'][i])
                #testi.append(i)
            else:
                # variables dans self.sources : ['visit','ccd','coord_ra','coord_dec','MJD-OBS',sourceFluxField+'_snr',sourceFluxField+'_flux',sourceFluxField+'_fluxSigma',sourceFluxField+'_mag',sourceFluxField+'_magerr','Nb_group','id','PSF-FWHM','FLUXMAG0','FLUXMAG0ERR', 'MeanGrpRa', 'MeanGrpDec', 'MedGrpSnr']
                self.RA_mean.append(np.mean(RA_grp))
                self.Dec_mean.append(np.mean(Dec_grp))
              #  print'ra meab',  np.mean(RA_grp),self.sources['MeanGrpRa'][i-1]

                self.posRMS.append(np.sqrt((np.std(RA_grp)* np.cos(np.mean(Dec_grp)))**2 + (np.std(Dec_grp))**2))
                self.mag_mean.append(np.mean(mag_grp))
                self.snr_med.append(np.median(snr_grp))
                self.mag_rms.append(np.std(mag_grp))
                self.magerr_med.append(np.median(magerr_grp))

                if additional:
                    self.medsnrlong += [np.median(snr_grp)]*len(RA_grp)
                    self.deltaRAcosdecs += list((np.array(RA_grp)-np.mean(RA_grp))* np.cos(np.mean(Dec_grp)))
                    self.deltaDecs += list(np.array(Dec_grp)-np.mean(Dec_grp))
                    self.deltamags += list(np.array(mag_grp)-np.mean(mag_grp))
                    self.RAcosDec_mean.append(np.mean(RA_grp) * np.cos(np.mean(Dec_grp)))
                    self.RAcosDec_RMS.append(np.std(RA_grp) * np.cos(np.mean(Dec_grp)))
                    self.Dec_RMS.append(np.std(Dec_grp))
                    self.Psf_fwhm_mean.append(np.mean( psf_fwhm_grp))
                if Extended:
                    self.extendedness += list(extendedness_grp)
                    self.maxextendedness += [ max(extendedness_grp)]*len(RA_grp)
            #    print 'test', test
            #    print 'testi', testi

                # initialisation du groupe suivant
                RA_grp = []
                Dec_grp = []
                
                mag_grp = []
                magerr_grp = []
                snr_grp = []
                Psf_fwhm_grp = []
                if Extended:
                    extendedness_grp = []

               # test = []
               # testi = []
                RA_grp.append(self.sources['coord_ra'][i])
                Dec_grp.append(self.sources['coord_dec'][i])
               # test.append(self.sources['Nb_group'][i])
                mag_grp.append(self.sources[sourceFluxField+'_mag'][i])
                magerr_grp.append(self.sources[sourceFluxField+'_magerr'][i])
                snr_grp.append(self.sources[sourceFluxField+'_snr'][i])
                psf_fwhm_grp.append(self.sources['PSF-FWHM'][i])
                if Extended:
                    extendedness_grp.append(self.sources['base_ClassificationExtendedness_value'][i])
             #   self.sources['Nb_group']

      
        self.RA_mean.append(np.mean(RA_grp))
        self.Dec_mean.append(np.mean(Dec_grp))
              #  print'ra meab',  np.mean(RA_grp),self.sources['MeanGrpRa'][i-1]

        self.posRMS.append(np.sqrt((np.std(RA_grp)* np.cos(np.mean(Dec_grp)))**2 + (np.std(Dec_grp))**2))
        self.mag_mean.append(np.mean(mag_grp))
        self.snr_med.append(np.median(snr_grp))
       
        self.mag_rms.append(np.std(mag_grp))
        self.magerr_med.append(np.median(magerr_grp))

        if additional:
            self.medsnrlong += [np.median(snr_grp)]*len(RA_grp)
            self.deltaRAcosdecs += list((np.array(RA_grp)-np.mean(RA_grp))* np.cos(np.mean(Dec_grp)))
            self.deltaDecs += list(np.array(Dec_grp)-np.mean(Dec_grp))
            self.deltamags += list(np.array(mag_grp)-np.mean(mag_grp))
            self.RAcosDec_mean.append(np.mean(RA_grp) * np.cos(np.mean(Dec_grp)))
            self.RAcosDec_RMS.append(np.std(RA_grp) * np.cos(np.mean(Dec_grp)))
            self.Dec_RMS.append(np.std(Dec_grp))
            self.Psf_fwhm_mean.append(np.mean( psf_fwhm_grp))
        if Extended:
            self.extendedness += list(extendedness_grp)
            self.maxextendedness += [ max(extendedness_grp)]*len(RA_grp)

        self.RA_mean = np.array(self.RA_mean)
        self.Dec_mean = np.array(self.Dec_mean)
        self.posRMS =  np.array(self.posRMS)
        self.mag_mean = np.array(self.mag_mean)
        self.snr_med = np.array(self.snr_med)
        self.mag_rms = np.array(self.mag_rms)
        self.magerr_med = np.array(self.magerr_med)

        if additional:
            self.deltaRAcosdecs = np.array(self.deltaRAcosdecs)
            self.deltaDecs = np.array(self.deltaDecs)
            self.deltamags = np.array(self.deltamags)
            self.RAcosDec_mean = np.array(self.RAcosDec_mean)
            self.RAcosDec_RMS = np.array(self.RAcosDec_RMS)
            self.Dec_RMS = np.array(self.Dec_RMS)
            self.Psf_fwhm_mean = np.array(self.Psf_fwhm_mean)
            self.medsnrlong = np.array(self.medsnrlong)

        if Extended:
            self.extendedness = np.array(self.extendedness)
            self.maxextendedness = np.array(self.maxextendedness)


class Validation(LoadDataValidation):

    def __init__(self, LOAD_val, additional=True, brightSnr=100., safeMaxExtended = 1.0):
        # variables dans self.sources : ['visit','ccd','coord_ra','coord_dec','MJD-OBS',sourceFluxField+'_snr',sourceFluxField+'_flux',sourceFluxField+'_fluxSigma',sourceFluxField+'_mag',sourceFluxField+'_magerr','Nb_group','id','PSF-FWHM','FLUXMAG0','FLUXMAG0ERR', 'MeanGrpRa', 'MeanGrpDec', 'MedGrpSnr']
        print "charge data"
        self.sources=LOAD_val.sources
        self.RA_mean = np.array(LOAD_val.RA_mean)
        self.Dec_mean = np.array(LOAD_val.Dec_mean)
        self.posRMS =  np.array(LOAD_val.posRMS)
        self.mag_mean = np.array(LOAD_val.mag_mean)
        self.snr_med = np.array(LOAD_val.snr_med)
        self.mag_rms = np.array(LOAD_val.mag_rms)
        self.magerr_med = np.array(LOAD_val.magerr_med)

        if additional:
            self.deltaRAcosdecs = np.array(LOAD_val.deltaRAcosdecs)
            self.deltaDecs = np.array(LOAD_val.deltaDecs)
            self.deltamags = np.array(LOAD_val.deltamags)
            self.RAcosDec_mean = np.array(LOAD_val.RAcosDec_mean)
            self.RAcosDec_RMS = np.array(LOAD_val.RAcosDec_RMS)
            self.Dec_RMS = np.array(LOAD_val.Dec_RMS)
            self.Psf_fwhm_mean = np.array(LOAD_val.Psf_fwhm_mean)
            self.medsnrlong = np.array(LOAD_val.medsnrlong)

        self.brightgrp, = np.where(self.medsnrlong>= brightSnr)
        self.brightsources = self.sources[self.brightgrp]

        if safeMaxExtended:
            self.extendedness = LOAD_val.extendedness
            self.maxextendedness = LOAD_val.maxextendedness 
            self.extendedgrp, = np.where(self.maxextendedness[self.brightgrp] < safeMaxExtended)
            self.brightsources = self.brightsources[ self.extendedgrp]


    def afficher(self):
        print 'validation'
    

    def calcPA1(self, numRandomShuffles=50):
        iqrPA1=[]
        for i in range(numRandomShuffles):

            self.doCalcPA1( verbose=False)
            print 'i, self.rmsSigma,self.iqrSigma',i,  self.rmsSigma, self.iqrSigma
            iqrPA1.append(self.iqrSigma)

        self.PA1 = np.mean(iqrPA1)
        print 'PA1=', self.PA1 # PA1 of validate_drp


    def doCalcPA1(self, verbose=False):

        self.diffmags = []
        self.magmean = []

        idgrps=set(self.brightsources['Nb_group'])

        for groupid in idgrps:
            inside_grp,= np.where(self.brightsources['Nb_group'] == groupid)
            mags_grp = self.brightsources[sourceFluxField+'_mag'][inside_grp]
            np.random.shuffle(mags_grp)
            diffmag = mags_grp[0] - mags_grp[1]
            self.diffmags.append((1000/math.sqrt(2)) * diffmag)
            self.magmean.append(np.mean( mags_grp))

        self.diffmags = np.array(self.diffmags)
        self.rmsSigma = math.sqrt(np.mean(self.diffmags**2))
        self.iqrSigma = np.subtract.reduce(np.percentile(self.diffmags, [75, 25])) / (scipy.stats.norm.ppf(0.75)*2)

        self.dmagsRMS = np.array(self.diffmags)




    def plotPA1(self,outputPrefix=""):
        """Plot the results of calculating the LSST SRC requirement PA1.
        
        Creates a file containing the plot with a filename beginning with `outputPrefix`.
        
        Parameters
        ----------
        pa1 : pipeBase.Struct
        Must contain:
        rms, iqr, magMean, magDiffs
        rmsUnits, iqrUnits, magDiffsUnits
        outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
        'Cfht_output_r_AM1_D_5_arcmin_17.0-21.5.png'
        for an AMx.name=='AM1' and AMx.magRange==[17, 21.5]
        """
        diffRange = (-100, +100)
        
        fig = plt.figure(figsize=(18, 12))
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.scatter(self.magmean, self.diffmags, s=10, color=color['bright'], linewidth=0)
        ax1.axhline(+self.rmsSigma , color=color['rms'], linewidth=3)
        ax1.axhline(-self.rmsSigma, color=color['rms'], linewidth=3)
        ax1.axhline(+self.iqrSigma, color=color['iqr'], linewidth=3)
        ax1.axhline(-self.iqrSigma, color=color['iqr'], linewidth=3)
        
        ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)
        ax2.hist(self.diffmags, bins=25, range=diffRange,
                 orientation='horizontal', histtype='stepfilled',
                 normed=True, color=color['bright'])
        ax2.set_xlabel("relative # / bin")
        
        yv = np.linspace(diffRange[0], diffRange[1], 100)
        ax2.plot(scipy.stats.norm.pdf(yv, scale=self.rmsSigma), yv,
                 marker='', linestyle='-', linewidth=3, color=color['rms'],
                 label="PA1(RMS) = %4.2f %s" % (self.rmsSigma, 'mmag'))
        ax2.plot(scipy.stats.norm.pdf(yv, scale=self.iqrSigma), yv,
                 marker='', linestyle='-', linewidth=3, color=color['iqr'],
                 label="PA1(IQR) = %4.2f %s" % (self.iqrSigma, 'mmag'))
        ax2.set_ylim(*diffRange)
        ax2.legend()
        #    ax1.set_ylabel(u"12-pixel aperture magnitude diff (mmag)")
        #    ax1.set_xlabel(u"12-pixel aperture magnitude")
        ax1.set_xlabel("psf magnitude")
        ax1.set_ylabel("psf magnitude diff (mmag)")
        for label in ax2.get_yticklabels():
            label.set_visible(False)

        plt.suptitle("PA1: %s" % outputPrefix.rstrip('_'))
        plotPath = "%s%s" % (outputPrefix, "PA1.png")
        plt.savefig(plotPath, format="png")
       # plt.close(fig)




     #  print 'self.diffmags, self.rmsSigma,self.iqrSigma', self.diffmags, self.rmsSigma, self.iqrSigma

    #    for i in range( numRandomShuffles):
    #        groupnb1=0
    #        groupnb2=0
    #        while groupnb1 == groupnb2:
    #            rand1=random.randint(0,len(self.brightsources))
    #            rand2=random.randint(0,len(self.brightsources))
    #            groupnb1=self.brightsources['Nb_group'][rand1]
    #            groupnb2=self.brightsources['Nb_group'][rand2]
    #            if groupnb1==groupnb2:
    #                print ' groupnb1 == groupnb2 !!'
    #                print 'groupnb1, groupnb2', groupnb1, groupnb2
    #        inside_grp1,= np.where(self.brightsources['Nb_group'] == groupnb1)
    #        inside_grp2,= np.where(self.brightsources['Nb_group'] == groupnb2)
    #        mags_grp1 = self.brightsources[sourceFluxField+'_mag'][inside_grp1]
    #        dmag_grp1 = mags_grp1-np.mean(mags_grp1)
    #        mags_grp2 = self.brightsources[sourceFluxField+'_mag'][inside_grp2]
    #        dmag_grp2 = mags_grp1-np.mean(mags_grp2)

            #diffmags.append(self.brightsources[sourceFluxField+'_mag'][rand2]-self.brightsources[sourceFluxField+'_mag'][rand1])

    #        diffmags.append(dmag_grp1-dmag_grp2)


    def monPA1(self,  numRandom=False):
        self.magmean = []
        self.dmagsRMS=[]
        if numRandom:
            idgrps=[random.randint(0,len(set(self.brightsources['Nb_group']))) for n in range(numRandom)]
            idgrps=self.brightsources['Nb_group'][idgrps]
        else:
            idgrps=set(self.brightsources['Nb_group'])

        for groupid in idgrps:
            inside_grp,= np.where(self.brightsources['Nb_group'] == groupid)
          #  snr_grp = self.brightsources[sourceFluxField+'_snr'][inside_grp]
          #  if np.median(snr_grp)>= brightSnr:
            mags_grp = self.brightsources[sourceFluxField+'_mag'][inside_grp]
            dmag_grp = mags_grp - np.mean(mags_grp)
            self.magmean.append(np.mean( mags_grp))
            self.dmagsRMS.append(np.std(dmag_grp))
                #self.dmagsIQR.append(np.subtract.reduce(np.percentile(dmag_grp, [75, 25])) / (scipy.stats.norm.ppf(0.75)*2))
        self.dmagsRMS = np.array(self.dmagsRMS)*1000
        print 'len(self.dmagsRMS)',len(self.dmagsRMS)
        self.PA1=np.median(self.dmagsRMS)

        print 'PA1', self.PA1 # PA1 I have understood in LSST srd...

        level="design"
        print '======================================================='
        print 'Comparison against *'+level+'* requirements (avec calcul magsRMS, pas celui dans validate_drp).'
        print 'Measured           Required      Passes        '      
        print 'PA1 : ', self.PA1, 'mmag < ', srdSpec['PA1'][level], 'mmag ==', self.PA1 < srdSpec['PA1'][level]
       
        magDiffs=list(self.dmagsRMS) #magDiffs
       # print('lenmagdiffs', len(magDiffs), 'magDiffs',magDiffs)
        PA2_spec = srdSpec['PA2']
        PF1_percentiles = 100 - np.asarray([srdSpec['PF1'][l] for l in srdSpec['levels']])
        PF1_percentile = 100 - np.asarray(srdSpec['PF1'][level])
        PA2_measured =  np.percentile(np.abs(magDiffs), PF1_percentile)
        PF1_measured = 100*np.mean(np.asarray(magDiffs) > srdSpec['PA2'][level])
        self.PF1=PF1_measured
        self.PA2=PA2_measured
        print 'PF1 : ', self.PF1, '%    < ', srdSpec['PF1'][level], '%    ==', self.PF1 < srdSpec['PF1'][level]
        print 'PA2 : ', self.PA2, 'mmag < ', srdSpec['PA2'][level], 'mmag ==', self.PA2 < srdSpec['PA2'][level]
      



        plt.figure()
        digits=1000.
        plt.hist(self.dmagsRMS, histtype='stepfilled')#, bins=20)
        plt.axvline(self.PA1, color=color['rms'], linewidth=3, label='PA1(RMS) = '+str(int(self.PA1*digits)/digits)+'mmag')

        plt.axvline(srdSpec['PA1']['design'], color='green', linewidth=3, label="PA1 'design' (LSST srd)= "+str(srdSpec['PA1']['design'])+'mmag')

        plt.xlabel('Mag RMS (mmag)')
        plt.ylabel('#/bin')
        plt.legend()
        plt.show()


    def afficher_phot_req(self, level="design"):
        print '======================================================='
        print 'Comparison against *'+level+'* requirements.'
        print 'Measured           Required      Passes        '      
        print 'PA1 : ', self.PA1, 'mmag < ', srdSpec['PA1'][level], 'mmag ==', self.PA1 < srdSpec['PA1'][level]
       
        magDiffs=list(self.diffmags) #magDiffs
       # print('lenmagdiffs', len(magDiffs), 'magDiffs',magDiffs)
        PA2_spec = srdSpec['PA2']
        PF1_percentiles = 100 - np.asarray([srdSpec['PF1'][l] for l in srdSpec['levels']])
        PF1_percentile = 100 - np.asarray(srdSpec['PF1'][level])
        PA2_measured =  np.percentile(np.abs(magDiffs), PF1_percentile)
        PF1_measured = 100*np.mean(np.asarray(magDiffs) > srdSpec['PA2'][level])
        self.PF1=PF1_measured
        self.PA2=PA2_measured
        print 'PF1 : ', self.PF1, '%    < ', srdSpec['PF1'][level], '%    ==', self.PF1 < srdSpec['PF1'][level]
        print 'PA2 : ', self.PA2, 'mmag < ', srdSpec['PA2'][level], 'mmag ==', self.PA2 < srdSpec['PA2'][level]
        ### to do
        print 'PA3 : ' # spatial uniformity of photometric zeropoints
        print 'PF2 : ' # 
        print 'PA4 : ' # 
        print 'PA5 : ' # band to band (flux ratio) photometric calibration  (besoin de plusieurs filtres a la fois)
        print 'PA6 : ' # overall external absolute photometry


class Validation_plots(LoadDataValidation):
    def __init__(self,sources):
        self.sources = sources

    def plotAstrometry(self,dist, mag, snr, fit_params=None, brightSnr=100,
                       outputPrefix=""):
        """Plot angular distance between matched sources from different exposures.
        
        Creates a file containing the plot with a filename beginning with `outputPrefix`.
        
        Parameters
        ----------
        dist : list or numpy.array
        Separation from reference [mas]
        mag : list or numpy.array
        Mean magnitude of PSF flux
        snr : list or numpy.array
        Median SNR of PSF flux
        fit_params : list or numpy.array
        Fit parameters to display
        brightSnr : float, optional
        Minimum SNR for a star to be considered "bright".
        outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
        'Cfht_output_r_check_astrometry.png'
        """
        bright, = np.where(np.asarray(snr) > brightSnr)
     
        numMatched = len(dist)
        dist_median = np.median(dist)
        bright_dist_median = np.median(np.asarray(dist)[bright])

        fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(18, 12))

        ax[0].hist(dist, bins=100, color=color['all'],
                   histtype='stepfilled', orientation='horizontal')

        ax[0].hist(np.asarray(dist)[bright], bins=100, color=color['bright'],
                   histtype='stepfilled', orientation='horizontal',
                   label='SNR > %.0f' % brightSnr)

        ax[0].set_ylim([0., 500.])
        ax[0].set_ylabel("Distance [mas]")
        ax[0].set_title("Median : %.1f, %.1f mas" %
                        (bright_dist_median, dist_median),
                        x=0.55, y=0.88)
        plotOutlinedLinesHorizontal(ax[0], dist_median, bright_dist_median)

        ax[1].scatter(snr, dist, s=10, color=color['all'], label='All')
        ax[1].scatter(np.asarray(snr)[bright], np.asarray(dist)[bright], s=10,
                      color=color['bright'],
                      label='SNR > %.0f' % brightSnr)
        ax[1].set_xlabel("SNR")
        ax[1].set_xscale("log")
        ax[1].set_ylim([0., 500.])
        ax[1].set_title("# of matches : %d, %d" % (len(bright), numMatched))

        w, = np.where(dist < 200)
        plotAstromErrModelFit(snr[w], dist[w], fit_params=fit_params, ax=ax[1])
  

        ax[1].legend(loc='upper right')

        ax[1].axvline(brightSnr, color='red', linewidth=4, linestyle='dashed')
        plotOutlinedLinesHorizontal(ax[1], dist_median, bright_dist_median)
       
        plt.suptitle("Astrometry Check : %s" % outputPrefix.rstrip('_'), fontsize=30)
      
        plotPath = outputPrefix+"check_astrometry.png"
        plt.savefig(plotPath, format="png")
        #plt.close(fig)
       


    def plotPhotometry(self,mag, snr, mmagerr, mmagrms, brightSnr=100,
                       fit_params=None,
                       filterName='Magnitude',
                       outputPrefix=""):
        """Plot photometric RMS for matched sources.
        
        Parameters
        ----------
        snr : list or numpy.array
        Median SNR of PSF flux
        mag : list or numpy.array
        Average Magnitude
        mmagerr : list or numpy.array
        Average Magnitude uncertainty [millimag]
        mmagrms ; list or numpy.array
        Magnitude RMS across visits [millimag]
        fit_params : list or numpy.array
        Fit parameters for photometry error model
        brightSnr : float, optional
        Minimum SNR for a star to be considered "bright".
        filterName : str, optional
        Name of the observed filter to use on axis labels.
        outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
        'Cfht_output_r_check_photometry.png'
        """

        bright, = np.where(np.asarray(snr) > brightSnr)

        numMatched = len(mag)
        mmagrms_median = np.median(mmagrms)
        bright_mmagrms_median = np.median(np.asarray(mmagrms)[bright])
        
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(18, 16))
        ax[0][0].hist(mmagrms, bins=100, range=(0, 500), color=color['all'], label='All',
                      histtype='stepfilled', orientation='horizontal')
        ax[0][0].hist(np.asarray(mmagrms)[bright], bins=100, range=(0, 500), color=color['bright'],
                      label='SNR > %.0f' % brightSnr,
                      histtype='stepfilled', orientation='horizontal')
        ax[0][0].set_ylim([0, 500])
        ax[0][0].set_ylabel("RMS [mmag]")
        ax[0][0].set_title("Median : %.1f, %.1f mmag" %
                           (bright_mmagrms_median, mmagrms_median),
                           x=0.55, y=0.88)
        plotOutlinedLinesHorizontal(ax[0][0], mmagrms_median, bright_mmagrms_median)
        
        ax[0][1].scatter(mag, mmagrms, s=10, color=color['all'], label='All')
        ax[0][1].scatter(np.asarray(mag)[bright], np.asarray(mmagrms)[bright],
                         s=10, color=color['bright'],
                         label='SNR > %.0f' % brightSnr)
        
        ax[0][1].set_xlabel("%s [mag]" % filterName)
        ax[0][1].set_ylabel("RMS [mmag]")
        ax[0][1].set_xlim([17, 24])
        ax[0][1].set_ylim([0, 500])
        ax[0][1].set_title("# of matches : %d, %d" % (len(bright), numMatched))
        ax[0][1].legend(loc='upper left')
        plotOutlinedLinesHorizontal(ax[0][1], mmagrms_median, bright_mmagrms_median)
        
        ax[1][0].scatter(mmagrms, mmagerr, s=10, color=color['all'], label=None)
        ax[1][0].scatter(np.asarray(mmagrms)[bright], np.asarray(mmagerr)[bright],
                         s=10, color=color['bright'],
                         label=None)
        ax[1][0].set_xscale('log')
        ax[1][0].set_yscale('log')
        ax[1][0].plot([1., 999.], [1., 999.],
                      linestyle='--', color='black', linewidth=2) #
        ax[1][0].set_xlabel("RMS [mmag]")
        ax[1][0].set_ylabel("Median Reported Magnitude Err [mmag]")
        brightSnrInMmag = 2.5*np.log10(1 + (1./brightSnr)) * 1000
        #print ' brightSnrInMmag, brightSnr', brightSnrInMmag, brightSnr
        ax[1][0].axhline(brightSnrInMmag, color='red', linewidth=4, linestyle='dashed',
                         label=r'SNR > %.0f = $\sigma_{\rm mmag} < $ %0.1f' % (brightSnr, brightSnrInMmag))
        ax[1][0].set_xlim([1, 500])
        ax[1][0].set_ylim([1, 500])
        ax[1][0].legend(loc='upper center')
        
        ax[1][1].scatter(mag, mmagerr, color=color['all'], label=None)
        ax[1][1].set_yscale('log')
        ax[1][1].scatter(np.asarray(mag)[bright], np.asarray(mmagerr)[bright],
                         s=10, color=color['bright'],
                         label=None,
                         )
        ax[1][1].set_xlabel("%s [mag]" % filterName)
        ax[1][1].set_ylabel("Median Reported Magnitude Err [mmag]")
        ax[1][1].set_xlim([17, 24])
        ax[1][1].set_ylim([1, 500])
        ax[1][1].axhline(brightSnrInMmag, color='red', linewidth=4, linestyle='dashed',
                         label=r'$\sigma_{\rm mmag} < $ %0.1f' % (brightSnrInMmag))
        
        ax2 = ax[1][1].twinx()
        ax2.scatter(mag, snr,
                    color=color['all'], facecolor='none',
                    marker='.', label=None)
        ax2.scatter(np.asarray(mag)[bright], np.asarray(snr)[bright],
                    color=color['bright'], facecolor='none',
                    marker='.', label=None)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("SNR")
        ax2.axhline(brightSnr, color='red', linewidth=2, linestyle='dashed',
                    label=r'SNR > %.0f' % (brightSnr))
        
        w, = np.where(mmagerr < 200)
        plotPhotErrModelFit(mag[w], mmagerr[w], fit_params=fit_params, ax=ax[1][1])
        ax[1][1].legend(loc='upper left')
        
        plt.suptitle("Photometry Check : %s" % outputPrefix.rstrip('_'), fontsize=30)
        plotPath = outputPrefix+"check_photometry.png"
        plt.savefig(plotPath, format="png")
      #  plt.show()
       # plt.close(fig)


    def plotVisitVsTime(self,
                        outputPrefix="OutputPlots/"):
        mjd = self.sources['MJD-OBS']
        visit = self.sources['visit']

        plt.figure(figsize=(10,8))
        time = mjd-min(mjd)
        plt.scatter(time, visit )
        plt.xlabel('t-tmin (mjd)')
        plt.ylabel('visit')
        plt.title('Visit in function of Time (mjd)')
        plotPath = outputPrefix + 'VisitVsTime.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure(figsize=(12,10))
        plt.scatter(mjd, visit )
        plt.xlabel('t (mjd)')
        plt.ylabel('visit')
        plt.title('Visit in function of Time (mjd)')
        plotPath = outputPrefix + 'VisitVsTime_mjd.png'
        plt.savefig(plotPath, format="png")
     #plt.show()



    def plotAstromPhotRMSvsTimeCcd(self,
                                   brightSnr=100, srcFluxField='base_PsfFlux',
                                   outputPrefix="",
                                   zoom=False, dico=False):
            
        sourceFluxField=srcFluxField
        sizelegend=12 # taille des legendes
        digits=1000. # precision des valeurs dans les legendes des histos
        
        if zoom: 
            outputPrefix=outputPrefix+"zoom_" # pour sauvegarder aussi les plots avec des zoom
        compt = 0
        goodmjd = self.sources['MJD-OBS']
        
        deltaRAcosdecs =  self.deltaRAcosdecs *radToMas
        deltaDecs =  self.deltaDecs *radToMas
        
        ccds = self.sources['ccd']
        sourcemag =  self.sources[sourceFluxField+'_mag']*1000
        sourcedmag = self.deltamags*1000
        sourcesnr =  self.sources[sourceFluxField+'_snr']

        FluxMag0s = self.sources['FLUXMAG0']
        FluxMag0Errs = self.sources['FLUXMAG0ERR']
        
        grpMeanRAcosdec =  self.RAcosDec_mean
        grpMeanDec = self.Dec_mean
        groupRMSracosdec = self.RAcosDec_RMS *radToMas
        groupRMSdec = self.Dec_RMS  *radToMas
        posRMS = self.posRMS *radToMas
        
        medsnr = self.snr_med
        medsnrlong =  self.medsnrlong
        Psf_fwhm = self.sources['PSF-FWHM']
        visits = self.sources['visit']
        grpMeanPsf_fwhm = self.Psf_fwhm_mean
  
        

        bright, = np.where(np.asarray(medsnr) > brightSnr)
        
        groupRMSracosdec_bright = np.array(groupRMSracosdec)[bright]
        groupRMSdec_bright = np.array(groupRMSdec)[bright]
        posRMS_bright = np.array(posRMS)[bright]
        
        nb_sigma=5
        bright_outliers, = np.where((np.asarray(posRMS_bright-np.median(posRMS_bright)) < nb_sigma*np.std(posRMS_bright))) and  np.where((np.asarray(groupRMSracosdec_bright-np.median(groupRMSracosdec_bright)) < nb_sigma*np.std(groupRMSracosdec_bright))) and np.where((np.asarray(groupRMSdec_bright-np.median(groupRMSdec_bright)) < nb_sigma*np.std(groupRMSdec_bright)))
        
        posRMS_bright_outliers = posRMS_bright[bright_outliers]
        groupRMSracosdec_bright_outliers = groupRMSracosdec_bright[bright_outliers]
        groupRMSdec_bright_outliers = groupRMSdec_bright[bright_outliers]
        
        grpMeanRAcosdec=np.array(grpMeanRAcosdec)
        
        deltaRAcosdecs = np.array(deltaRAcosdecs)
        deltaDecs = np.array(deltaDecs)
        Psf_fwhm = np.array(Psf_fwhm)
        
    #brightallsnr, = np.where(np.asarray(sourcesnr) > brightSnr)
        brightallsnr, = np.where(np.asarray(medsnrlong) > brightSnr) #pour avoir des plots comparables (meme cut sur les snr medianes)
    #    print '  brightallsnr',  brightallsnr
        #plt.close('all')
        plt.figure(figsize=(12,12))
        if zoom:
            bri,=np.where(np.array(medsnr) >= brightSnr)
            plt.hist(np.array(medsnr)[bri],bins=50, histtype ='stepfilled', alpha=0.8, color='b',label='Median='+str(int(np.median(np.array(medsnr)[bri])*digits)/digits)+' \nRMS='+str(int(np.std(np.array(medsnr)[bri])*digits)/digits)+' \nmin='+str(min(np.array(medsnr)[bri]))+' \nmax='+str(max(np.array(medsnr)[bri])))
            plt.title('SNR median distribution (focus on bright median SNR)')
        else :
            plt.hist(medsnr,bins=50, histtype ='stepfilled', alpha=0.8, color='b',label='Median='+str(int(np.median(medsnr)*digits)/digits)+' \nRMS='+str(int(np.std(medsnr)*digits)/digits)+' \nmin='+str(min(np.array(medsnr)))+' \nmax='+str(max(np.array(medsnr))))
            plt.title('SNR median distribution')
        plt.xlabel('median snr')
        plt.ylabel('# / bin')
        plt.legend(prop={'size':sizelegend})
        plotPath = outputPrefix+'SNRgroupMedian.png'
        plt.savefig(plotPath, format="png")
        
        
        plt.figure(figsize=(12,12))
        if zoom:
            bri,= np.where(np.array(sourcesnr) >= brightSnr)
            plt.hist(np.array(sourcesnr)[bri], bins=50, histtype ='stepfilled', alpha=0.8, color='b',label='Median='+str(int(np.median(np.array(sourcesnr)[bri])*digits)/digits)+' \nRMS='+str(int(np.std(np.array(sourcesnr)[bri])*digits)/digits)+' \nmin='+str(min(np.array(sourcesnr)[bri]))+' \nmax='+str(max(np.array(sourcesnr)[bri])))
            plt.title('SNR sources distribution (focus on bright median SNR)')
        else :
            plt.hist(sourcesnr, bins=50, histtype ='stepfilled', alpha=0.8, color='b',label='Median='+str(int(np.median(sourcesnr)*digits)/digits)+' \nRMS='+str(int(np.std(sourcesnr)*digits)/digits)+' \nmin='+str(min(np.array(sourcesnr)))+' \nmax='+str(max(np.array(sourcesnr))))
            plt.title('SNR sources distribution')
        plt.xlabel('source snr')
        plt.ylabel('# / bin')
        plt.legend(prop={'size':sizelegend})
        plotPath = outputPrefix+'SNRsources.png'
        plt.savefig(plotPath, format="png")
        # plt.show()



        plt.figure(figsize=(12,12))
        plt.title('PSF FWHM')
        plt.hist(Psf_fwhm, histtype ='stepfilled', alpha=0.8, color='b')
        histB = plt.hist(Psf_fwhm[brightallsnr] ,histtype ='stepfilled',alpha=0.8,color='r')
        plt.axvline(np.median(Psf_fwhm), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(Psf_fwhm)*digits)/digits)+'as')
        plt.axvline(np.median(Psf_fwhm[brightallsnr]), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(Psf_fwhm[brightallsnr] )*digits)/digits)+'as')
        plt.xlabel('psf fwhm (as)')
        plt.ylabel('# / bin')
        plt.legend(prop={'size':sizelegend})
        if zoom:
            plt.ylim(0., max( histB[0]))
        plotPath = outputPrefix+'PsfFwhmvshist.png'
        plt.savefig(plotPath, format="png")
        
        
        plt.figure(figsize=(12,12))
        plt.title('PSF FWHM vs deltaRAcosdecs')
        plt.scatter(deltaRAcosdecs, Psf_fwhm, color=color['all'], label='all')
        plt.scatter(deltaRAcosdecs[brightallsnr], Psf_fwhm[brightallsnr], color=color['bright'], label='bright')
        plt.legend()
        plt.xlabel('DeltaRaCosDec (mas)')
        plt.ylabel('psf fwhm (as)')
        if zoom:
            plt.xlim(-75.,75.)
            plt.ylim(min(Psf_fwhm),max(Psf_fwhm))
        plotPath = outputPrefix+'PsfFwhmvsdeltaRacosDec.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure(figsize=(12,12))
        plt.title('PSF FWHM vs deltaDecs')
        plt.scatter(deltaDecs, Psf_fwhm, color=color['all'], label='all')
        plt.scatter(deltaDecs[brightallsnr], Psf_fwhm[brightallsnr], color=color['bright'], label='bright')
        plt.legend()
        plt.xlabel('DeltaDec (mas)')
        plt.ylabel('psf fwhm (as)')
        if zoom:
            plt.xlim(-75.,75.)
            plt.ylim(min(Psf_fwhm),max(Psf_fwhm))
        plotPath = outputPrefix+'PsfFwhmvsdeltaDec.png'
        plt.savefig(plotPath, format="png")


       # plt.close('all')
        plt.figure()
        plt.title('racosdec')
        plt.scatter(grpMeanRAcosdec, groupRMSracosdec, color=color['all'], label='all')
        plt.scatter(grpMeanRAcosdec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
        plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')
        plt.legend()
        plt.xlabel('mean ra cosdec')
        plt.ylabel('rms ra cos dec')
        plotPath = outputPrefix+'RAcosDecRMSvsRacosDecMean.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure()
        plt.title('racosdec')
        plt.scatter(grpMeanRAcosdec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
        plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')
        
        plt.legend()
        plt.xlabel('mean ra cosdec')
        plt.ylabel('rms ra cos dec')
        
        plotPath = outputPrefix+'RAcosDecRMSvsRacosDecMean_bright.png'
        plt.savefig(plotPath, format="png")


# grpMeanShapex

   # plt.figure()
   # plt.title('test corr shape_x /RA')
   # plt.scatter(groupRMSracosdec,grpMeanShapex, color=color['all'], label='all')
    # #marche pas plt.scatter(groupRMSracosdec[bright], grpMeanShapex[bright], color=color['bright'], label='bright')#
    #plt.legend()
   # plt.xlabel('mean ra cosdec')
   # plt.ylabel('Shape_x')
        """
        plt.figure()
        plt.title('racosdec')
        plt.scatter(grpMeanDec, groupRMSracosdec, color=color['all'], label='all')
        plt.scatter(grpMeanDec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
        plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')
        plt.legend()
        
        plt.legend()
        plt.xlabel('dec')
        plt.ylabel('rms racos')
        """
  #  plt.show()
    
    #plot(x,y,"k.")
  #  y_av = movingaverage(groupRMSracosdec, 100)
   # plt.plot(grpMeanRAcosdec, y_av,"r")
  #  plt.xlim(0,1000)
   # plt.grid(True)
   # plt.show()


        plt.figure(figsize=(12,12))
        plt.title('RMS racosdec vs RMS dec')
        plt.scatter(groupRMSracosdec, groupRMSdec, alpha=0.8, color='b', label='all')
        plt.scatter(groupRMSracosdec_bright, groupRMSdec_bright, alpha=0.8, color='r', label='bright')
        plt.scatter(groupRMSracosdec_bright_outliers, groupRMSdec_bright_outliers, alpha=0.8, color='g', label='bright + 5 sigma outliers')
        plt.xlabel('RMS RAcosDec (mas)')
        plt.ylabel('RMS Dec (mas)')
        plt.xlim(0,max(groupRMSracosdec))
        plt.ylim(0,max(groupRMSdec))
        if zoom:
            plt.xlim(0.,55.)
            plt.ylim(0.,55.)
        plt.legend(prop={'size':sizelegend})
        plotPath = outputPrefix+'RAcosDecRMSvsDecRMS.png'
        plt.savefig(plotPath, format="png")
    #plt.show()


        plt.figure()
        plt.title('racosdec')
        plt.hist(groupRMSracosdec, bins=50, histtype ='stepfilled', alpha=0.8, color='b')# label='RMS='+str(int(np.std(groupRMSracosdec)*digits)/digits)+'mAcs\nMean='+str(int(np.mean(groupRMSracosdec)*digits)/digits)+'mAcs',alpha=0.5)#,histtype ='stepfilled',alpha=0.8,color='r')
        histB = plt.hist(groupRMSracosdec_bright,histtype ='stepfilled',alpha=0.8,color='r')
        plt.axvline(np.median(groupRMSracosdec), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(groupRMSracosdec)*digits)/digits)+'mas')
        plt.axvline(np.median(groupRMSracosdec_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(groupRMSracosdec_bright)*digits)/digits)+'mas')
        plt.hist(groupRMSracosdec_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
        plt.axvline(np.median( groupRMSracosdec_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( groupRMSracosdec_bright_outliers )*digits)/digits)+'mas')
        plt.xlabel('RMS RAcosDec')
        plt.ylabel('#/bin')
        plt.legend(prop={'size':sizelegend})
        if zoom:
            plt.xlim(0.,55.)
            plt.ylim(0., max( histB[0]))
        plotPath = outputPrefix+'RAcosDecRMS.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure()
        plt.title('dec')
        plt.hist(groupRMSdec, bins=50, histtype ='stepfilled', alpha=0.8, color='b')
        histB = plt.hist(groupRMSdec_bright,histtype ='stepfilled', alpha=0.8,color='r')
        plt.hist(groupRMSdec_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
        plt.axvline(np.median(groupRMSdec), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(groupRMSdec)*digits)/digits)+'mas')
        plt.axvline(np.median(groupRMSdec_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(groupRMSdec_bright)*digits)/digits)+'mas')
        plt.axvline(np.median( groupRMSdec_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( groupRMSdec_bright_outliers )*digits)/digits)+'mas')
        plt.xlabel('RMS Dec')
        plt.ylabel('#/bin')
        plt.legend(prop={'size':sizelegend})
        if zoom:
            plt.xlim(0.,55.)
            plt.ylim(0., max( histB[0]))
        plotPath = outputPrefix+'DecRMS.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure()
        plt.title('RMS distances (equiv plotAstrometry)')
        plt.hist(posRMS, bins=50, histtype ='stepfilled', alpha=0.8, color='b')
        histB = plt.hist(posRMS_bright,histtype ='stepfilled',alpha=0.8,color='r')
        plt.hist(posRMS_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
        plt.axvline(np.median(posRMS), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(posRMS)*digits)/digits)+'mas')
        plt.axvline(np.median(posRMS_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(posRMS_bright)*digits)/digits)+'mas')
        plt.axvline(np.median( posRMS_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( posRMS_bright_outliers )*digits)/digits)+'mas')

        plt.xlabel('RMS distance')
        plt.ylabel('#/bin')
        plt.legend(prop={'size':sizelegend})
        if zoom:
            plt.xlim(0.,55.)
            plt.ylim(0., max( histB[0]))
        plotPath = outputPrefix+'DistanceRMS.png'
        plt.savefig(plotPath, format="png")

        if zoom:
            return

        # quelques prints
        print('Nombre total de sources :',  len(Psf_fwhm))
        print('Nombre total de sources ayant passe le cut median de SNR (bright) :',  len(np.array(Psf_fwhm)[brightallsnr]))
        print('Nombre total de groupes :',  len(grpMeanPsf_fwhm))
        print('Nombre total de groupes ayant passe le cut median de SNR (bright) :',  len(np.array(grpMeanPsf_fwhm)[bright]))
        #  print('ccd',ccds)
        #  print("LONGUEUR DELTARAS", len(deltaRAcosdecs))
        #  print("LONGUEUR SOURCEMAG", len(sourcemag))
        #  print("LONGUEUR sourceSNR", len(sourcesnr))
        
        # ### astrometry and photometry vs time ####
        times=[]
        dic_time={}
        for i in range(len(deltaRAcosdecs)):
            dic_time[goodmjd[i]]={}
            dic_time[goodmjd[i]]['dra']=[]
            dic_time[goodmjd[i]]['ddec']=[]
            dic_time[goodmjd[i]]['sourcedmag']=[]
            dic_time[goodmjd[i]]['sourcesnr']=[]
            dic_time[goodmjd[i]]['FluxMag0s']=[]
            dic_time[goodmjd[i]]['FluxMag0Errs']=[]

        for i in range(len(deltaRAcosdecs)):
            dic_time[goodmjd[i]]['dra'].append(deltaRAcosdecs[i])
            dic_time[goodmjd[i]]['ddec'].append(deltaDecs[i])
            dic_time[goodmjd[i]]['sourcedmag'].append(sourcedmag[i])
            dic_time[goodmjd[i]]['sourcesnr'].append(sourcesnr[i])
            dic_time[goodmjd[i]]['FluxMag0s'].append(FluxMag0s[i])
            dic_time[goodmjd[i]]['FluxMag0Errs'].append(FluxMag0Errs[i])

        for time in dic_time.keys():
            times.append(time)
            dic_time[time]['dra_med'] = np.median(dic_time[time]['dra'])
            dic_time[time]['dra_rms'] = np.std(dic_time[time]['dra'])
            dic_time[time]['ddec_med'] = np.median(dic_time[time]['ddec'])
            dic_time[time]['ddec_rms'] = np.std(dic_time[time]['ddec'])
            dic_time[time]['sourcedmag_med'] = np.median(dic_time[time]['sourcedmag'])
            dic_time[time]['sourcedmag_rms'] = np.std(dic_time[time]['sourcedmag'])
            dic_time[time]['sourcesnr_med'] = np.median(dic_time[time]['sourcesnr'])
            dic_time[time]['sourcesnr_rms'] = np.std(dic_time[time]['sourcesnr']) # barre err inutile
            dic_time[time]['FluxMag0s_mean'] = np.mean(dic_time[time]['FluxMag0s'])
            dic_time[time]['FluxMag0s_med'] = np.median(dic_time[time]['FluxMag0s'])
            dic_time[time]['FluxMag0s_rms'] = np.std(dic_time[time]['FluxMag0s'])
            dic_time[time]['FluxMag0Errs_mean'] = np.mean(dic_time[time]['FluxMag0Errs'])
            dic_time[time]['FluxMag0Errs_med'] = np.median(dic_time[time]['FluxMag0Errs'])
            dic_time[time]['FluxMag0Errs_rms'] = np.std(dic_time[time]['FluxMag0Errs'])
            
            #  plt.errorbar(time,  dic_time[time]['dra_med'], xerr=None, yerr=dic_time[time]['dra_rms'],linestyle='',alpha=0.75, marker='o', color='b')
            #    plt.errorbar(time,  dic_time[time]['ddec_med'], xerr=None, yerr=dic_time[time]['ddec_rms'],linestyle='',alpha=0.75, marker='o', color='r') 
    #plt.figure(figsize=(10,8))
 
        sizelegend=12
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(16,9))
        # f.subplots_adjust(hspace=0) # pour enlever l'espace entre les graphiques (vertical)
        # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        for time in dic_time.keys():
            ax1.scatter(time, dic_time[time]['dra_rms'],alpha=0.85, marker='o', color='b')
            ax1.scatter(time, dic_time[time]['ddec_rms'],alpha=0.85, marker='o', color='r')
        ax1.scatter([], [] , color='b',marker='o',label='$\delta(RA) \cos(dec)$')
        ax1.scatter([], [] , color='r',marker='o',label='$\delta(Dec)$')
        ax1.legend(prop={'size':sizelegend})
        ax1.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
        ax1.set_ylabel('RMS (mas)')
    # ax1.set_xlabel('time (mjd)')
        ax1.set_title('RMS Dra/Ddec in function of time')
        for time in dic_time.keys():
            ax2.scatter(time, dic_time[time]['sourcedmag_rms'],alpha=0.85, marker='o', color='g') #
        ax2.scatter([], [], marker='o', color='g',label='source $\Delta$Magnitude RMS') #
        ax2.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
        ax2.set_ylabel('RMS $\Delta$Mag (mmag)')
    #plt.xlabel('time (mjd)')
        ax2.set_title('RMS mag in function of time')
        ax2.legend(prop={'size':sizelegend})
        for time in dic_time.keys():
            ax3.scatter(time, dic_time[time]['sourcesnr_med'],alpha=0.85, marker='o', color='k') #
        # ax3.errorbar(time,  dic_time[time]['sourcesnr_med'], xerr=None, yerr=dic_time[time]['sourcesnr_rms'],linestyle='',alpha=0.75, marker='o', color='b')
        ax3.scatter([], [], marker='o', color='b',label='source SNR Median') #
        ax3.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
        ax3.set_ylabel('SNR')
        ax3.set_xlabel('time (mjd)')
        ax3.set_title('Median SNR in function of time')
        ax3.legend(prop={'size':sizelegend})
        plotPath = outputPrefix+'AstrometryVsTime.png'
        plt.savefig(plotPath, format="png")
        
        # ### astrometry and photometry vs CCD ####
        dic_ccd={}
        for i in range(len(deltaRAcosdecs)):
            dic_ccd[ccds[i]]={}
            dic_ccd[ccds[i]]['dra']=[]
            dic_ccd[ccds[i]]['ddec']=[]
            dic_ccd[ccds[i]]['sourcedmag']=[]
            dic_ccd[ccds[i]]['sourcesnr']=[]
            dic_ccd[ccds[i]]['FluxMag0s']=[]
            dic_ccd[ccds[i]]['FluxMag0Errs']=[]

        for i in range(len(deltaRAcosdecs)):
            dic_ccd[ccds[i]]['dra'].append(deltaRAcosdecs[i])
            dic_ccd[ccds[i]]['ddec'].append(deltaDecs[i])
            dic_ccd[ccds[i]]['sourcedmag'].append(sourcedmag[i])
            dic_ccd[ccds[i]]['sourcesnr'].append(sourcesnr[i])
            dic_ccd[ccds[i]]['FluxMag0s'].append(FluxMag0s[i])
            dic_ccd[ccds[i]]['FluxMag0Errs'].append(FluxMag0Errs[i])
            
        for ccd in dic_ccd.keys():
            #ccds.append(ccd)
            dic_ccd[ccd]['dra_med'] = np.median(dic_ccd[ccd]['dra'])
            dic_ccd[ccd]['dra_rms'] = np.std(dic_ccd[ccd]['dra'])
            dic_ccd[ccd]['ddec_med'] = np.median(dic_ccd[ccd]['ddec'])
            dic_ccd[ccd]['ddec_rms'] = np.std(dic_ccd[ccd]['ddec'])
            dic_ccd[ccd]['sourcedmag_med'] = np.median(dic_ccd[ccd]['sourcedmag'])
            dic_ccd[ccd]['sourcedmag_rms'] = np.std(dic_ccd[ccd]['sourcedmag'])
            dic_ccd[ccd]['sourcesnr_med'] = np.median(dic_ccd[ccd]['sourcesnr'])
            dic_ccd[ccd]['sourcesnr_rms'] = np.std(dic_ccd[ccd]['sourcesnr'])# barre err inutile
            dic_ccd[ccd]['FluxMag0s_mean'] = np.mean(dic_ccd[ccd]['FluxMag0s'])
            dic_ccd[ccd]['FluxMag0s_med'] = np.median(dic_ccd[ccd]['FluxMag0s'])
            dic_ccd[ccd]['FluxMag0s_rms'] = np.std(dic_ccd[ccd]['FluxMag0s'])
            dic_ccd[ccd]['FluxMag0Errs_mean'] = np.mean(dic_ccd[ccd]['FluxMag0Errs'])
            dic_ccd[ccd]['FluxMag0Errs_med'] = np.median(dic_ccd[ccd]['FluxMag0Errs'])
            dic_ccd[ccd]['FluxMag0Errs_rms'] = np.std(dic_ccd[ccd]['FluxMag0Errs'])
 
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(16,9))
        for ccd in dic_ccd.keys():
            ax1.scatter(ccd, dic_ccd[ccd]['dra_rms'],alpha=0.85, marker='o', color='b')
            ax1.scatter(ccd, dic_ccd[ccd]['ddec_rms'],alpha=0.85, marker='o', color='r')
        ax1.scatter([], [] , color='b',marker='o',label='$\delta(RA) \cos(dec)$')
        ax1.scatter([], [] , color='r',marker='o',label='$\delta(Dec)$')
        ax1.legend(prop={'size':sizelegend})
        ax1.set_xlim(min(ccds)-(max(ccds)-min(ccds))/30., max(ccds)+(max(ccds)-min(ccds))/30.)
        ax1.set_ylabel('RMS (mas)')
        ax1.set_title('RMS Dra/Ddec in function of CCD')
        for ccd in dic_ccd.keys():
            ax2.scatter(ccd, dic_ccd[ccd]['sourcedmag_rms'],alpha=0.85, marker='o', color='g') #
        ax2.scatter([], [], marker='o', color='g',label='source $\Delta$Magnitude RMS') #
        ax2.set_ylabel('RMS $\Delta$Mag (mmag)')
        ax2.set_title('RMS mag in function of CCD')
        ax2.legend(prop={'size':sizelegend})
        for ccd in dic_ccd.keys():
            ax3.scatter(ccd, dic_ccd[ccd]['sourcesnr_med'],alpha=0.85, marker='o', color='k') #
         # ax3.errorbar(ccd,  dic_ccd[ccd]['sourcesnr_med'], xerr=None, yerr=dic_ccd[ccd]['sourcesnr_rms'],linestyle='',alpha=0.75, marker='o', color='b')
        ax3.scatter([], [], marker='o', color='b',label='SNR Median') #
        ax3.set_ylabel('SNR')
        ax3.set_xlabel('CCD #')
        ax3.set_title('Median SNR in function of CCD')
        ax3.legend(prop={'size':sizelegend})
        plotPath = outputPrefix+'AstrometryVsCcd.png'
        plt.savefig(plotPath, format="png")
        
        s=100 # taille points
        plt.figure( figsize=(14,5))
        plt.title('FluxMag0s vs ccd')
        plt.xlabel('ccd')
        plt.ylabel('FluxMag0s')
        for ccd in dic_ccd.keys():
            plt.scatter(ccd, dic_ccd[ccd]['FluxMag0s_med'],s=s,alpha=0.95, marker='+', color='k') 
            plt.scatter(ccd, dic_ccd[ccd]['FluxMag0s_mean'],s=s,alpha=0.95, marker='+', color='r') 
        #plt.scatter(ccd, dic_ccd[ccd]['FluxMag0Errs_med'],alpha=0.85, marker='o', color='k') 
        plt.scatter([], [], marker='+',s=s, color='k',label='Median')
        plt.scatter([], [], marker='+',s=s, color='r',label='Mean')
        plt.legend()
        plotPath = outputPrefix+'ZpVsCCD.png'
        plt.savefig(plotPath, format="png")
        
        plt.figure( figsize=(14,5))
        plt.title('FluxMag0s vs time')
        plt.xlabel('time')
        plt.ylabel('FluxMag0s')
        for time in dic_time.keys():
            plt.scatter(time, dic_time[time]['FluxMag0s_med'],s=s,alpha=0.95, marker='+', color='k') 
            plt.scatter(time, dic_time[time]['FluxMag0s_mean'],s=s,alpha=0.95, marker='+', color='r')
        plt.scatter([], [], marker='+', color='k',s=s,label='Median')
        plt.scatter([], [], marker='+', color='r',s=s,label='Mean')
        plt.legend()
        plotPath = outputPrefix+'ZpVsTime.png'
        plt.savefig(plotPath, format="png")
 




        

if __name__=="__main__":
  
  #  AAA=Comparaison_multiple(liste=[4,5,6,10,20,30,50,60])
  #  AAA.plotPos()

    #LDV = LoadDataValidation(file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/sources_ndarray_grp_16visits.pkl')

    LDV = LoadDataValidation(file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/Tests/sources_ndarray_grp_6visits.pkl')

    

    PLTV = Validation_plots(LDV.sources)
    # PLTV.plotVisitVsTime() # plot fonctionnel
    PLTV.createVariablesForPlots(additional=True, Extended=True)
    #PLTV.createVariablesForPlots(additional=False)
    print 'Data Loaded'
    PLTV.plotAstrometry( PLTV.posRMS*radToMas, PLTV.mag_mean, PLTV.snr_med) # plot fonctionnel, probleme d'initialisation des parametre du fit a regler... voir ce qui est fait dans validate_drp

    PLTV.plotPhotometry(PLTV.mag_mean, PLTV.snr_med, PLTV.magerr_med*1000, PLTV.mag_rms*1000) #idem 
                   # fit_params={'sigmaSysUnits': 'mmag', 'sigmaSys': 0.00067729702863904555, 'm5': 24.578444282189267, 'gammaUnits': '', 'm5Units': 'mag', 'gamma': 0.03914245347332538})
  

    PLTV.plotAstromPhotRMSvsTimeCcd()
    print 'programme termine'

    plt.close('all')
    ZZ=Validation(PLTV)
    ZZ.calcPA1()
