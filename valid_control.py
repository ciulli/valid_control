#import pylab as plt
import numpy as np
import matplotlib.pylab as plt
import scipy.stats

from scipy.optimize import curve_fit
#from ..kdtree import KDTree
import pickle



sourceFluxField='base_PsfFlux'

color = {'all': 'grey', 'bright': 'blue',
         'iqr': 'green', 'rms': 'red'}

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

    def __init__(self, liste=[4,5,6,20]):#liste=[4,5,6,10,20,30,50,60]):
        self.liste=liste
        for i in liste:
            num=str(i)
            print num
            exec("sourcesArray"+num+ " = load_pkl('/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation_phVisits_peudevisites/sources_ndarray_grp_"+num+"visits.pkl')")
       
            exec("self.RA"+num+" = sourcesArray"+num+ "['coord_ra']")
            exec("self.Dec"+num+" = sourcesArray"+num+ "['coord_dec']")


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


class LoadDataValidation:
    def __init__(self, file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/sources_ndarray_grp_16visits.pkl'):
        self.sources = load_pkl(file_path)
        
    def createVariablesForPlots(self):
     #   self.dist=[] #posRMS
      #  self.mag=[] #sourceFluxField
      #  self.snr=[] # group med snr
       # for i, j in enumerate(self.sources['Nb_group']):
          #  if self.sources['Nb_group'][i] == j:
          #      RA_grp.append(self.sources['coord_ra'][i]

        #for i in range(len(self.sources['Nb_group'])): 

        RA_grp = []
        Dec_grp = []
        mag_grp = []
        snr_grp = []

        self.RA_mean= []
        self.Dec_mean = []
        self.posRMS = []
        self.mag_mean = []
        self.snr_med = []

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
                snr_grp.append(self.sources[sourceFluxField+'_snr'][i])
                #testi.append(i)
            else:
                # c'est ici qu'il faut faire les choses
                # data_titles=['visit','ccd','coord_ra','coord_dec','MJD-OBS',sourceFluxField+'_snr',sourceFluxField+'_flux',sourceFluxField+'_fluxSigma',sourceFluxField+'_mag',sourceFluxField+'_magerr','Nb_group','id','PSF-FWHM','FLUXMAG0','FLUXMAG0ERR', 'MeanGrpRa', 'MeanGrpDec', 'MedGrpSnr']

                self.RA_mean.append(np.mean(RA_grp))
              #  print'ra meab',  np.mean(RA_grp),self.sources['MeanGrpRa'][i-1]

                self.posRMS.append(np.sqrt((np.std(  RA_grp))**2 + (np.std( Dec_grp))**2))
     
                self.mag_mean.append(np.mean(mag_grp))

                self.snr_med.append(np.median(snr_grp))
            #    print 'test', test
            #    print 'testi', testi

                # initialisation du groupe suivant
                RA_grp = []
                Dec_grp = []
             
                mag_grp = []
                snr_grp = []
               # test = []
               # testi = []
                RA_grp.append(self.sources['coord_ra'][i])
                Dec_grp.append(self.sources['coord_dec'][i])
               # test.append(self.sources['Nb_group'][i])
                mag_grp.append(self.sources[sourceFluxField+'_mag'][i])
                snr_grp.append(self.sources[sourceFluxField+'_snr'][i])
             #   self.sources['Nb_group']
        self.RA_mean = np.array(self.RA_mean)
        self.Dec_mean = np.array(self.Dec_mean)
        self.posRMS =  np.array(self.posRMS)
        self.mag_mean = np.array(self.mag_mean)
        self.snr_med = np.array(self.snr_med)


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


def fitAstromErrModel(snr, dist): # fonction originaire de check.py
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

def astromErrModel(snr, theta=1000, sigmaSys=10, C=1, **kwargs): # fonction originaire de check.py
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


class Validation_plots:
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
        plt.close(fig)
       


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


        

if __name__=="__main__":
  
  #  AAA=Comparaison_multiple()
  #  AAA.plotPos()

    LDV = LoadDataValidation(file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/sources_ndarray_grp_16visits.pkl')
    LDV.createVariablesForPlots()
    print 'Data Loaded'

PLTV=Validation_plots(LDV.sources)
# PLTV.plotVisitVsTime() # plot fonctionnel


  
PLTV.plotAstrometry( LDV.posRMS*radToMas, LDV.mag_mean, LDV.snr_med)

print 'programme termine'
