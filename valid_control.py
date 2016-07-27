import pylab as plt
import numpy as np

import scipy.stats
import pickle

#from ..kdtree import KDTree

sourceFluxField='base_PsfFlux'

color = {'all': 'grey', 'bright': 'blue',
         'iqr': 'green', 'rms': 'red'}

radToDeg = 180./np.pi
degToArcs = 3600.
radToArcs = radToDeg * degToArcs


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

 
class Validation:
    def __init__(self, file_path='/sps/lsst/dev/ciulli/Validation_lsst_lpc_git/test_validation/sources_ndarray_grp_16visits.pkl'):
        self.sources = load_pkl(file_path)
        

    def createVariablesForPlots(self):
        self.dist=[] #posRMS
        self.mag=[] #sourceFluxField
        self.snr=[] # group med snr
       # for i, j in enumerate(self.sources['Nb_group']):
          #  if self.sources['Nb_group'][i] == j:
          #      RA_grp.append(self.sources['coord_ra'][i]

        #for i in range(len(self.sources['Nb_group'])): 
        RA_grp = []
        Dec_grp = []
        test = []
        testi = []
        for i, j in enumerate(self.sources['Nb_group']):
            print 'i, j', i ,j
            if j == self.sources['Nb_group'][i+1]:
                RA_grp.append(self.sources['coord_ra'][i])
                Dec_grp.append(self.sources['coord_dec'][i])
                test.append(self.sources['Nb_group'][i])
                testi.append(i)
            else:
                RA_grp.append(self.sources['coord_ra'][i])
                Dec_grp.append(self.sources['coord_dec'][i])
                test.append(self.sources['Nb_group'][i])
                testi.append(i)
                print 'test', test
                print 'testi', testi
                RA_grp = []
                Dec_grp = []
                test = []
                testi = []

             #   self.sources['Nb_group']
     
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
     
    print 'programme termine'

    AAA=Comparaison_multiple()
    AAA.plotPos()

    BBB=Validation()
    BBB.plotVisitVsTime()
