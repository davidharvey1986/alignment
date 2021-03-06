'''
This main will loop through each of the clusters and determine
which galaxies i want to constrain

'''
import selectLenses as sl
import numpy as np
from matplotlib import pyplot as plt
def main( constrain='ellipticite', nConstrain=30):
    clusters = ['A1063', 'A2744', 'A370', 'MACSJ1149']
    
    for iCluster in clusters:
        print iCluster
        lensSensitivity = \
          sl.main( iCluster, constrain, nConstrain=30, retest=False, 
                    distCut=0.1, ellCut=0.2, nCutRadii=3, MLCut=1.)

        keys = lensSensitivity.keys()
        MeanSensitivity = \
          [ np.mean(lensSensitivity[i]['Delta']) for i in keys]

        plt.hist( MeanSensitivity, bins=np.linspace(0,2,20), \
                      label=iCluster, alpha=0.5 )


    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
