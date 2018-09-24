'''
This main will loop through each of the clusters and determine
which galaxies i want to constrain

'''
import selectLenses as sl

def main( constrain='ellipticite', nConstrain=30):
    clusters = ['A1063', 'A2744', 'A370', 'MACSJ1149']
    
    for iCluster in clusters:

        lensSensitivity = \
          sl.main( iCluster, constrain, nConstrain=30, retest=True,
                    distCut=0.1, ellCut=0.3, nCutRadii=4, MLCut=3.)

        keys = lensSensitivity.keys()
        MeanSensitivity = \
          [ np.mean(lensSensitivity[i]['Delta']) for i in keys]

        plt.hist( MeanSensitivity, bins=np.linspace(0,2,20), \
                      label=iCluster )


    plt.legend()



if __name__ == '__main__':
    main()
