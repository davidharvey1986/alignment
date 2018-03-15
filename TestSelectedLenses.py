
import getImageSensitivity as GIS
from offset import *
from CutGalaxiesWithSim import *
from matplotlib import pyplot as plt

def TestSelectedLenses( cluster, potentials, constrain, \
                            sourceList=None,\
                            FiducialImageList=None,
                            nRuns=20, retest=False):

    '''
    Test the selected lenses to see how much the images move

    '''
    if sourceList is None:
        sourceList = cluster+'_sources.dat'

    if (FiducialImageList is None):
        FiducialImageList = cluster+'_fid.img'

    TestPots = []
    for iRecPot in potentials:
        iPot = rec2pot( iRecPot )
        if (iRecPot['ConstrainFlag'] >= 10)  & \
          (iRecPot['GalaxyFlag'] == 1):
            baryons, darkMatter = \
              split_potential( iPot )
              
            TestPots.append( darkMatter )
            TestPots.append( baryons )
            
        elif (iRecPot['ConstrainFlag'] == 1) & \
          (iRecPot['GalaxyFlag'] == 1):
            baryons, darkMatter = \
              split_potential( iPot )

            #make the darkmatter have a identity
            #not the same
            darkMatter['identity']['str'] = \
              str(darkMatter['identity']['str'])[0:-3]

            TestPots.append( darkMatter )
            TestPots.append( baryons )
        else:
            TestPots.append( iPot )

    pklName= cluster+'_'+constrain+'_testimage.pkl'
    ImageSensitivity = \
      GIS.getImageSensitivity( cluster, sourceList, FiducialImageList,
                                   constrain,
                                   nRuns=nRuns, potentials=TestPots,\
                                   pklName=pklName, \
                                   rerun=retest )

    #plotTest(potentials, ImageSensitivity)

    return ImageSensitivity

def plotTest( potentials, images):

    
    
    for iPot in potentials:
        
        closestImageDelta = np.max(images[iPot['ClosestImage']]['Delta'])

        
        plt.plot( iPot['ClosestImageDist']/iPot['cut_radius'], closestImageDelta, '.')

    plt.ylim(0,2.0)
    plt.show()
                                   

    
