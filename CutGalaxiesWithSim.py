import offset as offset
from lensing import lenstool as lt
import os as os
import astro_tools_extra as at
import numpy as np
import ipdb as pdb
import getImageSensitivity as GIS
import getClosestLens as GCL
import getConstrainIDs as getConstrainIDs

def CutGalaxiesWithSim( cluster, potentials, constrain, \
                            nCutRadii=5.,\
                            distCut=0.3, \
                            nConstrain=30 ):
    '''
    #This script will run the second pass of selection process
    Method:
      1. The input is a list of potentials that can be
      written to a par file. 
      The 'ConstrainFlag' in these potentials represent
      each galaxy that COULD be constrained.
      I need to whittle down these ones to <20 galaxies

      1. Create a source list from the best.par
      2. Split them into baryon and dark matter components
         and create a fiducial image list
      4. Rotate the galaxies / stretch them / move them
         and reproject the sources back to the image plane
      5. Compare the image position before and after

    The final potentials will have the ConstrainFlag with the
    following meaning
    will have the following meaning
    0 : It has a too low ellipticity or M/L ratio AND it doesnt
        have any multiple image near it. The worst
    1 : It has a okay ellipticity and M/L but has no multiple
        image near it. Not good enough.
    10 : It has a okay ellipticity and M/L AND 
         it has a 1 multiple image near it. Good.
    100 : It has a okay ellipticity and M/L AND 
         it has a 2 multiple images near it. Very Good.
    1000 : It has a okay ellipticity and M/L AND 
         it has a 3 multiple image near it. Best.
    
    INPUTS : 
      CLUSTER : the names of the cluster
      POTENTIALS : a recarray of the potentials int eh clsuter
      CONSTRAIN : the paramter i want to constrain

    KEYWORDS : 
       nCutRadii : how many radii does an image need to be within to count
       distCut : how much does a multiple image have to move to counted
                 as significant movement.     
       nCOnstrain : Max number i want to constrain
    '''

    #1. Get the source positions
    sourceList = cluster+'_sources.dat'
    createSourceList( cluster, sourceList=sourceList )
    #2.Create the fiducial imageList
    FiducialImageList = cluster+'_fid.img'
    createImageList( cluster, potentials, sourceList, \
                         imageList=FiducialImageList )

    #3.Now determine how sensitive each image is to the model
    # by creating a new cluster model with rotated and shifted
    #  DM halos and creating a new image list
    #This returns a dictionary of the images with how sensitive they are:
    #I.e. the mean radial shift in the image
    imageSensitivity = GIS.getImageSensitivity( cluster, \
                                            sourceList, \
                                            FiducialImageList, \
                                            constrain, \
                                            nRuns=20)
    
    #Now append to the dictionart of images which lens is closest
    LensSensitivity = \
      GCL.getClosestLens( imageSensitivity, potentials )

    PotentialsWithClosestImage = \
      GCL.getClosestImage(  imageSensitivity, potentials )

      
    ConstrainPotentials = \
      getConstrainIDs.ImageWithinRadius( LensSensitivity, \
                                             PotentialsWithClosestImage,\
                                             nConstrain, \
                                             nCutRadii=nCutRadii,
                                             distCut=distCut)

    #ConstrainPotentials = \
     # getConstrainIDs.ClosestLenstoImage( LensSensitivity, \
      #     PotentialsWithClosestImage, \
       #     nConstrain)
    
    return ConstrainPotentials



def createSourceList( cluster, \
                          sourceList='source.dat', \
                          lenstoolDir = None):

    '''
    Create a fiducial source list from the best.par

    '''
    if lenstoolDir is None:
       lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

    if os.path.isfile( lenstoolDir+'/'+sourceList):
        return
        
 
       
    runMode, potentials = \
      lt.read_best( lenstoolDir+'/best.par')
    image = \
      lt.read_par( 'bestopt.par', par_dir=lenstoolDir, return_image=True)
      
    runMode['source']['int'] = 0
    runMode['image']['filename'] = image['multfile']['filename']
    runMode['image']['int'] = 1
    runMode['inverse']['int'] = 0
    runMode['restart']['int'] = 0
    
    source = lt.source()
    source['source']['int'] = 0

    image['image']['int'] = 0
   
    
    lt.write_par( lenstoolDir+'/CreateSources.par',
                     runmode=runMode,
                     potentiel=potentials,
                     image=image,
                     source=source)

    lt.run( 'CreateSources.par', outdir=lenstoolDir)

    os.system('mv '+lenstoolDir+'/source.dat '+\
                  lenstoolDir+'/'+sourceList)

    
    
    
def createImageList( cluster, potentials, sourceList,
                         lenstoolDir=None, imageList='image.dat'):
    '''
    Using the potentials, create a fiducial Image list
    that i can compare to 

    2. However, for each of the sources there will be many images
    so potentially 3 times as many imags as the original image list
    Therefore I will go through the new image list, find for 
    each image label the closest image to the original image
    and use that as the fiducial imagelist

    '''
    if lenstoolDir is None:
       lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster


    if os.path.isfile( lenstoolDir+'/'+imageList):
        return
    
    potentialList = []
    for iPot in potentials:
        #convert the recarray entry to a potential for lt scripts
        newPot = rec2pot( iPot )
        
        if (iPot['ConstrainFlag'] == 1) & \
          (iPot['GalaxyFlag'] == 1):
        
         
          baryons, darkMatter = offset.split_potential( newPot )

          potentialList.append( baryons)
          potentialList.append( darkMatter )
            
        else:

            potentialList.append( newPot )



    runmode, potentials = lt.read_par('bestopt.par', \
                                          par_dir=lenstoolDir)
    runmode['inverse']['int'] = 0
    runmode['source']['int'] = 1
    runmode['source']['filename'] = sourceList
    runmode['image']['int'] = 0

    source = lt.source()
    source['source']['int'] = 0
    image = lt.image()
    image['image']['int'] = 0
    

    lt.write_par(lenstoolDir+'/splitPotential.par',
                      runmode=runmode,
                     potentiel=potentialList,
                     image=image,
                     source=source)
    lt.run( 'splitPotential.par', outdir=lenstoolDir)

    #2.
    cleanImageList( lenstoolDir+'/'+runmode['image']['filename'],\
                        lenstoolDir+'/image.dat',
                        lenstoolDir+'/'+imageList)


def cleanImageList( fiducialList, newList, outputList):
    '''
    Take the new list and compare to the fiducial list
    by
    1. looping through all the image id's, finding the closest
       image to the fiducial image, and using that as the newimage
    
    '''
    FidCatHeader = open(fiducialList, 'rb').readline()
    
    
    dtypes = [('ID', '|S20'), ('RA', float), ('DEC', float), \
                  ('SemiMajor', float), ('SemiMinor', float), \
                  ('Angle', float), ('Redshift', float), ('Mag', float)]
    
    FidCat = np.loadtxt( fiducialList, dtype=dtypes)
    NewCat = np.loadtxt( newList, dtype=dtypes)

    outCat = open(outputList, 'wb')
    outCat.write( FidCatHeader )
    
    for iImage in FidCat:

        IDindex =  NewCat['ID'] == iImage['ID']
        NewCatRA = NewCat['RA'][ IDindex ]
        NewCatDEC = NewCat['DEC'][ IDindex ]

        Separation = at.ra_separation( NewCatRA, NewCatDEC, \
                                           iImage['RA'], iImage['DEC'], 
                                           abs=True)
        
        outCat.write('%s %0.7f %0.7f %0.1f %0.1f %0.1f %0.5f %0.1f \n' %  \
                         tuple(NewCat[IDindex][np.argmin(Separation)]))


    outCat.close()

        

        
        
        
        
    



def rec2pot( recEntry):
    '''
    Take a recarray entry and convert it into a potential
    that can be used in my lenstool scripts
    '''
    newPot = lt.potentiels.get_profile(recEntry['profil'])
    
    #Because the input potentials is not
    #int he standard format i need to readjust
    for iName in recEntry.dtype.names:
        if not iName in newPot.keys():
            continue
        metaKey = newPot[iName].dtype.names[-1]
        newPot[iName][metaKey] = recEntry[iName]
        
                         
    return newPot
