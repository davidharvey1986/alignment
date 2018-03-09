
import ipdb as pdb
import numpy as np
import astro_tools_extra as at
from numpy.lib.recfunctions import append_fields as append_rec

def getClosestLens( ImageSensitivity, potentials ):
    '''
    Loop through each ID in the ImageSensitivity Dictorionary
    and get the closest lens in the potentials
    '''

    #Get rid of the potentials that are massive
    #i.e. keep only the galaxy ones
    Index = np.arange(len(potentials))[ potentials['GalaxyFlag'] == 1]

    for iImage in ImageSensitivity.keys():

        ImageLensSeparation = \
          at.ra_separation( ImageSensitivity[iImage]['RA'], \
                            ImageSensitivity[iImage]['DEC'], \
                            potentials['ra'][Index], \
                            potentials['dec'][Index], abs=True)


        ImageSensitivity[iImage]['ClosestLens'] = \
          potentials['identity'][Index][np.argmin(ImageLensSeparation)]
                                

    return ImageSensitivity

def getClosestImage( ImageSensitivity, potentials ):
    '''
    Loop through each potential in the potentials recarray
    and get the closest image to that drom the ImageSensitivity dict
    '''

    #Get rid of the potentials that are massive
    #i.e. keep only the galaxy ones
    Index = np.arange(len(potentials))[ potentials['GalaxyFlag'] == 1]

    ImageID = ImageSensitivity.keys()
    ImageRA = np.array([ ImageSensitivity[i]['RA']  for i in ImageID ])
    ImageDEC = np.array([ ImageSensitivity[i]['DEC'] for i in ImageID ])
    

    ClosestImage = []
    ClosestImageDist = []
    
    
    for iPot in potentials:

        ImageLensSeparation = \
          at.ra_separation( ImageRA, ImageDEC, \
                            iPot['ra'], \
                            iPot['dec'], abs=True)


                        
        ClosestImage.append( ImageID[np.argmin(ImageLensSeparation)])
                
        ClosestImageDist.append(np.min(ImageLensSeparation))

    potentials = append_rec(potentials, 'ClosestImage', \
                                ClosestImage, usemask=False, \
                                asrecarray=True)
    potentials = append_rec(potentials, 'ClosestImageDist', \
                                ClosestImageDist, usemask=False, \
                                asrecarray=True)
    
    return potentials
