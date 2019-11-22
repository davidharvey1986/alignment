
import numpy as np
import ipdb as pdb
import astro_tools_extra as at
import lensing as lens
from copy import deepcopy as cp

def ImageWithinRadius( SensitivityDict, potentials, nConstrain,
                           nCutRadii=3, distCut=0.3):
    '''
    This script will take the potentials and the sensitivy dictionary
    and find those lenses that have iamgs within N*Cut Radius.
    
    If this is many, then it will find those that lie within 
    the mosts sensitive images

    1. Is it a galaxy? (as apposed form a large scale halo)
    2. Does it pass the M/L and ellipcitie cut?
    3. Does it have an IMage wihtin the radius oft he halo?

    Key 0 : It has a too low ellipticity or M/L ratio AND it doesnt
        have any multiple image near it. The worst
        1 : It has a okay ellipticity and M/L but has no multiple
        image within its radius
        10: It has a okay ellipticity and M/L and has a multiple
        image within its radius

    KEYWRODS :
        nCutRadii : Number of cut radii that can we search out to
                    for a image. Default = 2
       

    '''
    
    OldPotentials = cp( potentials )
    ImageID = SensitivityDict.keys()
    ImageRA = np.array([ SensitivityDict[i]['RA']  for i in ImageID ])
    ImageDEC = np.array([ SensitivityDict[i]['DEC'] for i in ImageID ])
    MeanSensitivity = \
      np.array([ np.mean(SensitivityDict[i]['Delta']) for i in ImageID ])
    
    kpc2arcsec = 206265./lens.ang_distance(potentials['z_lens'][-1])/1e3
    CutRadiusImageSep = (potentials['ClosestImageDist'] - \
            nCutRadii*potentials['cut_radius_kpc']*kpc2arcsec < 0)

        
    
    ImagesInRadius = potentials['ClosestImage'][CutRadiusImageSep]

    iPotSensitivity = np.array([])

    for iPot in potentials:
        
        #for each image get the lenses that are within nRadii
        #of that lens

        ImageDist = at.ra_separation( iPot['ra'], iPot['dec'], \
                                          ImageRA, ImageDEC, \
                                          abs=True )
                                          
        InRadius = ImageDist < iPot['cut_radius_kpc']*nCutRadii*kpc2arcsec
        
        
        if len(MeanSensitivity[ InRadius ]) == 0:
            continue
        if np.max(MeanSensitivity[ InRadius ]) >= distCut:
            
            iPot['ConstrainFlag'] *= 10**len(MeanSensitivity[ InRadius ])   
            if iPot['ConstrainFlag'] > 0:
                iPotSensitivity  = \
                  np.append( iPotSensitivity,  \
                                 np.max(MeanSensitivity[ InRadius ]))
    SecondCutConstrain = \
      len( potentials[ potentials['ConstrainFlag'] >= 10 ])
    
    if SecondCutConstrain > nConstrain:
        #Find the distcut that would remove all but nConstrain gals
        
        NewDistCut = np.sort(iPotSensitivity)[::-1][nConstrain]
        NewNConstrain = len(iPotSensitivity[iPotSensitivity>=NewDistCut])

        print('WARNING: Number to be constrained is large, using new DistCut of %f and new nConstrain %i' % (NewDistCut,NewNConstrain))
        potentials = \
          ImageWithinRadius( SensitivityDict, OldPotentials, \
                                 NewNConstrain, nCutRadii=nCutRadii, distCut=NewDistCut)

    return potentials



    

def ClosestLenstoImage( SensitivityDict, potentials, nConstrain,
                            distCut=0.3):
    '''
    This will find the lenses taht are closest to the multiple iamges
    SO the cut will happend

    1. Is it a galaxy? (as apposed form a large scale halo)
    2. Does it pass the M/L and ellipcitie cut?
    3. Is it the closest lens to a multiple iamge?
    4. Does that image move more than 0.5 arcseconds?


    So if i add 10 to each constrainable galaxy, ConstrainFlag
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
    
    '''

    meanDelta = [ np.mean(SensitivityDict[i]['Delta']) \
                      for i in SensitivityDict.keys() ]
                      
    ClosestLens = np.array([ SensitivityDict[i]['ClosestLens'] \
                        for i in SensitivityDict.keys() ])

    ClosestLensDelta = np.array([ np.max(SensitivityDict[i]['Delta']) \
                                      for i in SensitivityDict.keys() ])
    
    for i, iLens in enumerate(ClosestLens):
        if ClosestLensDelta[i] < distCut:
            continue
        
        if not iLens in potentials['identity']:
            raise ValueError('%s not in potentials' % iLens)
        
        IDindex = potentials['identity'] == iLens
        potentials['ConstrainFlag'][IDindex] *= 10


    SecondCutConstrain = \
      len( potentials[ potentials['ConstrainFlag'] >= 10 ])

    if SecondCutConstrain > nConstrain:
        raise ValueError('Need another cut in Image Sensitivity')
 

    return potentials
