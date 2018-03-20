
import HFFstellarMass as HSM
import lensing as lensing
import numpy as np
from numpy.lib.recfunctions import append_fields as append_rec
import ipdb as pdb
import CutGalaxiesWithSim as CutGalaxies
import ConstructParFile as CPF
import TestSelectedLenses as TSL

def main( cluster, constrain, nConstrain=30, retest=False,
              distCut=0.3, nCutRadii=4):
    '''
    Aim: To get the galaxy members we want to constrain, with a max
    number to constrain of nConstrain (default 30)

    #1st selection process
      1. Get the stellar mass of each halo, and determine
         the mass to light ratio of each galaxy
      2. Get the ellipticity of each halo and cut > 0.2

    #Second pass of selection process
    Method:
      1. Take the best fit model. Split the galaxies in to 
         baryonic and dark matter components.
      2. Take this model and project the arcs back to the 
         source position to a get a source list.
      3. Take the galaxies in the model and rotate / move / stretch
         run source positions through the model
      4. Measure the offset beyween the new image positions and
         the fiducial images
    '''
    
    
    AllMemberMasses = HSM.HFFstellarMass( cluster )
    
    MLRatioCat = MassToLightRatio( AllMemberMasses )

    #There seems to be a lot of galaxies to constrain so lets
    #now go to the second method
    nGalaxies = len(MLRatioCat[ MLRatioCat['GalaxyFlag'] == 1])

    print 'THERE ARE A TOTAL OF '+str(nGalaxies)+\
      ' POTENTIAL GALAXIES TO CONSTRAIN'

    FirstCutPots = CutLenses( MLRatioCat, ellCut = 0., MLCut = 0. )
      
    nFirstCut = len(FirstCutPots[ FirstCutPots['ConstrainFlag'] > 0])
    print 'AFTER INITAL CUTS WE HAVE '+str(nFirstCut)+\
      ' GALAXIES TO CONSTRAIN'

    if nFirstCut > nConstrain:   
        #This will return the potentials array with extra flags
        #on the 
        SecondCutPots = \
          CutGalaxies.CutGalaxiesWithSim( cluster, FirstCutPots, constrain,\
                                              nCutRadii=nCutRadii, \
                                              distCut=distCut)
        
        nSecondCut = \
          len(SecondCutPots[ SecondCutPots['ConstrainFlag'] >= 10])
        print 'AFTER SECOND CUTS WE HAVE '+str(nSecondCut)+\
          ' GALAXIES TO CONSTRAIN'
    else:
        SecondCutPots = FirstCutPots

    #Having appended some constrain flags, construct the par file
    #This is where i will cut out the galaxies that cant be constrainted
    #i.e. ConstrainFlag < 10
    TestImageSensitivity = \
      TSL.TestSelectedLenses( cluster, SecondCutPots,  constrain,
                                  nRuns=20, retest=retest)
                            
    CPF.ConstructParFile( cluster, SecondCutPots, constrain=constrain )

    return TestImageSensitivity



def CutLenses( cat, ellCut = 0.2, MLCut = 1./0.3 ):
    '''
    Apply the two initial cuts to the catalogue
    1. Ellipciity > 0.2
    2. ML ratio > 10.
    '''

    flags = np.ones( len( cat['ra'] ))
   
    flags[ cat['MLRatio'] < MLCut ] = 0
    flags[ cat['ellipticite'] < ellCut ] = 0

    
      
    cat = append_rec(cat, 'ConstrainFlag', flags, \
                         usemask=False, asrecarray=True)

    return cat

def MassToLightRatio( memberCat ):
    '''
    From the stellar mass and the lenstool fitted paramteres
    determine the stellar to mass ratio
    '''

    #Halo Mass for a piemd taken from
    #https://arxiv.org/pdf/astro-ph/0405607.pdf
    MLRatio = np.zeros(len(memberCat['RA'])) - 1.
    
    G = 4.302e-3
    rCutPc = memberCat['cut_radius']/206265.*\
      lensing.ang_distance( memberCat['z_lens'][0] )*1e6
    HaloMass = np.pi / G * memberCat['v_disp']**2*\
      rCutPc

    
    MLRatio =  HaloMass/ 10**memberCat['MSTAR']

    memberCat = append_rec(memberCat, 'MLRatio', MLRatio, \
                               usemask=False, asrecarray=True)
    
    return memberCat
