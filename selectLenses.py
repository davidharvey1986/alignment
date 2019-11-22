
import HFFstellarMass as HSM
import lensing as lensing
import numpy as np
from numpy.lib.recfunctions import append_fields as append_rec
import ipdb as pdb
import CutGalaxiesWithSim as CutGalaxies
import ConstructParFile as CPF
import TestSelectedLenses as TSL
from VerifyParFile import *

def main( cluster, constrain, nConstrain=30, retest=False,
              distCut=0.3, nCutRadii=4, MLCut=10., ellCut=0.):
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

    AllMemberMasses[1].data = MassToLightRatio( AllMemberMasses[1].data )
    
    MLRatioCat = AllMemberMasses[1].data

    #There seems to be a lot of galaxies to constrain so lets
    #now go to the second method
    nGalaxies = len(MLRatioCat[ MLRatioCat['GalaxyFlag'] == 1])

    print 'THERE ARE A TOTAL OF '+str(nGalaxies)+\
      ' POTENTIAL GALAXIES TO CONSTRAIN'

    AllMemberMasses[1].data = \
      CutLenses( MLRatioCat, ellCut = ellCut, MLCut = MLCut )
      
    FirstCutPots=  AllMemberMasses[1].data
    nFirstCut = len(FirstCutPots[ (FirstCutPots['ConstrainFlag'] > 0) & \
                        (FirstCutPots['GalaxyFlag'] == 1)])
    print 'AFTER INITAL CUTS WE HAVE '+str(nFirstCut)+\
      ' GALAXIES TO CONSTRAIN'

    if nFirstCut > nConstrain:   
        #This will return the potentials array with extra flags
        #on the 
        AllMemberMasses[1].data = \
          CutGalaxies.CutGalaxiesWithSim( cluster, AllMemberMasses, \
                                              constrain,\
                                              nCutRadii=nCutRadii, \
                                              distCut=distCut,\
                                              rerun=False)
        SecondCutPots = AllMemberMasses[1].data
        
        nSecondCut = \
          len(SecondCutPots[ SecondCutPots['ConstrainFlag'] >= 10])
        print 'AFTER SECOND CUTS WE HAVE '+str(nSecondCut)+\
          ' GALAXIES TO CONSTRAIN'
    else:
        
        FirstCutPots['ConstrainFlag'][FirstCutPots['ConstrainFlag'] == 1] = 10.
        AllMemberMasses[1].data = FirstCutPots
        nSecondCut = nFirstCut

    #Having appenSecondCutPotsded some constrain flags, construct the par file
    #This is where i will cut out the galaxies that cant be constrainted
    #i.e. ConstrainFlag < 10
    if nSecondCut > 0:
        TestImageSensitivity = \
        TSL.TestSelectedLenses( cluster, AllMemberMasses,  constrain,
                                    nRuns=20, retest=retest)
    else:
        TestImageSensitivity ='No potentials constrained so the sensitivity is irrelant'
    CPF.ConstructParFile( cluster, AllMemberMasses, constrain=constrain )
    VerifyParFile( cluster, constrain=constrain)
    
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
    rCutPc = memberCat['cut_radius_kpc']*1e3
    HaloMass = np.pi / G * memberCat['v_disp']**2*\
      rCutPc

    
    MLRatio =  HaloMass/ 10**memberCat['MSTAR']
    MLRatio[ memberCat['MSTAR'] == -1 ] = -1
    
    memberCat = append_rec(memberCat, 'MLRatio', MLRatio, \
                               usemask=False, asrecarray=True)

    return memberCat
