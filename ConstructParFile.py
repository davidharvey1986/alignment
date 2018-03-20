

from lensing import lenstool as lt
import numpy as np
import lensing as lens
from CutGalaxiesWithSim import rec2pot as rec2pot
import ipdb as pdb
import os as os
from copy import deepcopy as cp

def ConstructParFile( cluster, potentials, constrain='ellipticite' ):
    '''
    Construct the parameter file for the new reconstruction

    This needs to be based on the bestopt.par

    1. First get the potentials limits from the bestopt.par
       including the potfiles.

    2. Then get the cluster member catalogue. 
       Loop through this catalogue and
       a.) append to the potentials the stellar and dark matter pots
       b.) remove the galaxies from the cluster member catalogue

    3. Write out the file

    TO BE CAREFUL
    1. All consaining potentials need to be first in the list
       then the fixed baryons after
    2. Need to keep the first six halos as they are the big ones


    UPDATES:
    1. New idea is the fix the stellar mass to the SED mass, then for each galaxy
       in both the par file and the pot file, fix mass to the SHMR from
       Velander et al 2013. This means there will be no unphysical MLRatios.
     
    Issues: If i take the new total mass, which in a PIEMD is a function of 
            vdisp and cut radius, shold i keep the cut radius the same and just
            change the velocity dispersion?
    '''

    lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

    os.system('mkdir '+lenstoolDir+'/'+constrain)

    NewLesntoolDir = lenstoolDir+'/'+constrain
    runmode, potentiels, limits, image, potfile = \
      lt.read_par( 'bestopt.par', \
                       par_dir = lenstoolDir,\
                       return_all=True,)


    GalaxyPotfile =  [ i for i in potfile if 'member' in str(i['filein']['filename']) ]

    if len(GalaxyPotfile) > 1:
        raise ValueError('Expect exactly one potfile for the cluster member, found %i',\
                             len(GalaxyPotfile))
    else:
        GalaxyPotfile = GalaxyPotfile[0]
   
    MembersHeader = open(lenstoolDir+'/'+cluster+'_members.cat').readline()
    
    dtypes = [('ID', '|S20'), ('RA', float), ('DEC', float), \
                  ('SemiMajor', float), ('SemiMinor', float), \
                  ('Angle', float), ('Redshift', float), ('Mag', float)]
    
    clusterMembers = np.array( [], dtype = dtypes)

    
    inParFilePotsDM = [ rec2pot(i)  for i in potentials[0:len(limits)]]
    inParFilePotsBaryons = []
    inParFileLims = limits

    NewClusterMemberCat = open(NewLesntoolDir+'/new_member.cat', 'wb')
    NewClusterMemberCat.write(MembersHeader)

    nLensOpt = len(limits)
    nLens = 0
    for iPot in potentials:
        nLens+=1
        
        if (iPot['ConstrainFlag'] >= 10) & \
          (iPot['GalaxyFlag'] == 1):
            #First sort the baryons
            
            iBaryons = getBaryonPot( iPot )
            #Setting the baryonic mass to that of the SED
            #the dark matter mass is fixed at the SHMR
            #see function getDarkPot for more
            iDM, iDM_limit = getDarkPot( iPot )

            iDM['identity']['str'] = iDM['identity']['str']+'_dm'
            inParFilePotsDM.append(iDM)
            inParFileLims.append(iDM_limit)

            iBaryons['identity']['str'] = \
              iBaryons['identity']['str']+'_baryons'
            inParFilePotsBaryons.append(iBaryons)

            

            nLensOpt += 1
            #There are two lenses so add but only optimising 1
            nLens += 1
        elif iPot['GalaxyFlag'] == 1:

            #For each potfile galaxy, fix the stellar mass
            #to the SED mass, scale to the total mass
            #then set the cut radius and velocity dispersion
            #based on these from the potfile results
            iPot = FixStellarMass( iPot, GalaxyPotfile )
        
            TupleWrite = (iPot['identity'], iPot['ra'], iPot['dec'], \
                                   iPot['semiMajor'], iPot['semiMinor'],
                              iPot['angle'], iPot['z_lens'], iPot['mag'])
                                   
            NewClusterMemberCat.write('%s %0.7f %0.7f %0.1f %0.1f %f %0.5f %0.1f\n' % TupleWrite)                                  
            
        
                                          
    runmode['inverse']['int'] = 3
        
        
    NewClusterMemberCat.close()

    totalPots = inParFilePotsDM + inParFilePotsBaryons

    grille = lt.grille( nLens, nLensOpt)
    source = lt.source()
    source['source']['int'] = 0
    image['arcletstat']['option'] = 0
    lt.write_par( NewLesntoolDir+'/'+cluster+'_'+constrain+'.par', \
                      potentiel=totalPots, runmode=runmode,
                      potfile=potfile, limit=inParFileLims, image=image,
                      source=source, grille=grille)
    
    
def getDarkPot( iPot, constrain='ellipticite' ):
    '''
    Calulcate the dark matter potential. Take the stellar mass and
    then calculate the SHMR for the total mass. This will change the velocity dispersion
    assuming the cut radius is the same (which is isnt, so do i need to free this up?!)
    '''

    iDark = lt.potentiels.piemd()
    iLimit = lt.limits.piemd()
    
    for i in iDark.keys():
        metaKey = iDark[i].dtype.names[1]
        iDark[i][metaKey] = iPot[i]

    if constrain == 'ellipticite':
        iLimit['e1']['int'] = 1
        iLimit['e2']['int'] = 1

    DarkMass = Stellar2HaloMass( 10**iPot['MSTAR'], iPot['z_lens'] )

    
    oldVdisp = cp(iDark['v_disp']['float'])
    
    #Convert the old vdisp using the new total mass from the SHMR
    iDark['v_disp']['float'] = Mass2Vdisp( iPot, np.log10(DarkMass) )

    return iDark, iLimit

def Stellar2HaloMass( StellarMass, z ):
    '''
    Using the Velander et al 2013 stellar to halo mass relation 
    convert the stellar mass halo to the total mass 
    halo
    '''
    #Velander et al parameters for red galaxies
    h0 = 0.7 #Small hubble parameter
    M0 = 1.4e13/h0
        
    beta = 1.36
    Mfiducial = 2e11/h0**2   

    M200 = M0*(StellarMass/Mfiducial)**beta

    return M200

def Stellar2HaloMassLater( StellarMass, z):
    '''
    Using the Moster et al stellar to halo mass relation 
    convert the stellar mass halo to the total mass 
    halo

    https://arxiv.org/pdf/1205.5807.pdf
    '''

    M10 = 11.590
    M11 = 1.195

    N10 = 0.0351
    N11 = -0.0247

    B10 = 1.376
    B11 = -0.826

    G10 = 0.608
    G11 = 0.329
    
    N = N10 + N11*z 
    M1 = 10**(M10 + M11*z)
    beta = B10*(1+z)**(B11)
    gamma = G10*(1+z)**(G11)

    loMassIndex = 1./(1.-gamma)
    hiMassIndex = 1./(1.+beta)
                 
    M200  = (StellarMass*M1**(-gamma)/(2.*N))**loMassIndex + \
      (StellarMass*M1**(beta)/(2.*N))**hiMassIndex 

    return M200
    
def getBaryonPot( iPot ):
    '''
    Calculate all the necessary values for the Baryon potential
    '''

    iBaryons = lt.potentiels.piemd()
    for i in iBaryons.keys():
        metaKey = iBaryons[i].dtype.names[1]
        iBaryons[i][metaKey] = iPot[i]
        
    iBaryons['cut_radius']['float'] = iPot['semiMajor'] / \
      np.sqrt(1+iPot['ellipticite'])
    
    iBaryons['v_disp']['float'] = Mass2Vdisp( iPot, iPot['MSTAR'] )

    iBaryons['core_radius']['float'] = 0.
    
    return iBaryons



def Mass2Vdisp( iPot, mass ):
    '''
    Take the mass and readius and return a velocity 
    dispersion for the halo

    iPot is the potential structure and mass is the logarith of the nass
    of the halo

    '''
    
    
    #work out the FWHM in parsecs
    r_cut_pc = iPot['cut_radius']  / 206265. * \
      lens.ang_distance( iPot['z_lens'] ) *1e6
      
    G = 4.302e-3
    #Work out the velocity dispersion
    v_disp = np.sqrt(G*10**mass / ( r_cut_pc  * np.pi ))
    
    return v_disp


def FixStellarMass( iPot, potfile ):
    '''
    If i fix the stellar mass to the SED mass, then
    convert to a total mass. I then work out the 
    velocity disprsion and then the magntidue that
    corresponds to this. Then change the magnitude 

    Assuming the Faber-jackson relation of L proprto vdisp**4

    Pot is the potential structure with the variosu tags required
    of the indiviual galaxy

    Potfile is the constrained general terms from the fiduial
    run, i.e. the L*
    '''

    TotalMass = Stellar2HaloMass( 10**iPot['MSTAR'], iPot['z_lens'] )

    currentMass =  PotMass( iPot )

    newMag = getNewMag( TotalMass, potfile )
    
    
    if iPot['MSTAR'] != -1:
        iPot['mag'] = newMag

    return iPot
    
def getNewMag( mass, potfile ):
    '''
    Using the Faber and Jackson relatiom ,convert the
    sigma nack to magnotude

    '''
    G = 4.302e-3
    
    potfileMass = \
      np.pi / G * potfile['sigma']['lo']**2 * potfile['cutkpc']['lo']*1e3
   
    deltaLum = (mass / potfileMass)

    newMag = potfile['mag0']['float']-2.5*np.log10(deltaLum)

    return newMag
    



def PotMass( iPot ):
    '''
    Get the current mass of hte halo
    '''
    G = 4.302e-3

    r_cut_pc = iPot['cut_radius']  / 206265. * \
      lens.ang_distance( iPot['z_lens'] ) *1e6
      
    currentMass = np.pi / G * iPot['v_disp']**2 * r_cut_pc
    
    return currentMass
    
    
