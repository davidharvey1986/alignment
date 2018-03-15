

from lensing import lenstool as lt
import numpy as np
import lensing as lens
from CutGalaxiesWithSim import rec2pot as rec2pot
import ipdb as pdb
import os as os

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
    
    '''

    lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

    os.system('mkdir '+lenstoolDir+'/'+constrain)

    NewLesntoolDir = lenstoolDir+'/'+constrain
    runmode, potentiels, limits, image, potfile = \
      lt.read_par( 'bestopt.par', \
                       par_dir = lenstoolDir,\
                       return_all=True,)

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
    image['arcletstat']['int'] = 0
    lt.write_par( NewLesntoolDir+'/'+cluster+'_'+constrain+'.par', \
                      potentiel=totalPots, runmode=runmode,
                      potfile=potfile, limit=inParFileLims, image=image,
                      source=source, grille=grille)
    
    
def getDarkPot( iPot, constrain='ellipticite' ):
    '''
    Calulcate the dark matter potential 
    '''

    iDark = lt.potentiels.piemd()
    iLimit = lt.limits.piemd()
    
    for i in iDark.keys():
        metaKey = iDark[i].dtype.names[1]
        iDark[i][metaKey] = iPot[i]

    if constrain == 'ellipticite':
        iLimit['e1']['int'] = 1
        iLimit['e2']['int'] = 1
        
    return iDark, iLimit
        
    
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

     #work out the FWHM in parsecs
    r_cut_pc = iBaryons['cut_radius']['float']  / 206265. * \
      lens.ang_distance( iPot['z_lens'] ) *1e6
      
    G = 4.302e-3
    #Work out the velocity dispersion
    v_disp = np.sqrt(G*10**iPot['MSTAR'] / ( r_cut_pc  * np.pi ))
    iBaryons['v_disp']['float'] = v_disp
    iBaryons['core_radius']['float'] = 0.

    return iBaryons
