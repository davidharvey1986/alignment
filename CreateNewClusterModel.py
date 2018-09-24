
from lensing import lenstool as lt
import offset as offset
import os as os
from copy import deepcopy as cp

def CreateNewClusterModel( cluster, sourceList, constrain, \
                               fiducialModel='splitPotential.par',\
                               lenstoolDir=None,
                               imageList='image_tmp.dat',
                                p_offset=0., \
                                e_offset=0., \
                                a_offset=0., \
                               m_offset=0.,
                               potentials=None):
    '''
    Take a cluster model, with baryonic and dark matter
    componenets, rotate, stretch and move the dark matter
    potential, write it to a par file and 
    run it.

    This will require the splitPotential.par file to read in
    from the previous script. I.e. the fiducial cluster model.


    '''


    if lenstoolDir is None:
       lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

      
    runmode, pots = \
        lt.read_par( fiducialModel, par_dir=lenstoolDir)   
    champ = \
        lt.read_par( fiducialModel, par_dir=lenstoolDir, return_champ=True)
    
    if potentials is None:
        potentials = cp(pots)
        
    newPotentialList = []
    for iPot in potentials:
        if 'dm' in str(iPot['identity']['str']):
            if constrain == 'ellipticite':
                newDM = offset.ellipticity( iPot, 
                        e_offset )
            if constrain == 'position':
                newDM = offset.position( iPot, \
                        p_offset )
            if constrain == 'angle_pos':
                newDM = offset.angle_pos( iPot, 
                        a_offset )
                                                



            newPotentialList.append(newDM)
        else:
                        
            newPotentialList.append(iPot)

    source = lt.source()
    source['source']['int'] = 0

    image = lt.image()
    image['image']['int'] = 0

    runmode['image']['int'] = 0
    runmode['source']['int'] = 1
    runmode['source']['filename'] = sourceList
    runmode['verbose']['int'] = 0
    runmode['mass']['int'] = 0.
    lt.write_par(lenstoolDir+'/iOffset.par', \
                      runmode=runmode,
                     potentiel=newPotentialList,
                     image=image,
                     source=source, \
                     champ=champ)

    lt.run( 'iOffset.par', outdir=lenstoolDir)

    lt.source_to_wcs(filename=lenstoolDir+'/image.all',
                  outfile=lenstoolDir+'/'+imageList)

    
    
