'''

This scipt is designed to verify the par file that is constructed
from the ConstructParFile

It will compare the the expected image positions from these
with the true image positions

If it is massive we know it is wrong

'''

from  CompareImageLists import *
from CutGalaxiesWithSim import *
import os as os
def VerifyParFile( cluster, constrain='ellipticite'):
    '''
    Take the parfile, make the inverse =0 and then the output
    image.dat and compare to the the multi image file

    UPDATE: Dont use the old sources, just run the multiple images
    backwards and then forwards
    '''
    baseDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster
    lenstoolDir = baseDir+'/'+constrain

    runmode, potentiels, limits, image, source, potfile, grille, champ = \
      lt.read_par( cluster+'_'+constrain+'.par', \
                       par_dir = lenstoolDir,\
                       return_all=True)

    runmode['inverse']['int'] = 0
    runmode['image']['int'] = 0
    #runmode['image']['filename'] = cluster+'_mult.img'
    runmode['source']['filename'] = '../'+cluster+'_sources.dat'
    if source is None:
        source = lt.source()
        source['source']['int'] = 0
    image['arcletstat']['option'] = 0
    runmode['source']['int'] = 1
    runmode['verbose']['int'] = 0
    grille['nlens']['int'] = len( potentiels)
    grille['nlens_opt']['int'] = len( limits)
    #image['multfile']['filename'] = cluster+'_mult.img'


    newParFile = lenstoolDir+'/'+cluster+'_'+constrain+'_verify.par'
    nLensOpt = len(limits)
    nLens = len(potentiels)



    
    lt.write_par( newParFile, \
                      potentiel=potentiels[0:len(limits)], runmode=runmode,
                      potfile=potfile, limit=limits, image=image,
                      source=source, grille=grille, champ=champ)

    lt.run( cluster+'_'+constrain+'_verify.par', outdir=lenstoolDir)

    #this creates a list of images that the par file will produce
    #in the zeroth order.
    
    cleanImageList( baseDir+'/'+image['multfile']['filename'],\
                         lenstoolDir+'/image.dat',
                        lenstoolDir+'/verify_images.dat')
    
    imageDifferences = \
       CompareImageLists(lenstoolDir+'/verify_images.dat', \
                               baseDir+'/'+image['multfile']['filename'])

    
    for i in imageDifferences:
        print("%s %0.5f" % (i['ID'],i['Delta']))
    print("The distances are %0.5f" %  np.mean(imageDifferences['Delta']))
