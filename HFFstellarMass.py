
import numpy as np
from numpy.lib.recfunctions import append_fields as append_rec
import pyfits as fits
from lensing import lenstool as lt
import astro_tools_extra as at
import os as os
import ipdb as pdb
import idlsave as save
import getSherpaShapes as GSS
from CheckMultipleImageFile import *

def HFFstellarMass( cluster ):
    '''
    Get	the stellar masses for all the galaxies in ther cluster
    Using the astrodeep masses

    To do this I need to
    1. I need to get the galaxies out of the best.par using the members.cat
    2. Then correlate this with the astro deep catalogue which will give
       us the stellar masses


    NOTE:
    Change to Deepspace catalogues as these are uniform over the halos

    '''
    #All the files I will need for this script
  
    clusterDir = '/Users/DavidHarvey/Documents/Work/Mergers/data/HFF/'
    lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/clusters/'
    massDir = clusterDir+'/'+cluster
    modelDir = lenstoolDir+'/'+cluster

    massCat = massDir+'/DEEPSPACE/'+cluster+'_DEEPSPACE_masses.cat'
    massCatRADEC = massDir+'/DEEPSPACE/'+cluster+'_DEEPSPACE.cat.save'
    massRADECFits = massDir+'/DEEPSPACE/'+cluster+'_DEEPSPACE_all.fits'
    ClusterImage = massDir+'/'+cluster+'_F814W_drz_sci.fits'
    
    bestPar = modelDir+'/best.par'
    bestParFits = modelDir+'/best.fits'
    membersCat = modelDir+'/'+cluster+'_members.cat'
    membersCatSherpa = modelDir+'/'+cluster+'_SherpaMembers.cat'
    membersCatFits = modelDir+'/'+cluster+'_members.fits'
    DMcatFits = modelDir+'/'+cluster+'_members_model.fits'
    TOTcatFits = modelDir+'/'+cluster+'_members_final.fits'

    #0 Check that the member catalogue has been converted into
    #  sherpa shapes first.
    if not os.path.isfile(membersCatSherpa):
        GSS.getSherpaShapes( membersCat, ClusterImage,
                            membersCatSherpa, 
                            pixelScale=0.03)

    #1. Convert the member catalogue to fits file
    if not os.path.isfile( membersCatFits):
        clusterMembersToFits( membersCatSherpa, membersCatFits )

    #Check the multiple image file for free redshifts
    #that need to be updated to have their bestopt.par limits
    CheckMultipleImageFile( cluster )
    
    #2. Convert the par file to fits
    if not os.path.isfile(bestParFits):
        lt.best_to_fits( bestPar, outfile=bestParFits)

    #3. Match the catalogues to get the dark matter properties of the
    #   galaxies
    if not os.path.isfile(DMcatFits):
        combineMembersAndBest( membersCatFits, bestParFits, DMcatFits)

    #4. The astrodeep catalaogue does not have ra and dec in the
    #   mass file so add it
    if not os.path.isfile(massRADECFits):
        DeepSpaceAddMass( massCatRADEC, massCat, massRADECFits)
    
    #4. Match the catlogues with the astrodeep catalogue to get
    #   a final catalogue with stellar mass and all dark matter masses
    if not os.path.isfile(TOTcatFits):
        DMandSMCat = matchCat( massRADECFits, DMcatFits, TOTcatFits)

    return fits.open(TOTcatFits)

def DeepSpaceAddMass( radecFile, massCatFile, outCat):

    dtypes = [('id',int), ('z', float), ('ltau', float),\
                  ('metal', float), ('lage', float), ('Av',float), \
                  ('lmass', float), ('lsfr', float), ('lssfr', float),\
                  ('la2t', float), ('chi2', float)]

    massCat = np.loadtxt( massCatFile, dtype=dtypes)
    
    radec = save.read( radecFile )

    mass = np.array([])
    

    for iGal in radec['id']:
        index = iGal == massCat['id']
        mass = \
          np.append( mass, massCat['lmass'][ index] )


    columns = []
    for iKey in radec.keys():
        if iKey == 'readme':
            continue
        if radec[iKey].dtype == np.float:
            continue

        if radec[iKey].dtype == object:
            iFormat='A20'
        else:
            iFormat=radec[iKey].dtype

        iColumn = fits.Column( name =iKey, \
                                   array=radec[iKey],\
                                   format=iFormat)
        columns.append(iColumn)
       
    columns.append(fits.Column( name='MSTAR',
                                array=mass,
                                format=mass.dtype))

        
    TotalCat = fits.BinTableHDU.from_columns(columns)
    
    TotalCat.writeto( outCat, clobber=True)

def matchCat( catA, catB, outCat):
    '''
    Match catA with catB and change add a field, RA and DEC
    that are RA_1 and DEC_1
    '''
    #add indexes to the second (larger) catalogue, so i know which
    #ones are there after combining
    catBData = fits.open(catB)

    GalID = np.arange(len(catBData[1].data['RA']))
    
    catBData.writeto(catB, clobber=True)
    
    combinedCat = at.run_match(catA, catB)[1].data
    

    NewCols = []
    NewCols.append( fits.Column(name='ra', array=combinedCat['RA_1'], format='D'))
    NewCols.append( fits.Column(name='dec', array=combinedCat['DEC_1'], format='D'))

 
    
    if not 'GalaxyFlag' in combinedCat.columns.names:
        
        NewCols.append(fits.Column(name='GalaxyFlag', \
                                array=np.ones(len(combinedCat['DEC_1'])), \
                                       format='D'))
            
                                 
    combinedCat = fits.BinTableHDU.from_columns(combinedCat.columns + fits.ColDefs(NewCols)).data
    
    #Add back in the objects that missed out, i.e. non-galaxies
    columns = []
    
    for iField in combinedCat.dtype.names:
        
        if iField == 'identity':
            iFormat='A20'
        else:
            iFormat=np.dtype(float)

        if ('_2' in iField) | ('_1' in iField):
            continue
        

        if iField in  catBData[1].data.dtype.names:
            if catBData[1].data[iField].dtype == '|S20':
                iFormat = 'A20'
            
            iColumn = fits.Column( name =iField, \
                                       array=catBData[1].data[iField],\
                                       format=iFormat)
    
        else:
            if combinedCat[iField].dtype == '|S20':
                continue                
            array = np.zeros(len(GalID)) - 1
            
            array[combinedCat['ID_2'].astype(int)] = combinedCat[iField]
            iColumn = fits.Column( name =iField, array=array,\
                                       format=iFormat)

        
        columns.append(iColumn)
        
    TotalCat = fits.BinTableHDU.from_columns(columns)
    TotalCat = [fits.PrimaryHDU(), TotalCat]
    
    #append extra potentials in teh best.fits for
    #example external shear
    catBExtensions = fits.open( catB)
    if len(catBExtensions) >2:
        
        for iExt in xrange(len(catBExtensions)-2):
            TotalCat.append(catBExtensions[iExt+2])

    thdulist = fits.HDUList(TotalCat)
    
    thdulist.writeto(  outCat,  clobber=True )

    
def clusterMembersToFits( membersCat, outFits ):
    '''
    Convert a lenstool potfile for cluster members to a fits file

    '''
    #Get the header as
    #refInt = 0 is not relative RA and DEC
    #refInt = 3 is relative RA to the reference positions
    header = open(membersCat, 'rb').readline().split()

    try:
        refInt = np.float(header[1])
    except:
        raise ValueError('Header of the memebers cat is not right')
    
    refRA = np.float(header[2])
    refDEC = np.float(header[3])

             
    
    #1. Log all the clusterMembers and extract them from best.par
    dtypes = [('ID', float), ('deltaRA', float), ('deltaDEC', float), \
                  ('semiMajor', float), ('semiMinor', float), \
                  ('angle', float), ('mag', float), ('z', float) ]

    clusterMembers = np.loadtxt(membersCat, dtype=dtypes)
    
    if refInt == 3:
        ra = refRA - \
          clusterMembers['deltaRA']/np.cos(refDEC/180.*np.pi)/3600.
        dec = refDEC + \
          clusterMembers['deltaDEC']/3600.
    else:
        ra = clusterMembers['deltaRA']
        dec = clusterMembers['deltaDEC']

    clusterMembers = append_rec(clusterMembers, 'RA', ra, \
                                    usemask=False, asrecarray=True)
    clusterMembers = append_rec(clusterMembers, 'DEC', dec, \
                                    usemask=False, asrecarray=True)

    fits.writeto( outFits, clusterMembers, clobber=True )
    

def combineMembersAndBest( catAFileName, catBFileName, outCat):
    
    catA = fits.open(catAFileName)[1].data
    catB = fits.open(catBFileName)[1].data
    GalID = np.arange(len(catB['RA']))
    nNonGalaxies = len(catB) - len(catA)
    galaxyIDs = GalID[ nNonGalaxies:]

    GalaxyFlagArray = np.zeros(len(GalID)) - 1
    GalaxyFlagArray[galaxyIDs] = 1
    GalaxyFlag = []
    GalaxyFlag.append(fits.Column(name='GalaxyFlag', \
                                array=GalaxyFlagArray, \
                                       format='D'))
            
                                 
    catB = fits.BinTableHDU.from_columns(catB.columns + fits.ColDefs(GalaxyFlag)).data

    
    
    catB = append_rec(catB, 'ID', GalID, \
                          usemask=False, asrecarray=True)
    columns = []
    
    for iField in catB.dtype.names:
        if iField == 'identity':
            iFormat='A20'
        else:
            iFormat=np.dtype(float)
                           
        iColumn = fits.Column( name =iField, array=catB[iField],\
                                       format=iFormat)
                                  
        columns.append(iColumn)
        print iField
        
    for iField in catA.dtype.names:
        if iField == 'ID':
            continue
        if iField == 'identity':
            iFormat='A20'
        else:
            iFormat=np.dtype(float)
               

        array = np.zeros(len(GalID)) - 1
        array[galaxyIDs] = catA[iField]
        iColumn = fits.Column( name =iField, array=array,\
                                       format=iFormat)
                                  
        columns.append(iColumn)

        
    
    TotalCat = fits.BinTableHDU.from_columns(columns)
    TotalCat = [fits.PrimaryHDU(), TotalCat]

    #append extra potentials in teh best.fits for
    #example external shear
    catBExtensions = fits.open( catBFileName)
    if len(catBExtensions) >2:
        
        for iExt in xrange(len(catBExtensions)-2):
            TotalCat.append(catBExtensions[iExt+2])

    thdulist = fits.HDUList(TotalCat)
    
    thdulist.writeto(  outCat,  clobber=True )
