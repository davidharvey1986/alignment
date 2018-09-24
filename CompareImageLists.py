
from lensing import lenstool as lt
import numpy as np
import astro_tools_extra as at
import ipdb as pdb

def CompareImageLists( catAFileName, catBFileName, noFind=1.0 ):
    '''
    Compare two image lists, and get how much the each image
    has moved from one image to the next

    noFind = when the image disappears and cant be compared
    the defaul value
    
    '''
    dtypes = [('ID', '|S20'), ('RA', float), ('DEC', float), \
                  ('SemiMajor', float), ('SemiMinor', float), \
                  ('Angle', float), ('Redshift', float), ('Mag', float)]
    
    catA = np.loadtxt( catAFileName, dtype=dtypes)
    catB = np.loadtxt( catBFileName, dtype=dtypes)

    #Create a catalogue of ID labels with positional movement
    
    outCat = np.array([], dtype=[('ID', '|S20'), \
                                     ('Delta', float), \
                                     ('RA', float), \
                                     ('DEC', float), \
                                     ('Redshift', float)])
    for iImage in catA:
        IDindex =  catB['ID'] == iImage['ID']
        
        if np.all(IDindex == False):
      
            #This happens if the new model doesnt produce the image
            newEntry = np.array((iImage['ID'], np.min(noFind), \
                                iImage['RA'], iImage['DEC'], \
                                iImage['Redshift']), \
                             dtype=outCat.dtype)
            outCat = np.append( outCat, newEntry)
            continue
        
        NewCatRA = catB['RA'][ IDindex ]
        NewCatDEC = catB['DEC'][ IDindex ]

        Separation = at.ra_separation( NewCatRA, NewCatDEC, \
                                           iImage['RA'], iImage['DEC'], 
                                           abs=True)

        newEntry = np.array((iImage['ID'], np.min(Separation), \
                                iImage['RA'], iImage['DEC'], \
                                iImage['Redshift']), \
                             dtype=outCat.dtype)
        outCat = np.append( outCat, newEntry)
        
    return outCat

        
        
