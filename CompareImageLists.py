
from lensing import lenstool as lt
import numpy as np
import astro_tools_extra as at


def CompareImageLists( catAFileName, catBFileName ):
    '''
    Compare two image lists, and get how much the each image
    has moved from one image to the next
    
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

        
        
