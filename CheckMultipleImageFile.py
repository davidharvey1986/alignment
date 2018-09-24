'''
I need to check the multiple image file to ensure that there aremt
any 'zeros' in the mult file that correspond to a free parameter.

So i need to take the bestopt.par limits on the best constrained redshift
and put it in to a new multiple image file
'''
from lensing import lenstool as lt
import numpy as np

def CheckMultipleImageFile( cluster, OutputMultImageFile=None ):
    '''
    read the bestopt.par and the multiple image 
    file and correct all zeros for the best constrained 
    redshift
    '''
    if OutputMultImageFile is None:
        OutputMultImageFile=cluster+'_mult.img'

    baseDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

    runmode, potentiels, limits, image, source, potfile, grille, champ = \
      lt.read_par( 'bestopt.par', \
                       par_dir = baseDir,\
                       return_all=True)

    dtypes = [('ID', '|S20'), ('RA', float), ('DEC', float), \
                  ('SemiMajor', float), ('SemiMinor', float), \
                  ('Angle', float), ('Redshift', float), ('Mag', float)]


    MultipleImages = np.loadtxt(baseDir+'/'+image['multfile']['filename'],\
                                    dtype=dtypes)
    AllImageFamilies = \
      np.array([i.split('.')[0] \
                        for i in MultipleImages['ID']])

    for iMultiple in image.keys():
        if 'z_m_limit' in iMultiple:

            iImageFamily = np.int(image[iMultiple]['im_label'])
            matchImageFamilies = str(iImageFamily) == AllImageFamilies
            #Assuming that the Gauassian prior is 'lo'
            MultipleImages['Redshift'][ matchImageFamilies ] =\
              image[iMultiple]['lo']

    header = open(baseDir+'/'+image['multfile']['filename']).readline()
    NewMultImageFile = \
      open(baseDir+'/'+OutputMultImageFile, 'wb')
    NewMultImageFile.write(header)
    
    for iImage in MultipleImages:
        writeTuple = \
          tuple([ iImage[i] for i in MultipleImages.dtype.names])
        fmt = '%s %0.10f %0.10f %0.4f %0.4f %0.4f %0.4f %0.4f\n'
        
        NewMultImageFile.write(  fmt % writeTuple)
        
    NewMultImageFile.close()

