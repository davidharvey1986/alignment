
import astro_tools_extra as at
import numpy as np
import pyfits as pyfits
import ipdb as pdb

def getSherpaShapes( MemberCatalogue, ClusterImage,
                         outFileName, \
                         pixelScale=0.03):
    '''
    This script will take a cluster member catalogue
    find those galaxies in the ClusterImage, measure their
    shapes using SHerpa and rewrite the catalogue.
    '''

    image = pyfits.open( ClusterImage )[0].data
    
    dtypes = [('ID', '|S20'), ('deltaRA', float), ('deltaDEC', float), \
                  ('SemiMajor', float), ('SemiMinor', float), \
                  ('Angle', float), ('Mag', float), ('Redshift', float)]

    GalCat = np.loadtxt( MemberCatalogue, dtype=dtypes)

    
    CatHeader = open(MemberCatalogue, 'rb').readline()

    refRA = np.float(CatHeader.split()[2])
    refDEC = np.float(CatHeader.split()[3])
    
    outFile = open(outFileName, 'wb')
    outFile.write(CatHeader)
    
    for iGal in GalCat:
        iGalDict = Galaxy(iGal, ClusterImage, image, \
                        refRA, refDEC, pixelScale=pixelScale )
        
        iGalDict.getPostageSize()
        iGalDict.getPostageStamp(image )
        
        iGalDict.getEllipse()
      
        iIter = 0
        while (iGalDict['SherpaMajor'] > iGalDict['postageSize']):
            #If the estiamted major axis is greater than the
            #postage stamp then its wrong, make the postage stamp
            #bigger and redo
            iIter += 1
            iGalDict['postageSize'] *= 2.
            iGalDict.getPostageStamp(image )
            iGalDict.getEllipse()
                    
                

   
        TupleWrite = ( iGalDict['ID'], \
                           iGalDict['NewDeltaRA'], \
                           iGalDict['NewDeltaDEC'], \
                           iGalDict['SherpaMajor'], \
                           iGalDict['SherpaMinor'], \
                           iGalDict['SherpaAngle'], \
                           iGalDict['Mag'], \
                           iGalDict['Redshift'])
                
        outFile.write('%s %0.7f %0.7f %0.3f %0.3f %f %0.5f %0.1f\n' % TupleWrite)

    outFile.close()



          



    
class Galaxy(dict):

    def __init__( self, iGal, ClusterImage, Image, \
                      refRA, refDEC, pixelScale=0.03 ):

        
        
        for iName in iGal.dtype.names:
            self[iName] = iGal[iName]
            
        self['refRA'] = refRA
        self['refDEC'] = refDEC
        self['ClusterImage'] = ClusterImage
        self['pixelScale'] = pixelScale

        self['RA'] = -1.*self['deltaRA']/\
          (3600.*np.cos(self['refDEC']/180.*np.pi)) + self['refRA']
        self['DEC'] = self['deltaDEC']/3600. + self['refDEC']

    def __getattr__(self, attr):
        return self[attr]
    
    def getPostageSize( self ):
        '''
        Assuming the postage stamp size is the 2.* the semi major axis
        '''

        self['SexEll'] = \
          np.sqrt( (self['SemiMajor']**2 - self['SemiMinor']**2) / \
                    (self['SemiMajor']**2 + self['SemiMinor']**2))
                   
        self['postageSize'] = self['SemiMajor'] / \
          np.sqrt( 1. + self['SexEll'])/self['pixelScale']

          
    def getPostageStamp( self, image ):
        
        self['xImage'], self['yImage'] \
          = at.deg2pix( self['ClusterImage'], \
                            [self['RA']], [self['DEC']])

        xBegin = np.int(self['xImage']-self['postageSize'])
        yBegin = np.int(self['yImage']-self['postageSize'])
        xEnd = np.int(self['xImage']+self['postageSize'])
        yEnd = np.int(self['yImage']+self['postageSize'])

                        
        self['postageStamp'] = \
          image[ yBegin:yEnd, xBegin:xEnd ]
          
    def getEllipse( self ):
        '''
        Use Sherpa to get the ellipticity of the isophote
        '''
        
 
        self.sherpaIsophote( )

                                            
        self['SherpaMajor'] = \
          self['SherpaSize']*np.sqrt(1. + self['SherpaEll'])*self['pixelScale']
        self['SherpaMinor'] = \
          self['SherpaSize']*np.sqrt(1. - self['SherpaEll'])*self['pixelScale']




    def sherpaIsophote( self):
        '''
        use ciao to run the measurement
        '''
        print 'sherpa on ',self['ID']
        from sherpa.astro.ui import *
        
        pyfits.writeto('postage.fits', self['postageStamp'], \
                           clobber=True)
    
        load_image("postage.fits")
        

        set_source(sersic2d.g2)
        guess(g2)
        
        g2.n.thaw()
        g2.xpos.max =  self['postageStamp'].shape[1]
        g2.ypos.max = self['postageStamp'].shape[0]
        g2.xpos = self['postageSize']
        g2.ypos = self['postageSize']

        g2.ellip = self['SexEll']
        g2.theta = self['Angle']/180.*np.pi
        
        fit()
        
            
        
        NewX = self['xImage'] - self['postageSize'] + g2.xpos.val
        NewY = self['yImage'] - self['postageSize'] +  g2.ypos.val
        self['SherpaEll'] = g2.ellip.val
        self['SherpaAngle'] = g2.theta.val / np.pi *180.
        
        self['SherpaRA'], self['SherpaDEC'] =\
          at.pix2deg( self['ClusterImage'], \
                          [NewX], [NewY] )
        self['SherpaSize'] = g2.r0.val
        self['NewDeltaRA'] = at.ra_separation( self['SherpaRA'],
                                            self['refDEC'], \
                                            self['refRA'], \
                                            self['refDEC'] )
                                            
        self['NewDeltaDEC'] = (self['SherpaDEC'] - self['refDEC'])*3600.
