
import CreateNewClusterModel as newModel
import CompareImageLists as CIL
import numpy as np
import pickle as pkl
import os as os
import ipdb as pdb

def getImageSensitivity( cluster, sourceList, FiducialImageList,
                             constrain, rerun=False, 
                             nRuns=20, potentials=None,
                             pklName=None,
                             p_offset=0.1, \
                             e_offset=20, \
                             a_offset=20., \
                             m_offset=0.):
    '''
    Determine how sensitive each multiple image is to snall
    pertubations in the cluster model.

    To do this I will take the fiducial cluster model
    stretch, rotate or move the galaxies withint he cluster
    and see how the multiple images move

    clusrter : name of the cluster
    sourceList : fiducial source list
    fiducialImageList : the fiducial image list to compare to

    '''
    lenstoolDir = '/Users/DavidHarvey/Documents/Work/alignment/'+\
         'clusters/'+cluster

    if pklName is None:
        pklName = lenstoolDir+'/'+cluster+'_'+constrain+\
          '_imageSensitivity.pkl'
    else:
        pklName = lenstoolDir+'/'+pklName
    
    if (os.path.isfile( pklName)) & (not rerun):
        with open(pklName,'rb') as f:
            return pkl.load(f)
    
    
         
    NewImageList = 'image_tmp.dat'
    #Create a new model
    TotalImageDifference = {}
    
    for iRun in xrange(nRuns):
        newModel.CreateNewClusterModel( cluster, sourceList, constrain, \
                                            imageList=NewImageList,
                                            p_offset=p_offset, \
                                            e_offset=e_offset, \
                                            a_offset=a_offset, \
                                            m_offset=m_offset, \
                                            potentials=potentials)

        #4. Compare the image lists
        iImageDifference = \
          CIL.CompareImageLists( lenstoolDir+'/'+FiducialImageList, 
                                     lenstoolDir+'/'+NewImageList  )

        if iRun ==0:
            for iImage in iImageDifference:
                TotalImageDifference[iImage['ID']] = {}
                TotalImageDifference[iImage['ID']]['Delta'] = np.zeros(nRuns)

        
        for iImage in iImageDifference:
            TotalImageDifference[iImage['ID']]['Delta'][iRun] \
              = iImage['Delta']
            TotalImageDifference[iImage['ID']]['RA'] = iImage['RA']
            TotalImageDifference[iImage['ID']]['DEC'] = iImage['DEC']
            TotalImageDifference[iImage['ID']]['Redshift'] = \
              iImage['Redshift']

        
    pkl.dump(TotalImageDifference, open(pklName,'wb'))
        
    return TotalImageDifference
        

