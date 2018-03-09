'''
offset.py offsets the given potential in some way to see if
the potential alters the position of multiple images


position : offset the potential in x y space
ellipticity : offset the ellipticity
pos_angle : offset the position angle

'''
import astro_tools as at
import copy as cp
from lensing import lenstool as lt
import numpy as np
import ipdb as pdb

def position( potential, offset, centre=None, pos_rot=0. ):
    '''
    Separate the potential into a baryon and dark matter
    and then offset the dark matter potential

    centre : if i want to offset towards some common centre
            should be a tuple of two ra, dec coordinates

     pos_rot : if i want to add on something to the angle to
              change the moving direction

    '''



    if centre is not None:
        ra_diff = at.ra_separation( potential['ra']['float'], \
                                        potential['dec']['float'],
                                    centre[0], potential['dec']['float'])

        dec_diff = (potential['dec']['float'] - centre[1])*3600.

        #Get the angle, but then rotate it 180 so that the dark matter drags
        angle = np.arctan2( dec_diff, ra_diff ) + np.pi + pos_rot

        x_offset = offset*np.cos( angle )
        y_offset = offset*np.sin( angle )
    else:
        #If no centre is given just randomly select one

        angle = np.random.rand()*2.*np.pi
        x_offset = offset*np.cos( angle )
        y_offset = offset*np.sin( angle )
        
 

        
    #dark matter
    potential['x_centre']['float']  += x_offset
    potential['y_centre']['float']  += y_offset

    
    
    return potential





def ellipticity( potential, e_offset, centre=None ):
    '''
    Separate the potential into a baryon and dark matter
    and then offset the dark matter potential in ellipticity

    offset : offset in ellipticity ; always make it more elliptical
             percentrage increase in ellipticity
    

    '''

        
    #dark matter
    #np.random.normal doesnt like a variance of 0!
    if centre is None:
        if e_offset > 0.:
            rand_e_offset = np.random.normal(0., e_offset/100.)
        else:
            rand_e_offset = 0
    else:
        rand_e_offset = e_offset/100.


    a = np.sqrt( 1 + potential['ellipticite']['float'] )
    b = np.sqrt( 1 - potential['ellipticite']['float'] )
    
    a *= (1.0+rand_e_offset/np.sqrt(2.))
    b *= (1.0-rand_e_offset/np.sqrt(2.))
  
    potential['ellipticite']['float']  = (a**2 - b**2)/(a**2 +b**2)

    
    
    return potential

def angle_pos( potential, a_offset, centre=None):
    '''
    Separate the potential into a baryon and dark matter
    and then offset the dark matter potential

    the offset is the rms of the offset
    to do:
    centre : if i want to shift the position angle towards
    the centre 
    

    '''
   

    #np.random.normal doesnt like a variance of 0!
    if centre is None:
        if a_offset > 0.:
            rand_a_offset = np.random.normal( 0., a_offset)
        else:
            rand_a_offset = 0.
    else:
        rand_a_offset = a_offset

    
    potential['angle_pos']['float']  += rand_a_offset


    
    return potential


def mass( potential, m_offset, centre=None ):
    '''
    Offset the dark matter mass by some amount, such that the m/l ratio is
    not always the same

    do this by changing the r_cut of the dark matter halo
    '''


    rand_m_offset = 2
    while np.abs(rand_m_offset) > 1.:
        if centre is None:
            if m_offset > 0.:
                rand_m_offset = np.random.normal( 0., m_offset/100.)
            else:
                rand_m_offset = 0.
        else:
            rand_m_offset = m_offset/100.

   
    
    potential['cut_radius']['float']  *= 1+rand_m_offset
    #Conserve mass since m propto v_disp**2*r_cut
    potential['v_disp']['float'] *= np.sqrt(1./(1.+rand_m_offset))

    if potential['cut_radius']['float'] < 0:
        pdb.set_trace()

    return  potential

def all_offset( potential, \
                    p_offset, \
                    e_offset, \
                    a_offset, \
                    m_offset, \
                    centre=None):
    '''
    This script splits the potential and then adds a random
    ellipticity, angle position, positional offset and mass

    '''


    if p_offset > 0:
        rand_x_offset = np.random.normal( 0., p_offset/np.sqrt(2.))
        rand_y_offset = np.random.normal( 0., p_offset/np.sqrt(2.))
    else:
        rand_x_offset = 0.
        rand_y_offset = 0.

    if e_offset > 0:
        rand_e_offset = np.random.normal( 0., e_offset)
    else:
        rand_e_offset = 0.

    if a_offset > 0:
        rand_a_offset = np.random.normal( 0., a_offset)
    else:
        rand_a_offset = 0.

    if m_offset > 0:
        rand_m_offset = np.random.normal( 0., m_offset)
    else:
        rand_m_offset = 0.

    potential['x_centre']['float']  += rand_x_offset
    potential['y_centre']['float']  += rand_y_offset
    potential['ellipticite']['float']  += rand_e_offset
    potential['angle_pos']['float']  += rand_a_offset
    potential['cut_radius']['float'] *= 1+m_offset
    
    return potential


def split_potential( potential ):
    '''

    Split up the potential into a baryonic and dark matter part

    keep the total mass constant by having the cut radius
    for the baryons the same as the core radius for the dark mattter

    '''
    
    baryons = lt.potentiels.piemd()
    dark_matter = lt.potentiels.piemd()

    for ipar in baryons:
        baryons[ipar]= cp.deepcopy(potential[ipar])
        dark_matter[ipar] = cp.deepcopy(potential[ipar])

    #What shall i cut the baryons cut radius by?
    #originally it was 4, will try 10 now
    ML_ratio = 4.
    #Baryons 
    baryons['cut_radius']['float'] /=  ML_ratio
    baryons['identity'] = cp.deepcopy(potential['identity'])
    baryons['identity']['str'] = potential['identity']['str']+'_baryons'
        
    #dark matter
    
    dark_matter['core_radius']['float'] = baryons['cut_radius']['float']

    dark_matter['identity'] = cp.deepcopy(potential['identity'])
    dark_matter['identity']['str'] = potential['identity']['str']+'_dm'
    
    return baryons, dark_matter
