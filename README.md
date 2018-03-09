# alignment
Alignment code to simulate clusters

This code is designed to try and understand which galaxy members I
want to constrain within a cluster. To do this is follows mulitple
steps

It has to cuts at the moment. The first cut is argued to to be needed.

1. Convert the galaxy member catalogues in to fits files.
    Measure these galaxies using Sherpa and rewrite the catalogue
	with Sherpa.
2. Append to each galaxy in the member catalogue the stellar mass as
    from the DEEPSPACE photometric catalogues.
3. Determine the ML ratio of each galaxy, where the total mass is the 
    fiducial mass from the cluster model.
4. Now cut in ML ratio and ellipticity (1./0.3 and 0.3 resepctively). 
    THis is to change.

If the number of galaxies after this cut is < the ideal number to
constrain then cut using simulations and get the most sensitive images.
To do this we:

1. It takes the fiducial model (already constrained in the standard
   way) and projects back tot the source plane to create some sources.
2. It then splits the lenses in the fiducial model in to stellar and
   dark matter componnets.
3. It then reprojects the soruces back into the image plane to create
 a new fiducial list of images.
4. Once it has new images it then systematically changes the lenses
    within the model and identifies those images which move the most
    relative to the fiducial image list.
5. It then uses a rule to figure out which lenses cause this. At the 
    moment the rule is that if an image that moves > distCut is within
    nCutRadii factors of its cut radii within that image, then i choose
    that lens.


Once I have selected the lenses it then goes to the original file, it
splits all the lenses that i want to constrain into stellar and dark
matter. It creates a par file with limits for each one and then writes
them out. It then also creates a new galaxy catalogue of sherpa values
that no longer contains the lenses that i want to constrain.
This is put in ../clusters/CLUSTER/CONSTRAIN

The main function is selectLenses.py
import selectLenses as sl
sensitivityDict = sl.selectLenses( \
	cluster, constrain, nConstrain=30, retest=False,
	distCut=0.3, nCutRadii=4)

where cluster is the name of cluster (e.g. A2744, A1063)
constrain is the paramter i want to constrain (ellipticite, angle_pos)
nConstrain is the max number of galaxies I am allowed to constrain
retest is a boolean that says should i retest the galaxies that i have
chosen
distCut is the minimum distance an image must moved to be
sensitive
nCutRadii is the number of radii away from a lens that an image can
lie that is assumed to be within its circle of infuence.

