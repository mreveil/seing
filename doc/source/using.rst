HOW TO USE
============



To use SEING, a coordinate file and an option file are required. 
Example usage looks like this::

   /path/to/seing coordinate.xyz optionfile.in


TRAJECTORY FILE
----------------

Coordinates of each atom has to be provided in a coordinate file in the xyz format. 
Only the xyz file format for atomic coordinates is supported at the moment.
Trajectory files (i-e coordinate files with more than one frame) are also not supported at the moment
but support will be added soon.



OPTION FILE 
---------------

The option input file contains "key = value" pairs specifying the type of calculation to
perform, the input parameters for the method chosen, etc. Current keys and possible values are as follow:


General Options
******************
----------------------------

**type** (*optional*)
   The type of fingerprinting scheme to use.
   Possible values are: 

      * gaussian (*default*):
      * zernike:
      * bispectrum: 

**natomtypes** (*required*)
   The number of different species in the molecular system.
      * Integer values only


**atomtypes** (*required*)
   A space separated list of the abbreviated names of the different species in the system. Names should correspond
   to the ones used in the coordinate file.
      * Number of species provided has to match the value for the *natomtypes* option above

**strategy** (*optional*)
   How to account for more than one species.
   Possible values are:
      * augmented (*default*): the fingerprint size is increased with one subfingerprint for each different atom pair or triplets 
      * weighted: the size of the fingerprint remains the same (as in with just one atom type) but the contribution of each atom 
                  type is weighted based on a specified weight type. Please note: this strategy doesn't work with all fingerprints.

**weight_type** (*optional*)
   Defines how contributions are weighted for the *weighted* strategy explained above. Possible values are:
      * atomic number (*default*)
      * electronegativity

Derivatives Options
*******************
-------------------------------------

**calculate_derivatives** (*optional*)
   Whether or not to calculate fingerprint derivatives. If yes, derivatives are added to the fingerprint. Please see documentation of 
   your specific fingerprint for whether or not derivatives are supported and if so, how they are calculated and incorporated to the
   fingerprint vector or matrix.
      * yes
      * no (*default*)

**direction** (*optional*)
   The direction (x, y or z) for derivatives of the fingerprints to be calculated
      * 0 (*default*): Calculate only in the x direction
      * 1: Calculate only in the y direction
      * 2: Calculate only in the z direction
      * 3: Calculate in all three directions 

**nderivatives** (*optional*)
   The number of derivatives to calculate. *Default is one and is with respect to the center atom.* If a value greater than 1 provided,
   derivatives are calculated with respect to.


Output Options
****************
--------------------------------------

**output_file** (*optional*)
   Name of the output file to write the fingerprint in. Output file will be in current directory (where the coordinate and option file are).
   If the file already exists, the behavior of the program is determined by the *mode* keyword explained below. 
   * Default output name is *fingerprint_type+"_fingerprints.sg"*
    

**mode** (*optional*)
    Whether to append fingerprints to the given output file, if it already exists. If not, file will be overwritten
    * append
    * overwrite (*default*)


Neighbor Searching Options
***************************
--------------------------------


**cutoff** (*required*)
   Defines the cutoff value used to build the neighbor list.


**box_size** (*optional*)
   Defines the size of the simulation box in the following format: *xmin ymin zmin xmax ymax zmax*


Fingerprint-Specific Options
****************************
------------------------------

*Bispectrum*

**jmax**
   

*Behler-Parinello (Gaussian)*

**nmax** 
   
**nzetas** 

**zetas** 

**ngammas** 

**gammas** 

**netas** 

**etas** 

**netas2** 

**etas** 












