SEING v 0.1
============


SEING is a C/C++ package for fingerprint calculations suitable for machine learning studies 
of molecular systems. 

SEING was developed in the Clancy Group (http://clancygroup.cbe.cornell.edu/) at Cornell University.

Author: Mardochee Reveil (mr937@cornell.edu)

Fingerprints (in this context) are numerical representations of chemical environments designed 
to be invariant under property-perseving operations such as permutation of atoms of the same 
nature, geometric rotation, etc. For more information on fingerprints in general and the ones 
currently implemented in SEING, please see the official documentation and user-guide.

"SEING" is an old French word for signature.

DOCUMENTATION
--------------

The official documentation and user-guide can be found here: https://seing.readthedocs.io


INSTALLATION
-------------

SEING is built with minimal requirements and can be easily compiled with a suitable
C/C++ compiler. A generic Makefile is provided in src folder. As a starting point,
you can just type.

```
cd seing
mkdir bin			
cd src
make seing
```

If this doesn't work, changes might be necessary to adapt the makefile to your 
operating system and/or environment.

Please note that c++11 or later is required. I have successfully compiled with GCC 4.8.4 and 4.9.2 on CentOS.
If you are using GCC and run into C++ standard related errors, try adding -std=c++11 to your compile command.
I will look into enabling compatibility with earlier C++ standards in future releases.

LICENSE
----------

This program is free and open-source software distributed under the terms of the GNU GPL version 3 
(or later) which can be found here: www.gnu.org/licenses/gpl-3.0.en.html

Please note that SEING is provided WITHOUT WARRANTY OF ANY KIND, either expressed or implied, including,
but not limited to, the implied warranties of merchantability and fitness for a particular purpose. 
Please see the full terms of the GNU GPL license for more details. 


CONTRIBUTIONS
--------------

We welcome contributions to this project including implementation of new fingerprinting 
schemes, bug tracking and corrections, code optimization, documentation, etc. Please consult the 
"developer" section of the documentation for more information on how the code is organized. To make a contribution,
create your own branch, make your documented changes to the code and submit a pull request for code update.


USER SUPPORT
-------------

SEING is provided with no dedicated user support, however questions and suggestions are welcome and the author(s)
will do their best to provide answers in a timely fashion.


CITATION
----------

If you use this software in your research, please cite the appropriate paper(s) for your chosen fingerprinting method(s)
as well as the official SEING paper published here: http://pubs.rsc.org/en/content/articlelanding/2018/me/c8me00003d#!divAbstract

Mardochee Reveil and Paulette Clancy, "Classification of spatially resolved molecular fingerprints for machine learning applications and development of a codebase for their implementation." Molecular Systems Design & Engineering, published online 2/20/18. Invited paper, part of a themed issue on ‘Machine Learning and Data Science in Materials Design’ DOI: 10.1039/C8ME00003D
