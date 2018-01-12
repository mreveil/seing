=======================
Developer Information
=======================


To contribute the development of SEING, assuming you have cloned the Github repository,
and created a new branch, there are three (3) files to work with:

* **inputs (.cpp, .h)** : This is where all inputs are read and processed. You will want to add
  any input specific to your fingerprint scheme (such as a name for your fingerprint) to this 
  file.

* **genericlocalcalculator (.cpp, .h)** or **genericglobalcalculator (.cpp, .h)**: This class acts as 
  a way to switch between different fingerprints based on user inputs. Open those files
  and add in your own fingerprints as a new possibility. 

* **your_own_fingerprint (.cpp, .h)** files that implement codes specific to the new fingerprint.
  An existing fingerprint can be used as a starting point for code structure but really, only
  the *calculate_fingerprint* and *get_size* functions are required to be implemented for a 
  fingerprint to be valid. 

After implementation and validatin, please submit a pull request for addition to the master branch.