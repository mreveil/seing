#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <ctime>

#include "inputs.h"

using namespace std;

fingerprintProperties read_prop_file(string filename){

    string line, key, temp, value;
    ifstream myfile(filename.c_str());

    fingerprintProperties fpproperties;
  //  cout<<"Opening input file: "<<filename<<"\n";
    if (!myfile.is_open()) {
        throw "Cannot open input file";
        return fpproperties;
    }

    double dvalue;
    bool isnetas=false, isnetas2=false, isnzetas=false, isngammas = false;
    bool isnatomtypes=false, istypeprovided=false, isnmaxprovided=false;
    bool isjmaxprovided = false, error = false, isnderivativesprovided=false;
    bool isdirectionprovided=false, isderivativeprovided=false, isinnercutoffprovided=false;
    bool isstrategyprovided=false, isweighttypeprovided=false,isboxsizeprovided=false;
    bool ismodeprovided=false, isoutputfileprovided=false;
    int i = 0;
    double *zetas, *gammas, *etas, *etas2;

    while(getline(myfile,line)) {
        if ( (!line.empty()) and (line.at(0) != '#') ) {

            istringstream buffer(line);
            buffer >> key >> temp >> value;
            stringstream svalue(value);
            
            if (key== "cutoff") {
                 if (temp == "=") 
                     svalue >> fpproperties.cutoff;                                  
            }
            else if (key == "inner_cutoff") {
                 if (temp == "=") {
                      svalue >> fpproperties.inner_cutoff;
                      isinnercutoffprovided = true;
                 }
            } 
            else if (key == "box_size") {
                 if (temp == "=") {
                      isboxsizeprovided = true;
                 
                      fpproperties.box_size = new double[6];
                      svalue >> fpproperties.box_size[0];
                      i = 1;
                      while (buffer >> fpproperties.box_size[i]) {
                         i++;
                      }
                      if (i != 5) cerr<<"ERROR: Box size has to be specified with six values.\n";
                 }
            } 
            else if (key == "type") {
                 if (temp == "=") {
                      svalue >> fpproperties.type;
                      istypeprovided = true;
                 }
            } 
            else if (key == "strategy") { // strategy will weighted or augmented
                 if (temp == "=") {
                      svalue >> fpproperties.strategy;
                      isstrategyprovided = true;
                 }
            } 
            else if (key == "weight_type") { // weight_type can be atomic numbers or electronegativity
                 if (temp == "=") {
                      svalue >> fpproperties.weight_type;
                      isweighttypeprovided = true;
                 }
            } 
            else if (key == "direction") {
                 if (temp == "=") {
                      svalue >> fpproperties.direction;
                      isdirectionprovided = true;
                 }
            } 
            else if (key == "calculate_derivatives") {
                 if (temp == "=") {
                      svalue >> fpproperties.calculate_derivatives;
                      isderivativeprovided = true;
                 }
            } 
            else if (key == "output_file") {
                 if (temp == "=") {
                      svalue >> fpproperties.output_file;
                      isoutputfileprovided = true;
                 }
            } 
            else if (key == "mode") {
                 if (temp == "=") {
                      svalue >> fpproperties.mode;
                      ismodeprovided = true;
                 }
            } 
            else if (key == "nderivatives") {
                 if (temp == "=") {
                      svalue >> fpproperties.nderivatives;
                      isnderivativesprovided = true;
                 }
            } 
            else if (key == "nmax") {
                 if (temp == "=") {
                      svalue >> fpproperties.nmax;
                      isnmaxprovided = true;
                 }
            } 

            else if (key == "jmax") {
                 if (temp == "=") {
                      svalue >> fpproperties.jmax;
                      isjmaxprovided = true;
                 }
            } 


            else if (key == "natomtypes") {
                 if (temp == "=") {
                      svalue >> fpproperties.natomtypes;
                      isnatomtypes = true;                                  
                 }
            }                                      
            else if (key =="atomtypes") {
                if (temp == "=") {
                    if (isnatomtypes == false) {
                         error=true;cerr<<"ERROR: natomtypes has to be specified before atomtypes;";
                    } 
                    else { 
                        fpproperties.atomtypes= new string[fpproperties.natomtypes];
                           svalue >> fpproperties.atomtypes[0];
                           i = 1;
                           while (buffer >> fpproperties.atomtypes[i]) {
                             i++;
                        }
                        if (i != fpproperties.natomtypes)
                           cerr<<"ERROR: You have specified the wrong number of atomtypes\n";
                    }
                }
             }

            else if (key == "nzetas") {
                 if (temp == "=") {
                      svalue >> fpproperties.nzetas;
                      isnzetas = true;                                  
                 }
            }                                      
            else if (key =="zetas") {
                if (temp == "=") {
                    if (isnzetas == false) {
                         error=true;cerr<<"ERROR: nzetas has to be specified before zetas;";
                    } 
                    else { 
                        fpproperties.zetas = new double[fpproperties.nzetas];
                           svalue >> fpproperties.zetas[0];
                           i = 1;
                           while (buffer >> fpproperties.zetas[i]) {
                             i++;
                        }
                        if (i != fpproperties.nzetas)
                           cerr<<"ERROR: You have specified the wrong number of zeta values\n";
                    }
                }
             }
             else if (key =="ngammas"){
                 if (temp == "=") {
                     svalue >> fpproperties.ngammas;
                     isngammas = true;
                 }
             }
             else if (key =="gammas"){
                   if (temp == "=") {
                       if (isngammas == false) {
                           error=true;cerr<<"ERROR: ngammas has to be specified before gammas;";
                       } 
                       else{
                            fpproperties.gammas = new double[fpproperties.ngammas];
                            svalue >> fpproperties.gammas[0];
                            i = 1;
                            while (buffer >> fpproperties.gammas[i]) {
                                i++;
                            }
                            if (i != fpproperties.ngammas)
                                cerr<<"ERROR: You have specified the wrong number of gamma values\n";
                       }
                  }
              }
              else if (key =="netas"){
                  if (temp == "=") {
                      svalue >> fpproperties.netas;
                      isnetas = true;
                  }
              }
              else if (key == "etas"){
                  if (temp == "=") {
                      if (isnetas == false) {
                          error=true;cerr<<"ERROR: netas has to be specified before etas;";
                      } 
                      else{
                           fpproperties.etas = new double[fpproperties.netas];
                           svalue >> fpproperties.etas[0];
                           i = 1;
                           while (buffer >> fpproperties.etas[i]) {
                                i++;
                           }
                           if (i != fpproperties.netas)
                               cerr<<"ERROR: You have specified the wrong number of eta values\n";

                      }
                  }
              }
              else if (key =="netas2"){
                  if (temp == "=") {
                      fpproperties.netas2 = atoi(value.c_str());
                      isnetas2 = true;
                  }
              }
              else if (key =="etas2"){
                   if (temp == "=") {
                       if (isnetas2 == false) {
                            error=true;cerr<<"ERROR: netas2 has to be specified before netas2;"; 
                       }
                       else {
                           fpproperties.etas2 = new double[fpproperties.netas2];
                           svalue >> fpproperties.etas2[0];
                           i = 1;
                           while (buffer >> fpproperties.etas2[i]) {
                                i++;
                          }
                          if (i != fpproperties.netas2)
                           cerr<<"ERROR: You have specified the wrong number of etas2 values\n";

                      }
                   }
              }             
        }
    }


/////////////// Taking care of defaults //////////////////////////

    if (istypeprovided == false) {
        fpproperties.type = "gaussian";
        cout<<"No fingerprint type provided in input file. The default (gaussian) will be used.\n";
    }
    if (fpproperties.type == "zernike" and isnmaxprovided==false){
        fpproperties.nmax = 10;
        cout<<"No 'nmax' value provided for zernike fingerprint in input file. The default (10) will be used.\n";

    }
    if (fpproperties.type == "bispectrum" and isjmaxprovided==false){
        fpproperties.jmax = 10; 
        cout<<"No 'jmax' value provided for bispectrum fingerprint in input file. The default (10) will be used.\n";

    }
    if (fpproperties.type == "diamond" and isinnercutoffprovided==false){
        fpproperties.inner_cutoff = 1.0; 
        cout<<"No 'inner_cutoff' provided for diamond fingerprint in input file. The default (1.0) will be used.\n";

    }
    if (isderivativeprovided == false) {
        fpproperties.calculate_derivatives = "false";
        cout<<"Derivatives of the fingerprints will not be calculated.\n";
    }
    if (isderivativeprovided == true and isdirectionprovided == false) {
        fpproperties.direction = 0;
        cout<<"No direction provided. Derivatives will be calculated in the 'x' direction.\n";
    }
    if (isnderivativesprovided == false) {
        fpproperties.nderivatives = 1;
        if (isderivativeprovided == true)
            cout<<"No number of derivatives provided. Default of 1 will be used.\n";
    }
    if (isstrategyprovided == false) {
        fpproperties.strategy = "";
        cout<<"No strategy provided to treat multiple atom types. Default associated with type of fingerprint will be used.\n";
    }
    if (isweighttypeprovided == false) {
        fpproperties.weight_type = "atomic_number";
        if (fpproperties.strategy == "weighted")
            cout<<"No weight type provided. Default (atomic number) will be used.";
    }
    if (isboxsizeprovided == false) {
        fpproperties.is_box_size_provided = false;
    }
    else {
        fpproperties.is_box_size_provided = true;
    }
    if (ismodeprovided == false) {
        fpproperties.mode = "append";
        cout<<"Note: Fingeprints will be appended to output file";
    }
    if (isoutputfileprovided = false) {
        fpproperties.output_file = fpproperties.type+"_fingerprints.sg";
        cout<<"No output file name provided, default will be used.";

    }
///////////////  Some Input Validation /////////////////////////////// TODO: validate value given for mode, box size, etc...

    if (isderivativeprovided) {

        if (fpproperties.nderivatives < 1 or fpproperties.nderivatives > 100) {
            cout<<"INPUT ERROR: Please provide a valid number of derivatives (>= 1). Current value: "<<fpproperties.nderivatives<<"\n";
            error = true;
        }

        if (fpproperties.calculate_derivatives != "false" and fpproperties.calculate_derivatives != "true") {
            cout<<"INPUT ERROR: Please provide a valid value for calculate_derivatives (true or false)\n";
            error = true;
        }
    }
    if (error) throw "bad input";


///////////////// Finishing up //////////////////////

    cout<<"Done reading fingerprint file. This is what we have:\n";
    cout<<"\tFingerprint type is: "<<fpproperties.type<<"\n";
    cout<<"\tNumber of atom types: "<<fpproperties.natomtypes<<" (";
    for (int i=0; i<fpproperties.natomtypes;i++) cout<<" "<<fpproperties.atomtypes[i];
    cout<<"\tCalculate fingerprint derivatives: "<<fpproperties.calculate_derivatives<<"\n";
    if (fpproperties.calculate_derivatives == "true") cout<<"\t\tNumber of derivatives to include: "<<fpproperties.nderivatives<<"\n";
    cout<<")\nDone with fingerprint properties\n";
    return  fpproperties;
        
}
