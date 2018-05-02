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


    if (!myfile.is_open()) {
        throw "Cannot open input file";
        return fpproperties;
    }


    fpproperties.type = "";
    fpproperties.natomtypes = 0;
    fpproperties.strategy = "";
    fpproperties.weight_type = "";

    fpproperties.calculate_derivatives = "";
    fpproperties.ndirections = 0;
    fpproperties.nderivatives = 1;

    fpproperties.cutoff = 6.5;
    fpproperties.inner_cutoff = -1;
    fpproperties.is_box_size_provided = false;

    fpproperties.output_file = "";
    fpproperties.output_mode = "";

    fpproperties.nmax = -1;
    fpproperties.jmax = -1;

    fpproperties.nzetas = 0;
    fpproperties.ngammas = 0;
    fpproperties.netas = 0;
    fpproperties.netas2 = 0;   

    fpproperties.width = 0;
    fpproperties.dimensionality = 0;
    fpproperties.alpha = -1;       

    double dvalue;
    bool error = false;

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
                 }
            } 
            else if (key == "box_size") {
                 if (temp == "=") {
                      fpproperties.is_box_size_provided = true;
                 
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
                 }
            } 
            else if (key == "strategy") { // strategy will weighted or augmented
                 if (temp == "=") {
                      svalue >> fpproperties.strategy;
                 }
            } 
            else if (key == "weight_type") { // weight_type can be atomic numbers or electronegativity
                 if (temp == "=") {
                      svalue >> fpproperties.weight_type;
                 }
            } 
            else if (key == "ndirections") {
                 if (temp == "=") {
                      svalue >> fpproperties.ndirections;
                 }
            }
            else if (key == "directions"){
                  if (temp == "=") {
                      if (fpproperties.ndirections <= 0 || fpproperties.ndirections > 3) {
                          error=true;cerr<<"ERROR: Please provide a valid value for ndirections between 0 and 3;";
                      } 
                      else{
                           fpproperties.directions = new int[fpproperties.ndirections];
                           svalue >> fpproperties.directions[0];
                           i = 1;
                           while (buffer >> fpproperties.directions[i]) {
                                i++;
                           }
                           if (i != fpproperties.ndirections)
                               cerr<<"ERROR: You have specified the wrong number of directions for derivative calculations\n";

                      }
                  }
            }
            else if (key == "calculate_derivatives") {
                 if (temp == "=") {
                      svalue >> fpproperties.calculate_derivatives;
                 }
            } 
            else if (key == "output_file") {
                 if (temp == "=") {
                      svalue >> fpproperties.output_file;
                 }
            } 
            else if (key == "output_mode") {
                 if (temp == "=") {
                      svalue >> fpproperties.output_mode;
                 }
            } 
            else if (key == "nderivatives") {
                 if (temp == "=") {
                      svalue >> fpproperties.nderivatives;
                 }
            } 
            else if (key == "nmax") {
                 if (temp == "=") {
                      svalue >> fpproperties.nmax;
                 }
            } 

            else if (key == "jmax") {
                 if (temp == "=") {
                      svalue >> fpproperties.jmax;
                 }
            } 

            else if (key == "width") {
                 if (temp == "=") {
                      svalue >> fpproperties.width;
                 }
            } 

            else if (key == "alpha") {
                 if (temp == "=") {
                      svalue >> fpproperties.alpha;
                 }
            } 

            else if (key == "dimensionality") {
                 if (temp == "=") {
                      svalue >> fpproperties.dimensionality;
                 }
            } 

            else if (key == "natomtypes") {
                 if (temp == "=") {
                      svalue >> fpproperties.natomtypes;
                 }
            }                                      
            else if (key =="atomtypes") {
                if (temp == "=") {
                    if (fpproperties.natomtypes <= 0) {
                         error=true;cerr<<"ERROR: a valid value for 'natomtypes' has to be specified before atomtypes;";
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
                 }
            }                                      
            else if (key =="zetas") {
                if (temp == "=") {
                    if (fpproperties.nzetas <= 0 ) {
                         error=true;cerr<<"ERROR: A valid value for 'nzetas' has to be specified before 'zetas';";
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
                 }
             }
             else if (key =="gammas"){
                   if (temp == "=") {
                       if (fpproperties.ngammas <= 0) {
                           error=true;cerr<<"ERROR: A valid value for 'ngammas' has to be specified before 'gammas';";
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
                  }
              }
              else if (key == "etas"){
                  if (temp == "=") {
                      if (fpproperties.netas <= 0) {
                          error=true;cerr<<"ERROR: A valid value for 'netas' has to be specified before 'etas';";
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
                  }
              }
              else if (key =="etas2"){
                   if (temp == "=") {
                       if (fpproperties.netas2 <= 0) {
                            error=true;cerr<<"INPUT ERROR: A valid value for 'netas2' has to be specified before 'netas2';"; 
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


/////////////// Taking care of defaults and bad inputs //////////////////////////

    if (fpproperties.type == "") {
        fpproperties.type = "gaussian";
        cout<<"NOTE: No fingerprint type provided in input file. The default (gaussian) will be used.\n";    
    }

    if (fpproperties.type == "zernike" and fpproperties.nmax == -1){
        fpproperties.nmax = 10;
        cout<<"NOTE: No 'nmax' value provided for zernike fingerprint in input file. The default (10) will be used.\n";
    }

    if (fpproperties.type == "bispectrum" and fpproperties.jmax == -1){
        fpproperties.jmax = 10; 
        cout<<"NOTE: No 'jmax' value provided for bispectrum fingerprint in input file. The default (10) will be used.\n";
    }

    if (fpproperties.type == "diamond" and fpproperties.inner_cutoff == -1){
        fpproperties.inner_cutoff = 1.0; 
        cout<<"NOTE: No 'inner_cutoff' provided for diamond fingerprint in input file. The default (1.0) will be used.\n";
    }


    if (fpproperties.calculate_derivatives == "") {
        fpproperties.calculate_derivatives = "false";
        cout<<"NOTE: You didn't tell me if I need to calculate derivatives or not. By default, I don't.\n";
    }
    else if (fpproperties.calculate_derivatives == "false") {
        cout<<"NOTE: Derivatives of the fingerprints will not be calculated.\n";
    }
    else if (fpproperties.calculate_derivatives == "true") {
        cout<<"NOTE: Derivatives of the fingerprints will be calculated.\n";
        if (fpproperties.nderivatives == -1) {
            fpproperties.nderivatives = 1;
            cout<<"NOTE: No number of derivatives provided. Default of 1 will be used.\n";
        }
        else if (fpproperties.nderivatives > 5) {
            cout<<"NOTE: Too many number of derivatives specified. Max of 4 will be calculated.\n";
            fpproperties.nderivatives = 4;
        }
    }
    else {
         error = true;
         cout<< "INPUT ERROR: Invalid value for calculate_derivatives: "<<fpproperties.calculate_derivatives<<". It can only be 'true' or 'false';\n";
    }



    if (fpproperties.strategy == "") {
        cout<<"NOTE: No strategy provided to treat multiple atom types. Default associated with type of fingerprint will be used.\n";
    }
    else if (fpproperties.strategy != "augmented" and fpproperties.strategy != "weighted") {
         error = true;
         cout<<"INPUT ERROR: Invalid value for strategy. Only be 'augmented' or 'weighted' are allowed\n";
    }


    if (fpproperties.weight_type == "") {
        fpproperties.weight_type = "atomic_number";
        if (fpproperties.strategy == "weighted")
            cout<<"NOTE: No weight type provided. Default (atomic number) will be used.\n";
    }
    else if (fpproperties.weight_type != "atomic_number" and fpproperties.weight_type != "electronegativity" and fpproperties.weight_type != "constant") {
         error = true;
         cout<<"INPUT ERROR: Unsupported weight_type. Can only be 'atomic_number' or 'electronegativity' or 'constant'.\n";
    }


    if (fpproperties.output_mode == "") {
        fpproperties.output_mode = "append";
        cout<<"NOTE: No output_mode specified. Fingeprints will be appended to output file\n";
    }
    else if (fpproperties.output_mode != "append" and fpproperties.output_mode != "overwrite") {
         error = true;
         cout<<"INPUT ERROR: Unsupported output_mode. Can only be 'append' or 'overwrite'.\n";
    }


    if (fpproperties.output_file == "") {
        fpproperties.output_file = fpproperties.type+"_fingerprints.sg";
        cout<<"No output file name provided, default will be used.";

    }

    if (error) throw "FATAL: Bad input";

///////////////// Finishing up //////////////////////

    cout<<"Done reading fingerprint file. This is what we have:\n";
    cout<<"\tFingerprint type is: "<<fpproperties.type<<"\n";
    cout<<"\tNumber of atom types: "<<fpproperties.natomtypes<<" (";
    for (int i=0; i<fpproperties.natomtypes;i++) cout<<" "<<fpproperties.atomtypes[i];
    cout<<")\n";
    if (fpproperties.natomtypes > 1) {
        cout<<"\tStrategy to handle multiple species: "<<fpproperties.strategy;
        if (fpproperties.strategy == "weighted") cout<<" (Weight type:"<<fpproperties.weight_type<<")";
        cout<<"\n";
    }
    cout<<"\tCalculate fingerprint derivatives: "<<fpproperties.calculate_derivatives<<"\n";
    if (fpproperties.calculate_derivatives == "true") {
        cout<<"\tNumber of derivatives to include: "<<fpproperties.nderivatives<<"\n";
        cout<<"\tNumber of directions to include in derivatives: "<<fpproperties.ndirections<<"\n";
    }
   // cout<<")\nDone with fingerprint properties\n";
    return  fpproperties;
        
}
