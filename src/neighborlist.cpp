// Created by Mardochee Reveil on 9/27/18
// Class to create and manage neighbor list

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "neighborlist.h"

using namespace std;

NeighborList::NeighborList(void):natoms(0), totalbins(0), binlist(NULL),heads(NULL), neighbors(NULL), 
 atomsperbin(NULL), neighboringbins(NULL), neighborsperatom(NULL) {


}


NeighborList::NeighborList(AtomicSystem& asystem, double cutoff, int nxb, int nyb, int nzb, int maxneigh): 
 natoms(0), totalbins(0), binlist(NULL),heads(NULL), neighbors(NULL), atomicsystem(asystem), 
 atomsperbin(NULL), neighboringbins(NULL), neighborsperatom(NULL) {

    natoms = atomicsystem.get_n_atoms();
    maxneighbors = maxneigh;
    nxbins = nxb;
    nybins = nyb;
    nzbins = nzb;
    cutoffsq = cutoff*cutoff;

    totalbins = nxbins*nybins*nzbins;

    if (totalbins == 0) {

        if (nxbins == 0) nxbins = 1;
        if (nybins == 0) nybins = 1;
        if (nzbins == 0) nzbins = 1;
        totalbins = 1;

    }

    neighboringbins = new int*[totalbins];
    for (int i=0; i<totalbins; i++) neighboringbins[i] = new int[26];

    find_neighboring_bins();

    cout << "\n Total bins: " << totalbins << "\n";
    cout << "    Cutoff square is: " << cutoffsq << "\n";

    minxbox = atomicsystem.get_xmin();
    minybox = atomicsystem.get_ymin();
    minzbox = atomicsystem.get_zmin();

    xboxsize = atomicsystem.get_xsize();
    yboxsize = atomicsystem.get_ysize();
    zboxsize = atomicsystem.get_zsize();

    zbinsize = xboxsize/nxbins;
    ybinsize = yboxsize/nybins;
    zbinsize = zboxsize/nzbins;

    if (maxneighbors > natoms)
        maxneighbors = natoms; 

    neighbors = new int*[natoms];
    for (int i=0; i<natoms; i++)    neighbors[i] = new int[maxneighbors];
    neighborsperatom = new int[natoms];

    for (int i=0; i<natoms; i++) {

        neighborsperatom[i] = 0;
        for ( int j=0; j<maxneighbors; j++)
            neighbors[i][j] = -1;

    }

    binlist = new int[natoms];
    atomsperbin = new int[totalbins];
    heads = new int[totalbins];

}

void NeighborList::build() {

    bool success = initialize_binning();
    if (success == false) {
        throw "Failed to complete binning";
    }

    for (int i=0; i<natoms; i++) {
        neighborsperatom[i] = 0;
        for (int j=0; j<maxneighbors; j++)
            neighbors[i][j] = -1;
    }

    for (int atomid=0; atomid<natoms; atomid++) {
        find_neighbors(atomid);
    }

    double totalneighbors = 0;
    for (int i=0; i<natoms; i++){
        totalneighbors += neighborsperatom[i];
     /*   cout<<i<<" ("<<neighborsperatom[i]<<") ";
        for (int j=0; j<neighborsperatom[i]; j++)
            cout<<neighbors[i][j]<<" ";
        cout<<"\n";
        if (neighborsperatom[i] == 0) cout<<"Problem\n";*/
    }
    cout << "    Total neighbors for all atoms is: " << totalneighbors<< "\n";

}

NeighborList::~NeighborList() {

    delete [] neighbors;
    delete [] atomsperbin;
    delete [] neighborsperatom;
    delete [] binlist;
    delete [] heads;
    delete [] neighboringbins;
 
}

int* NeighborList::get_atoms_in_bin(int binnumber) {

    int * atomsinbin = new int[atomsperbin[binnumber]];
    int i = 0, atomid;

    if (heads[binnumber] == -1) cout << " WARNING: found an empty bin";

    else { 
        atomid = heads[binnumber];
        atomsinbin[i] = atomid;
        while(binlist[atomid] != -1) {
            i++;
            if ( i >= atomsperbin[binnumber]) cout<<"ERROR: there are more atoms in this bin than saved\n";
            atomid = binlist[atomid];
            atomsinbin[i] = atomid;
        }

    } 

    return atomsinbin;
}

int NeighborList::get_bin_number(double x, double y, double z) {

    int bin=0, binx=0, biny =0, binz=0, ibin;
    
    for (ibin=0; ibin<nxbins; ibin++)
        if ((x>=minxbox+ibin*xbinsize) && (x<minxbox + (ibin+1)*xbinsize)) {binx=ibin; break;}

    for (ibin=0; ibin<nybins; ibin++)
        if ((y>=minybox+ibin*ybinsize) && (y<minybox + (ibin+1)*ybinsize)) {biny=ibin; break;}

    for (ibin=0; ibin<nzbins; ibin++)
        if ((z>=minzbox+ibin*zbinsize) && (z<minzbox + (ibin+1)*zbinsize)) {binz=ibin; break;}

    bin = binx + biny*nxbins + binz*(nxbins*nybins);

    return bin;
}

bool NeighborList::initialize_binning() {

    bool success = true;

    for (int i=0; i<totalbins; i++) {
        atomsperbin[i] = 0;
        heads[i] = -1;
    }

    for (int i=0; i<natoms; i++) 
        binlist[i] = -1;

    // #pragma omp parallel for
    for (int atomid=0; atomid<natoms; atomid++) {

        Atom myatom = atomicsystem.get_atom(atomid);
        double x = myatom.get_x();
        double y = myatom.get_y();
        double z = myatom.get_z();
        int binnumber = get_bin_number(x,y,z); //bin number starts at zero
        int j = 0;
        if (heads[binnumber] == -1) {
           heads[binnumber] = atomid;
        }
        else {
            j = heads[binnumber];
            while(binlist[j] != -1) j = binlist[j];
            binlist[j] = atomid;
        }
        atomsperbin[binnumber] += 1; 

    }
    return success;
}

bool NeighborList::is_bin_valid(double x, double y, double z, int currentbin, int nextbin, int totalbins) {

    double distance1 = 0.0, distance2=0.0, distance3=0.0, distance4=0.0;
    double minxbin=0.0, minybin=0.0, minzbin=0.0, maxxbin=0.0, maxybin=0.0, maxzbin=0.0;

    int binsperlayer = nxbins*nybins, ixbin=0, iybin=0, izbin=0;

    bool isvalid = false;

    if (nextbin == currentbin) isvalid = true;

    else {

        izbin = nextbin/binsperlayer;
        iybin = (nextbin - (izbin*binsperlayer))/nxbins;
        ixbin = nextbin - (izbin*binsperlayer) - (iybin*nxbins);

        minxbin = minxbox + ixbin*xbinsize;
        minybin = minybox + iybin*ybinsize;
        minzbin = minzbox + izbin*zbinsize;

        maxxbin = minxbox + (ixbin+1)*xbinsize;
        maxybin = minybox + (iybin+1)*ybinsize;
        maxzbin = minzbox + (izbin+1)*zbinsize;
       
        distance1 = pow(x-minxbin,2);
        distance2 = pow(x-maxxbin,2);
        distance3 = pow(maxxbin-xboxsize-x,2);
        distance4 = pow(minxbin+xboxsize-x,2);

        if (distance1 <= cutoffsq || distance2 <= cutoffsq || distance3 <= cutoffsq || distance4 <= cutoffsq) isvalid = true;

        if (isvalid == false) {

            distance1 = pow(y-minybin,2);
            distance2 = pow(y-maxybin,2);
            distance3 = pow(maxybin-yboxsize-y,2);
            distance4 = pow(minybin+yboxsize-y,2); 

            if (distance1 <= cutoffsq || distance2 <= cutoffsq || distance3 <= cutoffsq || distance4 <= cutoffsq) isvalid = true;
        }

        if (isvalid == false) {

            distance1 = pow(z-minzbin,2);
            distance2 = pow(z-maxzbin,2);
            distance3 = pow(maxzbin-zboxsize-z,2);
            distance4 = pow(minzbin+zboxsize-z,2); 

            if (distance1 <= cutoffsq || distance2 <= cutoffsq || distance3 <= cutoffsq || distance4 <= cutoffsq) isvalid = true;
        }
    }

    return isvalid;

}


void NeighborList::find_neighbors(int atomid) {


    int nneighbors=0, currentbin=0, atomid2;
    double x1, y1, z1, x2,y2,z2, distancesq=0.0;
    bool binisvalid = true;

    int *myneighbors = new int[maxneighbors];
    for (int i=0; i<maxneighbors; i++) myneighbors[i] = -1;
    
    Atom myatom1 = atomicsystem.get_atom(atomid), myatom2;

    x1 = myatom1.get_x();
    y1 = myatom1.get_y();
    z1 = myatom1.get_z();

    currentbin = get_bin_number(x1,y1,z1);
    int maxbins;
    bool visitallbins = false;
    int binsubset[13] = {0,7,8,15,16,23,5,6,13,14,21,22,25}; //only 13 neighboring bins at a time

    if (totalbins < 27) {
        maxbins = totalbins-1;
        visitallbins = true;
    } 
    else{
        maxbins = 13;
    }


    for (int i=0; i<=maxbins; i++) {
        int binid, binpos = binsubset[i];
        if (visitallbins == false) {
            if (i<maxbins) binid = neighboringbins[currentbin][binpos];
            if (i==maxbins) binid = currentbin;
        }
        else { binid = i; //cout<<"Working on bin #"<<i<<"\n";
         }
        if (binisvalid == true && atomsperbin[binid] > 0) {

            int *bin = get_atoms_in_bin(binid);
            
            for (int pos=0; pos<atomsperbin[binid]; pos++) {

                atomid2 = bin[pos];

                if (atomid2 != atomid){
                    distancesq = atomicsystem.get_square_distance(atomid,atomid2);
                  //  if (atomid==0) cout<<"0 "<<atomid2<<" "<<sqrt(distancesq)<<"\n";

                    if (distancesq <= cutoffsq) {
                        bool repeat = false;
                        for (int r=0; r<neighborsperatom[atomid]; r++)
                            if (atomid2 == neighbors[atomid][r]) {
                                repeat = true; 
                                break;
                            }
                        if (repeat == false) {

                            neighbors[atomid][neighborsperatom[atomid]] = atomid2;
                            neighborsperatom[atomid] += 1;

                            neighbors[atomid2][neighborsperatom[atomid2]] = atomid;
                            neighborsperatom[atomid2] += 1;

                        }

                        if (neighborsperatom[atomid] == maxneighbors-1) {
                            cout << "WARNING: maximum number of neighbors reached for atom id " << atomid << " ("<<maxneighbors<< ")\n";
                            goto exitloop;
 
                        }
                    }
                    
                }     
 
            }
        }
    }

    exitloop:     
    int a = 0;
}


void NeighborList::find_neighboring_bins() {

    int jxbin=0, jybin=0, jzbin=0, binnumber;
    int binsperlayer = nxbins*nybins;

    for (int bin=0; bin<totalbins; bin++) {

        int izbin = bin/binsperlayer;
        int iybin = (bin - (izbin*binsperlayer))/nxbins;
        int ixbin = bin - (izbin*binsperlayer) - (iybin*nxbins);

        //add neibhboring bins clockwise starting from the same layer as the current bin, then moving down, then up

        for (int layer=0; layer<3; layer++) {

            if (layer == 0) jzbin = izbin;
            if (layer == 1) jzbin = (izbin-1 >= 0) ? izbin-1 : nzbins-1;
            if (layer == 2) jzbin = (izbin + 1 < nzbins) ? izbin + 1: 0;

            jxbin = (ixbin + 1 < nxbins) ? ixbin +1 : 0;
            binnumber = jxbin + iybin*nxbins + jzbin*(nxbins*nybins);
            neighboringbins[bin][layer*8+0] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            jybin = (iybin - 1 >= 0) ? iybin- 1 : nybins - 1;
            binnumber = jxbin + jybin*nxbins + jzbin*(nxbins*nybins);
            neighboringbins[bin][layer*8+1] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            binnumber = ixbin + jybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+2] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            jxbin = (ixbin-1 >= 0)? ixbin-1 : nxbins-1;
            binnumber = jxbin + jybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+3] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            binnumber = jxbin + iybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+4] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            jybin = (iybin + 1 < nybins) ? iybin + 1 : 0;
            binnumber = jxbin + jybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+5] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            binnumber = ixbin + jybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+6] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";

            jxbin = (ixbin + 1 < nxbins) ? ixbin + 1: 0;
            binnumber = ixbin + jybin*nxbins + jzbin*nxbins*nybins;
            neighboringbins[bin][layer*8+7] = binnumber;

            if (binnumber >= totalbins) cout<<bin<<" "<<layer<<" 0 "<<binnumber<<" "<<ixbin<<" "<<iybin<<" "<<izbin<<" "<<jzbin<<"\n";
            
        }

        jzbin = (izbin-1 >= 0) ? izbin - 1 : nzbins - 1;
        binnumber = ixbin + iybin*nxbins + jzbin*nxbins*nybins;
        neighboringbins[bin][24] = binnumber;

        jzbin = (izbin + 1 < nzbins) ? izbin +1 : 0;
        binnumber = ixbin + iybin*nxbins + jzbin*(nxbins*nybins);
        neighboringbins[bin][25] = binnumber;

    }

}

int* NeighborList::get_neighboring_bins(int binid) {return neighboringbins[binid];}

int* NeighborList::get_neighbors(int id) {return neighbors[id];}

int* NeighborList::get_sorted_neighbors(int id){
     int n = 0;
     int tempid;
     int * sortedneighbors;
     sortedneighbors = new int[neighborsperatom[id]];
     for (int i=0; i<neighborsperatom[id]; i++) sortedneighbors[i] = neighbors[id][i];

     while (n<neighborsperatom[id]){
        tempid = sortedneighbors[n];
        double temp = atomicsystem.get_square_distance(id,sortedneighbors[n]);
        for (int i=n; i<neighborsperatom[id]; i++){        
            double d1 = atomicsystem.get_square_distance(id,sortedneighbors[i]);
            if (d1 < temp) {tempid=sortedneighbors[i];sortedneighbors[i]=sortedneighbors[n];sortedneighbors[n]=tempid;}
        }
        n++;
    }
    return sortedneighbors;
}

int* NeighborList::get_sorted_neighbors(int atomid, string atomtype) { //returns only neighbor atoms that are of the specified type

    vector<int> neigh;
    int *sortedneighbors = get_sorted_neighbors(atomid);

    for (int i=0; i<neighborsperatom[atomid]; i++) {
        Atom myatom = atomicsystem.get_atom(sortedneighbors[i]);
        if (myatom.get_atom_type() == atomtype)
            neigh.push_back(sortedneighbors[i]);        
    }
    int * neighborssubset = new int[neigh.size()];
    for (int i=0; i<neigh.size(); i++) 
        neighborssubset[i] = neigh[i];

    return neighborssubset;
} 

int* NeighborList::get_sorted_neighbors(int atomid, vector<string> atomtypes) { //returns only neighbor atoms that are of the specified type

    vector<int> neigh;
    int *sortedneighbors = get_sorted_neighbors(atomid);
  //  cout<<"Getting neighbors of types: "<<atomtypes[0]<<" and "<<atomtypes[1]<<" out of "<<neighborsperatom[atomid]<<" neighbors\n";
    for (int i=0; i<neighborsperatom[atomid]; i++) {
        Atom myatom = atomicsystem.get_atom(sortedneighbors[i]);
        if ( std::find( atomtypes.begin(), atomtypes.end(), myatom.get_atom_type() ) != atomtypes.end() ){
            neigh.push_back(sortedneighbors[i]);        
           // cout<<sortedneighbors[i]<<" ";
        }
    }
    int * neighborssubset = new int[neigh.size()];
    for (int i=0; i<neigh.size(); i++) 
        neighborssubset[i] = neigh[i];
  //  cout<<"\nTotal found: "<<neigh.size()<<"\n";
    return neighborssubset;
} 

int NeighborList::get_n_neighbors(int atomid, string atomtype) { //returns only neighbor atoms that are of the specified type

    int n = 0;
    for (int i=0; i<neighborsperatom[atomid]; i++) {
        Atom myatom = atomicsystem.get_atom(neighbors[atomid][i]);
        if (myatom.get_atom_type() == atomtype)
           n++;
    }  

    return n;
} 

int NeighborList::get_n_neighbors(int atomid, vector<string> atomtypes){

    int n = 0;
    for (int i=0; i<neighborsperatom[atomid]; i++) {
        Atom myatom = atomicsystem.get_atom(neighbors[atomid][i]);
        if ( std::find( atomtypes.begin(), atomtypes.end(), myatom.get_atom_type() ) != atomtypes.end() )
           n++;
    }  

    return n;
}

int NeighborList::get_atoms_per_bin(int binid) {return atomsperbin[binid];}
int NeighborList::get_n_neighbors(int id) {return neighborsperatom[id];}
