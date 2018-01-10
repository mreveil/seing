import os


def read_coords(filename):
    coords = []
    with open(filename,'r') as inFile:
        for line in inFile:
            data = line.strip().split()
            if len(data) == 4 and '{' not in line:
                coords.append(data)
    return coords

def writexyz(coords,filename):
    with open(filename,'w') as outfile:
        outfile.write(str(len(coords))+"\n\n")
        for c in coords:
            outfile.write(c[0]+"\t"+c[1]+"\t"+c[2]+"\t"+c[3]+"\n")

    return

#move to directory and write xyz file

directory = '/fs/home/mr937/Diffusion+MachineLearning/data/FORCE_CALCULATIONS/'
i = 0
print "Starting preprocessing with SEING..."
for root, subdirs, files in os.walk(directory):
    for filename in files: 
        
        if 'scf.in' in filename:
            i += 1
            print "Processing file #",i
            coords = read_coords(os.path.join(root,filename))
            writexyz(coords,os.path.join(root,"coordinates.xyz"))

            ## call seing
            os.chdir(root)
            print "Current directory:",os.getcwd()
            os.system("/fs/home/mr937/SEING/bin/seing coordinates.xyz fingerprints.sg")

#Go through folders again read forces and fingerprints
all_fingerprints, all_forces = [],[]
for root, subdirs, files in os.walk(directory):
    for filename in files: 
        
        if 'fingerprints' in filename:
            fingerprints += readfingerprints(filename);
        if 'scf.out' in filename:
            all_forces += read_forces(os.path.join(root,filename))

#
#write CSV file
