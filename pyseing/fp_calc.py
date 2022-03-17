from _pyseing import *
import time

def CalculateFP(xyzfilename, propfilename):
    start = time.time()
    fpproperties = read_prop_file(propfilename)
    asys = AtomicSystem(xyzfilename, True, True, True, 2.01778)
    print(asys.get_n_atoms())
    print(fpproperties.box_size)
    # xmin = fpproperties.box_size[0] 
    # ymin = fpproperties.box_size[1] 
    # zmin = fpproperties.box_size[2] 
    # xmax = fpproperties.box_size[3] 
    # ymax = fpproperties.box_size[4] 
    # zmax = fpproperties.box_size[5]

    fpgen = FingerprintGenerator(asys, fpproperties)
    print('Wrinting fingerprint to to file...', fpproperties.output_file+'\n')
    fpgen.write2file(fpproperties.output_file, fpproperties.output_mode)
    end = time.time()
    print('Calculation time =', str(end - start))

