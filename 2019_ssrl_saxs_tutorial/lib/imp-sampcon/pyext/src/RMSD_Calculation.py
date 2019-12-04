
try:
    import pyRMSD.RMSDCalculator
    #from pyRMSD.matrixHandler import MatrixHandler
    from pyRMSD.condensedMatrix import CondensedMatrix
    pyR=True
except:
    pyR=False

import numpy as np
import sys, os, glob

import IMP
import IMP.atom
import IMP.rmf
import RMF

def get_pdbs_coordinates(path, idfile_A, idfile_B):

    pts = []
    conform = []
    num = 0
    masses = []
    radii = []
    
    models_name = []
     
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')
    for str_file in sorted(glob.glob("%s/analysis/sample_A/*.pdb" % path),key=lambda x:int(x.split('/')[-1].split('.')[0])):
        print >>f1, str_file, num
        models_name.append(str_file)
        
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh) 
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        if num == 0:
            masses = [IMP.atom.Mass(p).get_mass() for p in mps]
            radii  = [IMP.core.XYZR(p).get_radius() for p in mps]
        conform.append(pts)

        pts = []
        num = num + 1

        
    for str_file in sorted(glob.glob("%s/sample_B/*.pdb" % path),key=lambda x:int(x.split('/')[-1].split('.')[0])):
        print >>f2, str_file, num
        models_name.append(str_file)
        
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh)
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        conform.append(pts)
        pts = []   
        num = num + 1
        
    return np.array(conform), masses, radii, models_name

def get_rmfs_coordinates(path, rmf_A, rmf_B, subunit_name):

    rmf_fh = RMF.open_rmf_file_read_only(rmf_A)
    n_models = rmf_fh.get_number_of_frames()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_B)
    n_models += rmf_fh.get_number_of_frames()
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    IMP.rmf.load_frame(rmf_fh, 0)
    m.update()
    pts = 0

    if subunit_name:
        s0 = IMP.atom.Selection(h, resolution=1,molecule=subunit_name)
    else:
        s0 = IMP.atom.Selection(h, resolution=1)


    for leaf in s0.get_selected_particles():
        p=IMP.core.XYZR(leaf)
        pts+=1
    
    conform = np.empty([n_models,pts,3])

    num = 0
    masses = []
    radii = []
    ps_names = []
    n_models = [] 
    m = IMP.Model()
    models_name = []
    mod_id = 0
    for rmf_file in [rmf_A, rmf_B]:
        #print >>rmf_file
        models_name.append(rmf_file)

        rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
        h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
        n_models.append(rmf_fh.get_number_of_frames())
        print("Opening RMF file:", rmf_file, "with", n_models[-1], "frames")    
        for f in range(rmf_fh.get_number_of_frames()):
            if f%100==0:
                pass
                #print("  -- Opening frame", f, "of", rmf_fh.get_number_of_frames())
            IMP.rmf.load_frame(rmf_fh, f)
            m.update()
            pts = 0

            if subunit_name:
                s0 = IMP.atom.Selection(h, resolution=1,molecule=subunit_name)
            else:
                s0 = IMP.atom.Selection(h, resolution=1)
                
            
            for leaf in s0.get_selected_particles():
                
                p=IMP.core.XYZR(leaf)
                pxyz = p.get_coordinates()
                conform[mod_id][pts][0] = pxyz[0]
                conform[mod_id][pts][1] = pxyz[1]
                conform[mod_id][pts][2] = pxyz[2]
                pts+=1
                
                if num == 0 and rmf_file==rmf_A:
                    masses.append(IMP.atom.Mass(leaf).get_mass())
                    radii.append(p.get_radius())
                    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
                    copy_number = "X"
                    # Need to find the copy number from the molecule
                    # In PMI, this is three levels above the individual residues/beads
                    mol_p = IMP.atom.Hierarchy(p).get_parent().get_parent().get_parent()
                    if IMP.atom.Copy().get_is_setup(mol_p):
                        copy_number = str(IMP.atom.Copy(mol_p).get_copy_index())
                    
                    if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
                        residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()
                        
                        ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead))+"_"+copy_number)
                            
                    else:
                        residue_in_bead = str(IMP.atom.Residue(leaf).get_index())
                        
                        ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead+"_"+copy_number)
            mod_id+=1
            #conform.append(pts)
            pts = []
            num = num + 1
        del rmf_fh        
    return ps_names, masses, radii, conform, models_name, n_models

def get_rmsds_matrix(conforms, mode, sup, cores):
    
    if pyR:
        print("Mode:",mode,"Superposition:",sup,"Number of cores:",cores)

        if(mode=="cpu_serial" and not sup):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)

        elif(mode=="cpu_omp" and not sup):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)
            calculator.setNumberOfOpenMPThreads(int(cores))

        elif(mode=="cpu_omp" and sup):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", conforms)
            calculator.setNumberOfOpenMPThreads(int(cores))

        elif(mode=="cuda" and sup):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_MEM_CALCULATOR", conforms)  

        else:
            print("Wrong values to pyRMSD ! Please Fix")
            exit()

        rmsd = calculator.pairwiseRMSDMatrix()
        rmsd_matrix=CondensedMatrix(rmsd)
        inner_data = rmsd_matrix.get_data()

    else:
        inner_data = np.empty(int((len(conforms)*len(conforms)-len(conforms))/2))
        print(inner_data.shape)
        # Create list of upper-triangleNeed upper 
        l = 0
        for x in range(len(conforms)):
            if x%10==0:
                print("Conform", x)
            for y in range(x+1, len(conforms)):

                d1 = []
                d0 = []
                for i in range(len(conforms[x])):
                    #print(conforms[x][i])
                    d1.append(IMP.algebra.Vector3D(conforms[x][i]))
                    d0.append(IMP.algebra.Vector3D(conforms[y][i]))

                xform = IMP.algebra.get_transformation_aligning_first_to_second(d1, d0)
                rmsd = IMP.algebra.get_rmsd_transforming_first(xform, d1, d0)
                inner_data[l]=rmsd
                l+=1

    np.save("Distances_Matrix.data", inner_data)

    return inner_data
