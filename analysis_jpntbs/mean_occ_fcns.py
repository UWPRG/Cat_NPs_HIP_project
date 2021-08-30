#Protein/polmyer analysis functions

import MDAnalysis as mda
import MDAnalysis.analysis.distances as maa_dist
import statsmodels as stats
import math
import numpy as np
import pandas as pd


def aa_frmcount(prot_atoms, g2_atoms, dmax, universe, start, end):
    """This function will output a dictionary of AA protein residue number and its corresponding frame count and occupancy. 
    
    Inputs consisted of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms selected ONLY from the Universe 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups 
    - universe: MD Analysis variable containing information for the enitre trajectory
    - start: Start of the trajectory (first frame number)
    - end: End of the trajectory (last frame number)
    
    Output is consisted of: 
    
    - aa_dict: A dictionary whose keys are residues that contacted group two atoms within the trajectory and
               values are an array containing the frame count and occupancy. Here, occupancy is defined as the 
               frame count of each residue divided by the total number of frames in trajectory block. Frame count 
               is the number of frames that a group two atom was less than dmax distance from any protein atom. 
    """
    
    aa_dict = {} # Initialize a dictionary to store the output. 
    
    #laa = np.zeros(shape=len(universe.select_atoms("protein").residues)) 
    # Initialize a numpy array filled with zeros whose dim is no. of protein residues 
    laa = np.zeros(shape=len(prot_atoms.residues))
    
    # Store protein residues as given from MD Analysis
    #br = np.array(universe.select_atoms("protein").residues)  
    br = np.array(prot_atoms.residues)
    
    # This for loop goes through the trajectory block
    for ts in universe.trajectory[start:end]: 
        
        count = 0  # start counter at zero 
        
        # Use the get_protresd_list function to get residues that contact the group two atoms 
        bsres = get_protresd_list(prot_atoms, g2_atoms, dmax, universe) 
        
        # If the output array is empty, then no residues were within 4 Angstroms of group two atoms 
        if bsres.size == 0: 
            pass
        elif bsres.size != 0: # If array is not empty
    
            count += 1 # add to counter
            
            # Go through each residue 
            for i in bsres.flat:
        
                res_ind = np.where(br == i)[0]  # Get index for stored residue array that matches the residues from output in get_protresd_list
            
                laa[res_ind[0]] = laa[res_ind[0]] + count # update frame count to the specific residure 
                
    fin_res = np.where(laa != 0)[0] # Get indexes for residues with frame counts that are not zero
    
    # Calculate occupancy for each residue and store in a dictionary 
    for i in fin_res.flat:
        aa_dict[str(list(prot_atoms.residues[i:i+1])[0])] = [laa[i], laa[i]/(end - start)]
    
    return aa_dict 


def AA_list_org(lorg_list):
    
    """List elements need have 'GLY  XX' as string format, where XX reps the number of GLY residues. Output is a
    sorted list of 'AA XX' according to the below order.  """
    
    aromatic_res = ['PHE', 'TRP', 'TYR', 'HID','HIE']
    hydrophobic_res = ['ALA', 'ILE', 'LEU', 'VAL', 'PRO','PHE', 'TRP','MET','TYR']
    polar_res = ['ASN', 'CYS','CYX','ASH', 'GLH','GLN', 'SER', 'THR','GLY','HIE','HID']
    neg_res = ['ASP', 'GLU']
    pos_res = ['ARG', 'HIP', 'LYS']

    all_res = [pos_res, neg_res, polar_res, hydrophobic_res, aromatic_res]
    #Change order of residues before making the bar graph
    # (1) Positively charged
    # (2) Negatively charged
    # (3) Polar residues 
    # (4) Hydrophobic residues 
    
    # This chunk of code sorts the counts of each AA that have 1001 or 1002 frame count based 
    # on the AA order in all_res
    arr_list = []

    for row in all_res:
        for i in range(len(lorg_list)):
            for j in range(len(row)):
                if lorg_list[i][0:3].find(row[j]) != -1:
                    arr_list.append(lorg_list[i])
                    
    #This chunk of code splits the list arr_list to makes the AA: count of 1001 or 1002 frames data plottable 
    f_list = []
    fn_list = []
    for i in range(len(arr_list)):
        f_list.append(arr_list[i][0:3])
        fn_list.append(int(arr_list[i][5:]))
        
    return f_list, fn_list


def read_xvg(file_name):
    with open(file_name) as f:
        f_lines = []
        names = []
        for ln in f:
            if ln.startswith("#"):
                pass
            elif ln.startswith("@ s"):
                name_line = ln[0:]
                fix_line = name_line.split('"')[1::2]
                names.append(fix_line[0])
            elif ln.startswith("@"):
                pass
            else:
                f_lines.append(ln[0:])
    count = len(open(file_name).readlines())
    skip_rows = count-len(f_lines)
    df = pd.read_csv(file_name, delim_whitespace=True, names=names, skiprows=skip_rows)
    return df
