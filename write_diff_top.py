import numpy as np

def extract_lines(top_file, ang_ind, dih_ind, out_file):

    with open(top_file, 'r') as file:
        lines = file.readlines()
    

    angles = {tuple(ang) for ang in ang_ind}
    dihedrals = {tuple(dih) for dih in dih_ind}
    

    extracted_angles = []
    extracted_dihedrals = []
    inside_section = None
    i = 0
    total_lines = len(lines)

    while i < total_lines:
        line = lines[i].strip()
        
   
        if line.startswith("[ angles ]"):
            inside_section = "angles"
            i += 1
            continue
        elif line.startswith("[ dihedrals ]"):
            inside_section = "dihedrals"
            i += 1
            continue
        elif line.startswith("["):
            inside_section = None

        if inside_section == "angles" and not line.startswith(';'):
            data = line.split()
            if len(data) >= 4:
                current_angle = tuple(int(data[j]) for j in range(3))
                if current_angle in angles:
                    extracted_angles.append(line)

        elif inside_section == "dihedrals" and not line.startswith(';'):
              data = line.split()
              if len(data) >= 5:
                  current_dihedral = tuple(int(data[j]) for j in range(4))
                  funct = int(data[4])
                  if current_dihedral in dihedrals:
                      extracted_dihedrals.append(line)
                      if funct == 1 and i + 1 < total_lines:
                          next_line = lines[i + 1].strip()
                          if not next_line.startswith(';'):
                              extracted_dihedrals.append(next_line)
                              i += 1  
        i += 1  
    
   
    with open(out_file, 'w') as file:
        file.write('[ angles ]\n')
        file.writelines(f"{line}\n" for line in extracted_angles)
        file.write('[ dihedrals ]\n')
        file.writelines(f"{line}\n" for line in extracted_dihedrals)


# top_file = 'opened.AA.top'
# out_file = 'opened_angles_dihedrals.txt'

# ang_ind = np.load('sel_ang_1based.npy')
# dih_ind = np.load('sel_dih_1based.npy')

# extract_lines(top_file, ang_ind, dih_ind, out_file)
