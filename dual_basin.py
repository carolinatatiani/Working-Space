# by Carolina Tatiani
# email: carolina.tatiani@unesp.br
# Last Modified: 03/07/2025

import write_xml_dual_basin2 as write_xml
import get_diff as gd
import write_diff_top as wd
import xml.etree.ElementTree as ET
from lxml import etree
import numpy as np
import mdtraj as md 

output_filename = "dual_basin.xml"
xsd_file = "OpenSMOG.xsd"

ref='closed.AA'
obj='opened.AA'

# Load the topologies
top1 = ref + '.top'
top2 = obj + '.top'

# Load the PDB files
pdb1 = md.load(ref + '.pdb')
pdb2 = md.load(obj + '.pdb')

# Parse the topology files to get angles and dihedrals indices
angles_indices, dihedrals_indices,angles,dihedrals = gd.parse_top_file(top1)


print(angles_indices[0])
print(angles[0])
#Compute the angles and dihedrals for both PDB files
# and calculate the differences using the indices from the topology file

ang1 = np.rad2deg(md.compute_angles(pdb1, angles_indices)).T
ang2 = np.rad2deg(md.compute_angles(pdb2, angles_indices)).T
angd = gd.diff(ang1, ang2)


#Check and print the angles that are greater than 20% of the second angle

indices_ang_percent = []
for i in range(len(ang1)):
   if np.any(ang1[i] > 1.5 * ang2[i]):
      indices_ang_percent.append(i)
np.save('sel_ang_1based.npy', angles_indices[indices_ang_percent])
np.save('angvalues.npy', angles[indices_ang_percent])


# Delete duplicate indices, convert from radians to degrees, and compute the difference 
#dihedrals_indices= np.unique(dihedrals_indices, axis=0)
dih1 = np.rad2deg(md.compute_dihedrals(pdb1, dihedrals_indices)).T
dih2 = np.rad2deg(md.compute_dihedrals(pdb2, dihedrals_indices)).T
dihd = gd.diff (dih1, dih2)

print(len(dihedrals_indices))
print(len(dihedrals))

#Check and print the dihedrals that are greater than 54.5 degrees
print(len(np.where(dihd>54.5)[0]))
ind_dih=dihedrals_indices[np.where(dihd>54.5)[0]]
np.save('sel_dih_1based.npy',ind_dih+1) # +1 to convert from 0-based to 1-based indexing




# # Read XML of each model to  select the unique interactions
# xml1 = write_xml.parse_xml(f'{ref}.xml')
# xml2 = write_xml.parse_xml(f'{obj}.xml')

# c, u1, u2 = write_xml.classify_contacts(xml1, xml2)

# np.save("common", c)
# np.save("closed", u1)
# np.save("open", u2)


# # Create root element
# root = ET.Element("OpenSMOGforces")

# # Add interactions to the root element
# write_xml.write_single_gaussian(root, u1, 'unique-closed', 0.625)
# write_xml.write_double_gaussian(root, c, 'common-contacts')
# write_xml.write_single_gaussian(root, u2, 'unique-open', 1.375)
# #write_single_gaussian(root, uniqueC, 'unique-closed', 1)#len(uniqueC)/len(u2)*1.36)#*1.)

# # Finalize the XML tree
# tree = ET.ElementTree(root)
# ET.indent(tree, '  ')
# tree.write(output_filename, encoding="ISO-8859-1", xml_declaration=True)

# # Validate the XML file
# with open(output_filename, 'rb') as f:
#     xml_bytes = f.read()

# if write_xml.validate_xml(xml_bytes, xsd_file):
#     print("XML is valid against the XSD schema.")
# else:
#     print("XML is NOT valid against the XSD schema.")


# top2 = 'opened.AA.top'
# out_file = 'opened_angles_dihedrals.txt'

# ang_ind = np.load('sel_ang_1based.npy')
# dih_ind = np.load('sel_dih_1based.npy')

# wd.extract_lines(top2, ang_ind, dih_ind, out_file)