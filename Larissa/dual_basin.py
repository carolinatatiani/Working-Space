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


#Calculate the dihedrals and their differences
dih1 = np.rad2deg(md.compute_dihedrals(pdb1, dihedrals_indices)).T
dih2 = np.rad2deg(md.compute_dihedrals(pdb2, dihedrals_indices)).T
dihd = gd.diff (dih1, dih2)

#Check and print the dihedrals that the diference is less than 54.5 degrees
# This is the threshold used in OpenSMOG to select dihedrals for the dual basin model
# The threshold is 54.5 degrees, which is approximately 1.0 radians
common= dihedrals[np.where(dihd<54.5)[0]]

# Remove duplicates from the common dihedrals
unique_rows, idx = np.unique(common[:, :4], axis=0, return_index=True)
common = common[idx]

common[:,:4] = common[:,:4].astype(int) + 1  # Convert to 1-based indexing
np.save('dihcommon.npy', common)

tmp1, tmp2 ,tmp3,du1 = gd.parse_top_file(top1)
tmp1,tmp2,tmp3,du2 = gd.parse_top_file(top2)
# Select the unique dihedrals from both topologies
# and remove the dihedrals that are not in the common list
du1=du1[np.where(dihd>54.5)[0]]
unique_rows, idx = np.unique(du1[:, :4], axis=0, return_index=True)
du1 = du1[idx]
du1[:,:4] = du1[:,:4].astype(int) + 1 # Convert to 1-based indexing

du2=du2[np.where(dihd>54.5)[0]]
unique_rows, idx = np.unique(du2[:, :4], axis=0, return_index=True)
du2 = du2[idx]
du2[:,:4] = du2[:,:4].astype(int) + 1 # Convert to 1-based indexing

np.save('dihunique1.npy', du1)
np.save('dihunique2.npy', du2)

# Read XML of each model to  select the unique interactions
xml1 = write_xml.parse_xml(f'{ref}.xml')
xml2 = write_xml.parse_xml(f'{obj}.xml')

c, u1, u2 = write_xml.classify_contacts(xml1, xml2)

np.save("common", c)
np.save("closed", u1)
np.save("open", u2)


# Create root element
root = ET.Element("OpenSMOGforces")

# Add interactions to the root element
write_xml.write_single_gaussian(root, u1, 'unique-closed', 0.625)
write_xml.write_double_gaussian(root, c, 'common-contacts')
write_xml.write_single_gaussian(root, u2, 'unique-open', 1.375)
write_xml.write_cosine_squared(root, common, 'common-dihedrals', w=1.0)
write_xml.write_cosine_squared(root, du1, 'unique-closed-dihedrals', w=1)
write_xml.write_cosine_squared(root, du2, 'unique-open-dihedrals', w=1)
#write_single_gaussian(root, uniqueC, 'unique-closed', 1)#len(uniqueC)/len(u2)*1.36)#*1.)

# Finalize the XML tree
tree = ET.ElementTree(root)
ET.indent(tree, '  ')
tree.write(output_filename, encoding="ISO-8859-1", xml_declaration=True)

# Validate the XML file
with open(output_filename, 'rb') as f:
    xml_bytes = f.read()

if write_xml.validate_xml(xml_bytes, xsd_file):
    print("XML is valid against the XSD schema.")
else:
    print("XML is NOT valid against the XSD schema.")


# top2 = 'opened.AA.top'
# out_file = 'opened_angles_dihedrals.txt'

# ang_ind = np.load('sel_ang_1based.npy')
# dih_ind = np.load('sel_dih_1based.npy')

# wd.extract_lines(top2, ang_ind, dih_ind, out_file)