
    variable _pathGRO "/home/bio/Desktop/jupyter_test/antibody_test/output_20241126_112516/configuration/cycle1_MD"
    variable _FileNameGRO "system_Compl_MD"
    variable _pathPDB "/home/bio/Desktop/jupyter_test/antibody_test/output_20241126_112516/configuration/REPOSITORY"
    variable _FileNamePDB "HLA_BiAB_protein_50ns"
    variable _FileNamePDB_OUT "cycle1_starting_protein"
    
# This script will load a GRO file, strip water and ions, save it as a new_PDB file. Then upload the starting PDB file
# with the chain information, load the new_PDB file on it and save it again. This will update the chain information
# on the new_PDB

# I NEED THE FOLLOWING VARIABLES TO BE SET:
#       variable _pathGRO path
#       variable _FileNameGRO name
#       variable _pathPDB path
#       variable _FileNamePDB name
#       variable _FileNamePDB_OUT name

mol new ${_pathGRO}/${_FileNameGRO}.gro
set sel_protein [atomselect top "not water and not ions"]
$sel_protein writepdb ${_pathGRO}/temp_protein.pdb
mol delete 0

mol new ${_pathPDB}/${_FileNamePDB}.pdb
mol addfile ${_pathGRO}/temp_protein.pdb

set sel_protein [atomselect top "all" frame 1]
$sel_protein writepdb ${_pathGRO}/${_FileNamePDB_OUT}.pdb

puts "PDB file saved: ${_pathGRO}/${_FileNamePDB_OUT}.pdb "
exit
