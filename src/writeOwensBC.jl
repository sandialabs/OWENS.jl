function writeOwensBC(fileRoot, nodes, cmkType, cmkValues)
#writeBCFile writes a boundary condition file for the OWENS Toolkit
# **********************************************************************
# *                   Part of SNL VAWTGen                              *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   writeOwensBC(meshout)
#
#   This function writes a boundary condition file for the OWENS Toolkit
#
#      input:
#      fileRoot     = string containing input prefix of file name

#      output:     (NONE)

# if extension is used in the filename, remove it
if exist(strfind(fileRoot,'.'),'var')
    fileRoot = fileRoot(1:find(fileRoot,'.')-1);
end

# open the BC file to save the boundary conditions to
BCfile = [fileRoot '.bc'];    #construct file name
fid = fopen(BCfile,'w');     #open boundary condition file

# first line should list the number of BCs to follow
totalBCs = 0; 
for nn = 1:length(nodes)
    totalBCs = totalBCs + numel(find(cmkValues{nn} !=0));
end 
fprintf(fid, '#i\n', totalBCs);

# write out the boundary conditions into the file
for nn = 1:length(nodes)
    [row, col, val] = find(cmkValues{nn});    
    for ii = 1:length(row)
        fprintf(fid,'#i #s #i #i #.2f\n',nodes(nn), cmkType{nn}, row(ii), col(ii), val(ii));
    end
end

# close boundary condition file
fclose(fid); 

end