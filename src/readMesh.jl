mutable struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
end

function readMesh(filename)
#readMesh  reads mesh file and stores data in mesh object
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [mesh] = readMesh(filename)
#
#   This function reads the mesh file and stores data in the mesh object.
#
#      input:
#      filename      = string containing mesh filename

#      output:
#      mesh          = object containing mesh data

fid = open(filename,"r")   #open mesh file

# temp = fscanf(fid,'#i',2)   #read in number of nodes and number of elements
line = readline(fid)
temp = split(line)

numNodes = parse(Int,temp[1])
numEl = parse(Int,temp[2])

nodeNum = zeros(numNodes,1)
x = zeros(numNodes,1)
y = zeros(numNodes,1)
z = zeros(numNodes,1)

conn = zeros(numEl,2)
elNum = zeros(numEl,1)

for i=1:numNodes            # read in node number and node coordinates
    line = readline(fid)
    temp = split(line)
    nodeNum[i] = parse(Float64,temp[1])
    x[i] = parse(Float64,temp[2])
    y[i] = parse(Float64,temp[3])
    z[i] = parse(Float64,temp[4])
end

for i=1:numEl               # read in element number and connectivity list
    line = readline(fid)
    temp = split(line)
    elNum[i] = parse(Float64,temp[1])

    conn[i,1] = parse(Float64,temp[3])
    conn[i,2] = parse(Float64,temp[4])
end

close(fid)  #close mesh file

mesh = Mesh(nodeNum,
numEl,
numNodes,
x,
y,
z,
elNum,
conn)

return mesh

end
