import numpy as np
import pdb
import time

def findParent(id, gridSize, debug):

    nIndicesInRefLvl = list()
    for refLvl in np.arange(1,10):
        nIndicesInRefLvl.append(gridSize * 2 ** ((refLvl - 1) * 3))
    
    for i in np.arange(len(nIndicesInRefLvl)):
        if id <= sum(nIndicesInRefLvl[:i+1]):
            refLvl = i
            break
    if refLvl == 0:
        if id > 0:
            print("cell {:3d}".format(id)+" does not have a parent")
            pass
            
        return 0, refLvl       
    
    id2 = id - sum(nIndicesInRefLvl[:refLvl])
    ix = (id2 - 1) % (xdim * 2 ** refLvl) + 1
    iy = (id2 - 1) / (xdim * 2 ** refLvl) % (ydim * 2 ** refLvl) + 1
    iz = (id2 - 1) / (xdim * 2 ** refLvl * ydim * 2 ** refLvl) + 1
    parentId = (int(np.ceil(iz / 2.0) - 1) * xdim * 2 ** (refLvl - 1) * ydim * 2 ** (refLvl - 1) +
                int(np.ceil(iy / 2.0) - 1) * xdim * 2 ** (refLvl - 1) +
                int(np.ceil(ix / 2.0)) +
                sum(nIndicesInRefLvl[:refLvl-1]))
    if debug:
        print("id = {:3d}".format(id)+", id2 = {:3d}".format(id2)+
              ", col = {:2d}".format(ix)+", row = {:2d}".format(iy)+
              ", plane = {:2d}".format(iz)+", parentId = {:2d}".format(parentId)+
              ", refLvl = {:1d}".format(refLvl))
    else:
        print("cell {:3d}".format(id)+" is the child of cell {:2d}".format(parentId))
        pass

    return parentId, refLvl

def getChildren(children, parentId, dimension = 0, up = True, left = True):

    down  = not up
    right = not left
        
    N = 8

    myChildren = list()

    # Select 2/8 children per parent according to the logical parameters up, down, left, right.
    # The names depict sides of the four children seen when looking along the direction of the
    # pencil.

    #  ----    ---- 
    # /   /|  /   /| 
    # ---- |  ---- | 
    # |UU| |  |LR| | 
    # |DD|/   |LR|/  
    # ----    ----    
    #
    if dimension == 0:
        if up and left:
            i1 = 0
            i2 = 1
        if down and left:
            i1 = 2
            i2 = 3
        if up and right:
            i1 = 4
            i2 = 5
        if down and right:
            i1 = 6
            i2 = 7

    if dimension == 1:
        if up and left:
            i1 = 0
            i2 = 2
        if down and left:
            i1 = 1
            i2 = 3
        if up and right:
            i1 = 4
            i2 = 6
        if down and right:
            i1 = 5
            i2 = 7

    if dimension == 2:
        if up and left:
            i1 = 0
            i2 = 4
        if down and left:
            i1 = 1
            i2 = 5
        if up and right:
            i1 = 2
            i2 = 6
        if down and right:
            i1 = 3
            i2 = 7

    if parentId in children.keys():
        myChildren.extend(children[parentId][i1::N])
        myChildren.extend(children[parentId][i2::N])
    else:
        # If no children were found, return the parent
        myChildren.extend(parentId)

    #print(up,left,myChildren)
    return myChildren

def buildPencils(pencils,initialPencil,idsIn,dimension = 0,path = list()):

    # pencils - list of completed pencils
    # initalPencil - list of ids that have already been added to the pencil being built
    # idsIn - candidate cell ids to be added to the pencil being built (unless they contain refinement)
    # dimension - dimension along which the pencils are built
    # path - the steps (up/down, left/right) taken while building the current pencil

    # Global arrays that are accessed read-only
    # isRefined - global array that contains how many times each cell has been refined
    # refLvls - global array that contains the refinement level of each cell
    # children - global array that contains the children of each refined cell
    import copy

    # (Hard) Copy the input ids to a working set of ids
    ids = copy.copy(idsIn)

    # (Hard) Copy the already computed pencil to the output list
    idsOut = copy.copy(initialPencil)
    
    # Walk along the input pencil
    for i,id in enumerate(ids):

        i1 = i + 1
        # Check if the current cell contains refined cells
        if isRefined[id] > 0:
            
            # Check if we have encountered this refinement level before and stored
            # The path this builder followed
            if len(path) > refLvls[id]:

                # Get children using the stored path
                myChildren = getChildren(children,id,dimension,
                                         path[refLvls[id]][0],path[refLvls[id]][1])
                
                # Add the children to the working set
                ids[i1:i1] = myChildren
                
            else:

                # Spawn new builders to construct pencils at the new refinement level
                for up in [True, False]:
                    for left in [True, False]:
                        
                        # Store the path this builder has chosen
                        myPath = copy.copy(path)
                        myPath.append((up,left))

                        # Get children along my path
                        myChildren = getChildren(children,id,dimension,up,left)
                        myIds = ids[i1:]

                        # The current builder will continue along the bottom-right path
                        if not up and not left:
                            
                            # Add the children to the working set. Next iteration of the
                            # main looop (over ids) will start on the first child
                            ids[i1:i1] = myChildren
                            path = myPath
                            #print('building pencil for'+str(ids[i1:]))
                            pass

                        # Other paths will spawn a new builder
                        else:

                            # Create a new working set by adding the remainder of the old
                            # working set to the current children.
                            myChildren.extend(myIds)
                            
                            buildPencils(pencils,idsOut,myChildren,dimension,myPath)

        # Add unrefined cells to the pencil directly
        else:

            idsOut.append(id)
            
            pass

    pencils.append(idsOut)
    return pencils

    #print(idsOut)
    #print(pencils)

import argparse

parser = argparse.ArgumentParser(description='Create pencils on a refined grid.')
parser.add_argument('--dimension', metavar = 'N', type=int, nargs=1,
                    default=[0], help='Dimension (x = 0, y = 1, z = 2)')
parser.add_argument('--filename', metavar = 'fn', type=str, nargs=1,
                    default=['test.vtk'], help='Input vtk file name')
parser.add_argument('--debug', metavar = 'd', type=int, nargs=1,
                    default=[0], help='Debug printouts (no = 0, yes = 1)')
args = parser.parse_args()

if args.dimension[0] > 0 and args.dimension[0] <= 2:
    dimension = args.dimension[0]
else:
    dimension = 0

debug = bool(args.debug[0])

#filename = 'test.vtk'
filename = args.filename[0]
fh = open(filename)
lines = fh.readlines()
fh.close()

ids = list()

xdim = 1
ydim = 1
zdim = 1
for i,line in enumerate(lines):
    if 'DATASET UNSTRUCTURED_GRID' in line:
        n = int(lines[i+1].split()[1])
        for j in np.arange(n):
            xyz = lines[i+j+2].split()
            xdim = max(xdim,float(xyz[0]))
            ydim = max(ydim,float(xyz[1]))
            zdim = max(zdim,float(xyz[2]))
    if 'SCALARS id int' in line:
        n = int(lines[i-1].split()[1])
        for j in np.arange(n):
            ids.append(int(lines[i+j+2]))

xdim = int(xdim)
ydim = int(ydim)
zdim = int(zdim)

print('grid dimensions are {:2d} x {:2d} x {:2d}'.format(xdim,ydim,zdim))
gridSize = xdim*ydim*zdim

#debug = True

t1 = time.time()

parents = dict()
children = dict()
refLvls = dict()
hasChildren = list()

for id in ids:

    # Find the parent of cell id
    parentId, refLvl = findParent(id,gridSize,debug)
    
    parents[id] = parentId
    refLvls[id] = refLvl

    # Parents are not stored in the id array by default, let's add them
    # For completeness
    if not parentId in ids and parentId > 0:
        ids.append(parentId)

    # Make a list of cells that have been refined at least once
    if parentId > 0:
        if not parentId in hasChildren:
            children[parentId] = list()
            hasChildren.append(parentId)

        # Make a list of children for each cell        
        children[parentId].append(id)

# Sort the id and children lists, this is needed when adding cells to pencils
# to get the order right
for key in children.keys():
    children[key].sort()
ids.sort()

# Second pass to count how many times each cell has been refined
isRefined = dict()
for id in ids:
    isRefined[id] = 0
    if refLvls[id] > 0:
        parentId = parents[id]
        while parentId is not 0:
            isRefined[parentId] = refLvls[id] - refLvls[parentId]
            parentId = parents[parentId]

# Begin sorting, select the dimension by which we sort
# dimension = 0
# dimensions = ('x','y','z')
print
print('Building pencils along dimension {:1d}'.format(dimension))
print

#sortedIds = list()
mapping = dict()
for id in ids:
    # Sort the unrefined mesh ids following Sebastians c++ code
    if dimension == 0:

        dims = (zdim, ydim, xdim)
        
        idMapped = id
    
    if dimension == 1:

        dims = (zdim, xdim, ydim)
        
        x_index = (id-1) % xdim
        y_index = ((id-1) / xdim) % ydim
        idMapped = id - (x_index + y_index * xdim) + y_index + x_index * ydim
        
    if dimension == 2:

        dims = (ydim, xdim, zdim)
        
        x_index = (id-1) % xdim
        y_index = ((id-1) / xdim) % ydim
        z_index = ((id-1) / (xdim * ydim))
        idMapped = 1 + z_index + y_index * zdim + x_index * ydim * zdim

    if refLvls[id] == 0:
        mapping[idMapped] = id

# Create pencils of unrefined cells, store the level of refinement for each cell
unrefinedPencils = list()
for i in np.arange(dims[0]):
    for j in np.arange(dims[1]):
        ibeg = 1 + i * dims[2] * dims[1] + j * dims[2]
        iend = 1 + i * dims[2] * dims[1] + (j + 1) * dims[2]
        myIsRefined = list()
        myIds = list()
        for k in np.arange(ibeg,iend):
            myIds.append(mapping[k])
            myIsRefined.append(isRefined[mapping[k]])
        unrefinedPencils.append({'ids' : myIds,
                                 'refLvl' : myIsRefined})

# Refine the unrefined pencils that contain refined cells

pencils = list()

# Loop over the unrefined pencils
for unrefinedPencil in unrefinedPencils:

    pencils = buildPencils(pencils,[],unrefinedPencil['ids'],dimension)
    
t2 = time.time()

print('I have created the following pencils:')
print
for pencil in pencils:
    print(pencil)
    
print
print('Execution time was {:.4f} seconds'.format(t2-t1))
