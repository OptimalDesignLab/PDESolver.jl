# SimpleMeshFunction2.jl
# Created for the new simpleMesh type

@doc """
### createMesh

It creates a simple mesh with triangle elements in a rectangular domain in 2D. In This code The horizontal axis is referred to as X-axis and the vertical axis is referred to as the Y-axis.

**Inputs**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nnpe`    : Number of nodes per element edge (same for all edges)
*  `numDofPerNode` : Number of degrees of freedom per node

**Outputs**

*  `xnodes` : Number of nodes along X-axis in one row
*  `ynodes` : Number of nodes along Y-axis in one column
*  `nnp`    : Number of nodal points
*  `nel`    : Number of elements mesh
*  `dofs`   : Element information matrix

We assume a counter-clockwise numbering scheme for the the elements and edges with the numbering starting from the base of the geometry.
"""->

function createMesh(lengthx,lengthy,nedx,nedy,nnpe,numDofPerNode)
  xnodes = (nnpe-1)*nedx + 1 # Total number of nodes on all adjacent element edges along X-axis
  ynodes = (nnpe-1)*nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
  nel = 2*nedx*nedy # number of elements
  nsd = 2  # number of spatial dimensions

  if nnpe == 2
    nen = 3 # Number of element nodes
    dofs = Array{Int}(numDofPerNode, nen, nel)
    nnp = xnodes*ynodes # number of nodal points
    m = 0
    for i = 1:2:nel
      for j = 1:numDofPerNode
        dofs[j,1,i] = m + j
      end
      #
      for j = 1:numDofPerNode
        dofs[j,2,i] = m + numDofPerNode + j
      end
      #m += numDofPerNode*xnodes
      for j = 1:numDofPerNode
        dofs[j,3,i] = m + numDofPerNode*xnodes + j
      end
      m += numDofPerNode
      for k = 1:nedy
        if i == 2*nedx*k - 1
          m +=numDofPerNode
        end
      end
    end
    m = numDofPerNode
    for i = 2:2:nel
      for j = 1:numDofPerNode
        dofs[j,1,i] = m + numDofPerNode*xnodes + j
      end
      for j = 1:numDofPerNode
        dofs[j,2,i] = m + numDofPerNode*(xnodes-1) + j
      end
      for j = 1:numDofPerNode
        dofs[j,3,i] = m + j
      end
      m += numDofPerNode
      for k = 1:nedy
        if i == 2*nedx*k
          m +=numDofPerNode
        end
      end
    end

    return dofs, nnp, nel, nen
  elseif nnpe == 3
    nen = 7 # Number of element nodes
    nnp = xnodes*ynodes + nel;
    dofs = Array{Int}(numDofPerNode, nen, nel)
    m = 0
    for j = 1:nedy
      for i = 1:2:2*nedx
        index = i+2*(j-1)*nedx
        for k = 1:numDofPerNode
          dofs[k,:,index] = [m+k;
                             m+2*numDofPerNode + k;
                             m+2*xnodes*numDofPerNode + k;
                             m+numDofPerNode + k;
                             m+(xnodes+1)*numDofPerNode + k;
                             m+xnodes*numDofPerNode + k;0];
        end
        m = m+2*numDofPerNode;
        if index == 2*nedx*j - 1
          m = m + (xnodes+1)*numDofPerNode;
        end
      end
    end
    m = 2*numDofPerNode;
    for j = 1:nedy
      for i = 2:2:2*nedx
        index = i + 2*(j-1)*nedx;
        for k = 1:numDofPerNode
          dofs[k,:,index] = [m+2*xnodes*numDofPerNode+k;
                             m+(2*xnodes-2)*numDofPerNode+k;
                             m+k;
                             m+((2*xnodes)-1)*numDofPerNode+k;
                             m+(xnodes-1)*numDofPerNode+k;
                             m+xnodes*numDofPerNode+k;0];
        end
        m = m+2*numDofPerNode;
        if index == 2*nedx*j
          m = m + (xnodes+1)*numDofPerNode;
        end
      end
    end
    m = xnodes*ynodes*numDofPerNode;
    for j = 1:nedy
      for i = 1:2*nedx
        index = i + 2*(j-1)*nedx;
        for k = 1:numDofPerNode
          dofs[k,7,index] = m + k;
        end
        m = m+numDofPerNode;
      end
    end

    return dofs, nnp, nel, nen
  elseif nnpe == 4
    nen = 12 # Number of element nodes
    nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 3*nel;
    dofs = Array{Int}(numDofPerNode, nen, nel)
    m = 0
    nmen = 2*nedx +1; # Number of midpoint element nodes in a row
    for j = 1:nedy
      k = 2;
      for i = 1:2:2*nedx
        index = i+2*(j-1)*nedx;
        intv = (j-1)*(xnodes+2*nmen);     # Intermediate variable
        for l = 1:numDofPerNode
          dofs[l,:,index] = [m+l;
                             m+(nnpe-1)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen)*numDofPerNode+l;
                             m+numDofPerNode+l;
                             m+2*numDofPerNode+l;
                             (xnodes+k+intv-1)*numDofPerNode+l;
                             (xnodes+nmen+k+intv-1)*numDofPerNode+l;
                             (xnodes+(k-1)+intv-1)*numDofPerNode+l;
                             (xnodes+nmen+(k-1)+intv-1)*numDofPerNode+l;0;0;0];
        end
        m = m+(nnpe-1)*numDofPerNode;
        k = k+2;
        if index == 2*nedx*j - 1
          m = m+((nnpe-2)*(2*nedx+1)+1)*numDofPerNode;
        end
      end
    end
    m = 3*numDofPerNode;
    for j = 1:nedy
      k = 2;
      for i = 2:2:2*nedx
        index = i + 2*(j-1)*nedx;
        intv = (j-1)*(xnodes+2*nmen);     # Intermediate variable
        for l = 1:numDofPerNode
          dofs[l,:,index] = [m+(xnodes+(nnpe-2)*nmen)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen-3)*numDofPerNode+l;
                             m+l;
                             m+(xnodes+(nnpe-2)*nmen-1)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen-2)*numDofPerNode+l;
                             (xnodes+nmen+k+intv-1)*numDofPerNode+l;
                             (xnodes+k+intv-1)*numDofPerNode+l;
                             (xnodes+nmen+k+intv)*numDofPerNode+l;
                             (xnodes+k+intv)*numDofPerNode+l;0;0;0];
        end
        m = m+(nnpe-1)*numDofPerNode;
        k = k+2;
        if index == 2*nedx*j
          m = m+((nnpe-2)*(2*nedx+1)+1)*numDofPerNode;
        end
      end
    end

    m = (xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1))*numDofPerNode;
    for j = 1:nedy
      for i = 1:2*nedx
        index = i + 2*(j-1)*nedx;
        for l = 1:numDofPerNode
          dofs[l,10,index] = m+l;
          dofs[l,11,index] = m+numDofPerNode+l;
          dofs[l,12,index] = m+2*numDofPerNode+l;
        end          
        m =  m+3*numDofPerNode; 
      end
    end
  
  return dofs, nnp, nel, nen
 elseif nnpe == 5
    nen = 18 # Number of element nodes
    nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 6*nel;
    dofs = Array{Int}(numDofPerNode, nen, nel)
    m = 0
    nmen = 2*nedx +1; # Number of midpoint element nodes in a row
    for j = 1:nedy
      k = 2;
      for i = 1:2:2*nedx
        index = i + 2*(j-1)*nedx;
        intv = (j-1)*(xnodes+(nnpe-2)*nmen);     # Intermediate variable
        for l = 1:numDofPerNode
          dofs[l,:,index] = [m+l;
                             m+(nnpe-1)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen)*numDofPerNode+l;
                             m+numDofPerNode+l;
                             m+2*numDofPerNode+l;
                             m+3*numDofPerNode+l;
                             (xnodes+k+intv-1)*numDofPerNode+l;
                             (xnodes+nmen+k+intv-1)*numDofPerNode+l;
                             (xnodes+(2*nmen)+k+intv-1)*numDofPerNode+l;
                             (xnodes+k-2+intv)*numDofPerNode+l;
                             (xnodes+nmen+k-2+intv)*numDofPerNode+l;
                             (xnodes+(2*nmen)+k-2+intv)*numDofPerNode+l;
                             0;0;0;0;0;0];
        end
        m = m+(nnpe-1)*numDofPerNode;
        k = k+2;
        if index == 2*nedx*j - 1
          m = m+((nnpe-2)*(2*nedx+1)+1)*numDofPerNode;
        end
      end
    end
    m = 4*numDofPerNode;
    for j = 1:nedy
      k = 2;
      for i = 2:2:2*nedx
        index = i + 2*(j-1)*nedx;
        intv = (j-1)*(xnodes+(nnpe-2)*nmen);     # Intermediate variable
        for l = 1:numDofPerNode
          dofs[l,:,index] = [m+(xnodes+(nnpe-2)*nmen)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen-4)*numDofPerNode+l;
                             m+l;
                             m+(xnodes+(nnpe-2)*nmen-1)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen-2)*numDofPerNode+l;
                             m+(xnodes+(nnpe-2)*nmen-3)*numDofPerNode+l;
                             (xnodes+(2*nmen)+k+intv-1)*numDofPerNode+l;
                             (xnodes+nmen+k+intv-1)*numDofPerNode+l;
                             (xnodes+k+intv-1)*numDofPerNode+l;
                             (xnodes+(2*nmen)+k+intv)*numDofPerNode+l;
                             (xnodes+nmen+k+intv)*numDofPerNode+l;
                             (xnodes+k+intv)*numDofPerNode+l;
                             0;0;0;0;0;0];
        end
        m = m+(nnpe-1)*numDofPerNode;
        k = k+2;
        if index == 2*nedx*j
          m = m+((nnpe-2)*(2*nedx+1)+1)*numDofPerNode;
        end
      end
    end 

    m = (xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1))*numDofPerNode;
    for j = 1:nedy
      for i = 1:2*nedx
        index = i + 2*(j-1)*nedx;
        for l = 1:numDofPerNode
          dofs[l,13,index] = m+l;
          dofs[l,14,index] = m+2*numDofPerNode+l;
          dofs[l,15,index] = m+4*numDofPerNode+l;
          dofs[l,16,index] = m+3*numDofPerNode+l;
          dofs[l,17,index] = m+5*numDofPerNode+l;
          dofs[l,18,index] = m+numDofPerNode+l;
        end
        m = m+6*numDofPerNode;
      end
    end
    return dofs, nnp, nel, nen
  end # ends the if-else statement
end # Ends the function

@doc """
### vtxLocation

It creates 3D array storing the X,Y & Z coordinates of the vertices of every element in addition to calculating the element edge length along the X and Y direction 

**Inputs**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nel`     : Number of elements in the mesh

**Outputs**

*  `vtxCoord` : 3D array storing the X,Y & Z coordinates of the vertices. First index stores the coordinates, 2nd index has the local vertex numbering, 3rd index is the element number
*  `elem_edge_lengthx : Length of element edge along the X-axis
*  `elem_edge_lengthy : Length of element edge along the Y-axis

"""->

function vtxLocation(lengthx, lengthy, nedx, nedy, nel)
  elemEdgeLengthx = lengthx/nedx;
  elemEdgeLengthy = lengthy/nedy;
  vtxCoord = zeros(Float64, 3, 3, nel) # Dim1: xyz coord, Dim2: Local vtx numbering, Dim3 Element no
  xvtx = nedx + 1 # Total number of nodes on all adjacent element edges along X-axis 
  yvtx = nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
  m = 0
  n = 0
  for j = 1:nedy
    for i = 1:2:2*nedx
      index = i + 2*(j-1)*nedx;
      vtxCoord[1,1,index] = m*elemEdgeLengthx
      vtxCoord[2,1,index] = n*elemEdgeLengthy
      vtxCoord[1,2,index] = (m+1)*elemEdgeLengthx
      vtxCoord[2,2,index] = n*elemEdgeLengthy
      vtxCoord[1,3,index] = m*elemEdgeLengthx
      vtxCoord[2,3,index] = (n+1)*elemEdgeLengthy
      m += 1
      if index == 2*nedx*j - 1
        m = 0
        n += 1
      end
    end
  end
  m = 1
  n = 1
  for j = 1:nedy
    for i = 2:2:2*nedx
      index = i + 2*(j-1)*nedx;
      vtxCoord[1,1,index] = m*elemEdgeLengthx
      vtxCoord[2,1,index] = n*elemEdgeLengthy
      vtxCoord[1,2,index] = (m-1)*elemEdgeLengthx
      vtxCoord[2,2,index] = n*elemEdgeLengthy
      vtxCoord[1,3,index] = m*elemEdgeLengthx
      vtxCoord[2,3,index] = (n-1)*elemEdgeLengthy
      m += 1
      if index == 2*nedx*j
        m = 1
        n += 1
      end
    end
  end

  return vtxCoord, elemEdgeLengthx, elemEdgeLengthy
end # Ends the function

@doc """
### boundaryNodeInfo

It generates information about nodes on the boundary for implementing boundary conditions.

**Inputs**

*  `nedx`   : Number of element edges along X-axis
*  `nedy`   : Number of element edges along Y-axis
*  `nnpe`   : Number of nodes per element edge (same for all edges)

**Outputs**

*  `NodeEdgex` : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  `NodeEdgey` : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.

"""->

function boundaryNodeInfo(nedx,nedy,nnpe)   
  if nnpe == 2
    xnodes::Int = nedx + 1 # Total number of nodes on all adjacent element edges along X-axis
    ynodes::Int = nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
    NodeEdgex = zeros(2,xnodes)
    NodeEdgey = zeros(2,ynodes)
    for i = 1:xnodes
      NodeEdgex[1,i] = i
      NodeEdgex[2,i] = xnodes*nedy + i
    end
    for j = 1:ynodes
      NodeEdgey[1,j] = j*xnodes
      NodeEdgey[2,j] = 1+(j-1)*xnodes
    end
    return NodeEdgex, NodeEdgey

  elseif nnpe == 3
    xnodes = (nedx*2)+1; # Total number of nodes on all adjacent element edges along X-axis
    ynodes = (nedy*2)+1; # Total number of nodes on all adjacent element edges along Y-axis
    NodeEdgex = zeros(2,xnodes)
    NodeEdgey = zeros(2,ynodes)
    for i = 1:xnodes
      NodeEdgex[1,i] = i;
      NodeEdgex[2,i] = xnodes*nedy + (nnpe-2)*(2*nedx+1)*nedy + i;
    end
    for j = 1:ynodes
      NodeEdgey[1,j] = j*xnodes
      NodeEdgey[2,j] = 1+(j-1)*xnodes
    end
    return NodeEdgex, NodeEdgey

  elseif nnpe == 4
    xnodes = (nnpe-1)*nedx + 1; # Total number of nodes on all adjacent element edges along X-axis
    ynodes = (nnpe-1)*nedy + 1; # Total number of nodes on all adjacent element edges along Y-axis
    NodeEdgex = zeros(2,xnodes)
    NodeEdgey = zeros(2,ynodes)
    intv = 2*nedx + 1 # Intermediate variable which stores the value of number of midpoint nodes in a row
    for i = 1:xnodes
      NodeEdgex[1,i] = i;
      NodeEdgex[2,i] = xnodes*nedy + (nnpe-2)*(intv)*nedy + i;
    end
    m = 1;
    for j = 1:nedy
      NodeEdgey[1,m] = xnodes + (j-1)*(xnodes+(nnpe-2)*(intv));
      NodeEdgey[1,m+1] = (xnodes + intv) + (j-1)*(xnodes+(nnpe-2)*(intv));
      NodeEdgey[1,m+2] = (xnodes + 2*intv) + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m] = 1 + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m+1] = 1 + xnodes + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m+2] = 1 + xnodes + intv  + (j-1)*(xnodes+(nnpe-2)*intv);
      if j == nedy
        NodeEdgey[1,m+3] = 2*(xnodes+intv) + (j-1)*(xnodes+(nnpe-2)*intv);
        NodeEdgey[2,m+3] = 1 + xnodes + 2*intv + (j-1)*(xnodes+(nnpe-2)*intv);
      end
      m += 3;
    end
    return NodeEdgex, NodeEdgey

  elseif nnpe == 5
    xnodes = (nnpe-1)*nedx + 1; # Total number of nodes on all adjacent element edges along X-axis
    ynodes = (nnpe-1)*nedy + 1; # Total number of nodes on all adjacent element edges along Y-axis
    NodeEdgex = zeros(2,xnodes)
    NodeEdgey = zeros(2,ynodes)
    intv = 2*nedx + 1
    for i = 1:xnodes
      NodeEdgex[1,i] = i;
      NodeEdgex[2,i] = xnodes*nedy + (nnpe-2)*(2*nedx+1)*nedy + i;
    end
    m = 1;
    for j = 1:nedy
      NodeEdgey[1,m] = xnodes + (j-1)*(xnodes+(nnpe-2)*(intv));
      NodeEdgey[1,m+1] = (xnodes + intv) + (j-1)*(xnodes+(nnpe-2)*(intv));
      NodeEdgey[1,m+2] = (xnodes + 2*intv) + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[1,m+3] = (xnodes + 3*intv) + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m] = 1 + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m+1] = 1 + xnodes + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m+2] = 1 + xnodes + intv  + (j-1)*(xnodes+(nnpe-2)*intv);
      NodeEdgey[2,m+3] = 1 + xnodes + 2*intv + (j-1)*(xnodes+(nnpe-2)*intv);
      if j == nedy
        NodeEdgey[1,m+4] = 2*xnodes + 3*intv + (j-1)*(xnodes+(nnpe-2)*intv);
        NodeEdgey[2,m+4] = 1 + xnodes + 3*intv + (j-1)*(xnodes+(nnpe-2)*intv);
      end
      m += 4;
    end
    return NodeEdgex, NodeEdgey 
  end # Ends the if statement
end # Ends function boundaryNodeInfo

@doc """
### boundaryEdges

This function identifies the edges and elements on the oundary. It is of the type Boundary

**Inputs**

*  `nedx`   : Number of element edges along X-axis
*  `nedy`   : Number of element edges along Y-axis
*  `nel`    : Number of elements in the mesh

**Outputs**

*  `boundaries` : Array containing the edges and elements on the boundary
*  `nedges`     : Total number of edges in the mesh

"""->

function boundaryEdges(nedx,nedy,nel)
  nedges = nedx*(nedy+1) + nedx*nedy + (nedx+1)*nedy;
  boundaries = Array{Boundary}(2*(nedx+nedy))
  m = 1
  for i = 1:nedx
    boundaries[i] = Boundary(m,1)
    m += 2
  end
  m = nel - 2*nedx+2
  for i = 1:nedx
    boundaries[nedx+i] = Boundary(m,1)
    m += 2
  end
  m = 1
  for i = 1:nedy
    boundaries[2*nedx+i] = Boundary(m,3)
    m += 2*nedx
  end
  m = 2*nedx
  for i = 1:nedy
    boundaries[2*nedx+nedy+i] = Boundary(m,3)
    m += 2*nedx
  end

  # println(boundaries)
  return boundaries, nedges
end

@doc """ 
### interiorEdgeinfo

This function identifies all the element interfaces (edges in 2D and faces in 3D) in the grid.
It is of the type Interface.

**Inputs**

*  `nedx`   : Number of element edges along X-axis
*  `nedy`   : Number of element edges along Y-axis
*  `nel`    : Number of elements in the mesh

**Outputs**

* `interfaces` : Array of type Interface containing the interior edge info

"""->

function interiorEdgeinfo(nedx,nedy,nel)
  nedges = nedx*(nedy+1) + nedx*nedy + (nedx+1)*nedy;
  nIntFace = 0
  interfaces = Array{Interface}(nIntFace)
  for j = 1:nedy
    for i = 1:2:2*nedx
      index = i + 2*(j-1)*nedx;
      if index == 1 + 2*(j-1)*nedx
        if index == 1
          nIntFace += 1
          interfaces = resize!(interfaces,nIntFace)
          interfaces[nIntFace] = Interface(index,index+1,2,2)
        else
          nIntFace += 2
          interfaces = resize!(interfaces,nIntFace)
          interfaces[nIntFace-1] = Interface(index,index-2*nedx+1,1,1)
          interfaces[nIntFace] = Interface(index,index+1,2,2)
        end
      elseif index < 2*nedx
        nIntFace += 2
        interfaces = resize!(interfaces,nIntFace)
        interfaces[nIntFace-1] = Interface(index,index-1,3,3)
        interfaces[nIntFace] = Interface(index,index+1,2,2)
      else
        nIntFace += 3
        interfaces = resize!(interfaces,nIntFace)
        interfaces[nIntFace-2] = Interface(index,index-2*nedx+1,1,1)
        interfaces[nIntFace-1] = Interface(index,index-1,3,3)
        interfaces[nIntFace] = Interface(index,index+1,2,2)
      end
    end
  end

  return interfaces
end

@doc """
### boundaryEdgeInfo

This function identifies the mesh edges on the boundary of the geometry.

**Inputs**

*  `nedx`   : Number of element edges along X-axis
*  `nedy`   : Number of element edges along Y-axis

**Outputs**

*  `HBedges` : A 2*nedx array which stores the edge numbers on the boundary along the X-axis. The first row stores the edges along the base of the rectangle. The second row stores edges on the top of the rectangle.
*  `VBedges` : A 2*nedy array which stores the edge numbers on the boundary along the Y-axis. The first row stores the vertical edges along the right side of the rectangle. The second row stores the edges on the left side of the rectangle.
*  `nedges`  : Total number of mesh edges in the geometry.

"""->

function boundaryEdgeInfo(nedx,nedy)
  # Boundary Edges
  nedges = nedx*(nedy+1) + nedx*nedy + (nedx+1)*nedy;
  HBedges = zeros(2,nedx);  # Horizontal Boundary edges
  VBedges = zeros(2,nedy);  # Vertical Boundary edges
  for i = 1:nedx
    HBedges[1,i] = i;
    HBedges[2,i] = nedx*nedy + i;
  end  
  for i = 1:nedy
    VBedges[1,i] = nedx*(nedy+1) + nedx*nedy + 1 + (i-1)*(nedx+1);
    VBedges[2,i] = nedx*(nedy+1) + nedx*nedy + i*(nedx+1);
  end
  return HBedges, VBedges, nedges
end   # Ends function boundaryEdgeInfo
