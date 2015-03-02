module SimpleMeshType
# Function to create a simple Mesh type.

@doc """
### createMesh

It creates a simple mesh with triangle elements in a rectangular domain in 2D. In This code The horizontal axis is referred to as X-axis and the vertical axis is referred to as the Y-axis.

**Inputs**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis

**Outputs**

*  `xnodes` : Number of nodes along X-axis in one row
*  `ynodes` : Number of nodes along Y-axis in one column
*  `nnp`    : Number of nodal points
*  `nel`    : Number of elements mesh
*  `IEN`    : Element information matrix
*  `vtx_loc`: Location coordinates of the vertices

We assume a counter-clockwise numbering scheme for the the elements and edges with the numbering starting from the base of the geometry.
"""->

function createMesh(lengthx,lengthy,nedx,nedy,nnpe)

if nnpe == 2
	xnodes = nedx + 1
	ynodes = nedy + 1
	nnp = xnodes*ynodes	# number of nodal points
	nel = 2*(xnodes - 1)*(ynodes - 1)
	nsd = 2             # Number of spatial dimensions
	IEN = zeros(3,nel)
	  
	# Populating the IEN Matrix
	m = 1
	for j = 1:nedy
	    for i = 1:2:2*nedx
	      index = i+2*(j-1)*nedx
	      IEN[:,index] = [m;m+1;m+xnodes]
	      m = m+1
	      if index == 2*nedx*j - 1
	          m = m+1
	      end
	    end
	end
	m = 2
	for j = 1:nedy
	    for i = 2:2:2*nedx
	        index = i + 2*(j-1)*nedx
	        IEN[:,index] = [m;m+xnodes;m+(xnodes-1)]
	        m = m+1
	        if index == 2*nedx*j
	            m = m+1
	        end
	    end
	end
	vtx_loc = zeros(2,nnp);
        elem_edge_lengthx = lengthx/nedx;
        elem_edge_lengthy = lengthy/nedy;
        for j = 1:(ynodes)
            for i = 1:(xnodes)
                vtx_loc[1,i+(j-1)*xnodes] = (i-1)*elem_edge_lengthx;
                vtx_loc[2,i+(j-1)*xnodes] = (j-1)*elem_edge_lengthy;
            end
        end
	return xnodes,ynodes,nnp,nel,IEN
	break

elseif nnpe == 3
	xnodes = (nedx*2)+1;
  ynodes = (nedy*2)+1;
  nel = 2*nedx*nedy;
  nnp = xnodes*ynodes + nel;
  IEN = zeros(7,nel);
  m = 1;
  for j = 1:nedy
    for i = 1:2:2*nedx
      index = i+2*(j-1)*nedx;
      IEN[:,index] = [m;m+2;m+(2*xnodes);m+1;m+xnodes+1;m+xnodes;0];
      m = m+2;
      if index == 2*nedx*j - 1
        m = m+xnodes+1;
      end
    end
  end
  m = 3;
  for j = 1:nedy
    for i = 2:2:2*nedx
      index = i + 2*(j-1)*nedx;
      IEN[:,index] = [m;m+2*xnodes;m+(2*xnodes)-2;m+xnodes;m+(2*xnodes)-1;m+xnodes-1;0];
      m = m+2;
      if index == 2*nedx*j
        m = m+xnodes+1;
      end
    end
  end
  m = xnodes*ynodes + 1;
  for j = 1:nedy
    for i = 1:2*nedx
      index = i + 2*(j-1)*nedx;
      IEN[7,index] = m;
      m = m+1;
    end
  end

  # Determining the locations
  elem_edge_lengthx = lengthx/nedx;
  elem_edge_lengthy = lengthy/nedy;
  vtx_loc = zeros(2,nnp);
  m = 0;
  n = 0;
  for j = 1:2:ynodes
    for i = 1:2:xnodes
      vtx_loc[1,i+(j-1)*xnodes] = m*elem_edge_lengthx;
      vtx_loc[2,i+(j-1)*xnodes] = n*elem_edge_lengthy;
      m = m+1;
    end
    m = 0;
    n = n+1;
  end
  break

elseif nnpe==4
	xnodes = (nnpe-1)*nedx + 1;
  ynodes = (nnpe-1)*nedy + 1;
  nel = 2*nedx*nedy;
  nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 3*nel;
  IEN = zeros(12,nel);
  m = 1;
  nmen = 2*nedx +1; # Number of midpoint element nodes in a row
  for j = 1:nedy
    k = 2;
    for i = 1:2:2*nedx
      index = i+2*(j-1)*nedx;
      intv = (j-1)*(xnodes+2*nmen);     # Intermediate variable
      IEN[:,index] = [m;m+(nnpe-1);m+xnodes+(nnpe-2)*nmen;m+1;m+2;xnodes+k+intv;xnodes+nmen+k+intv;xnodes+(k-1)+intv;xnodes+nmen+(k-1)+intv;0;0;0];
      m = m+(nnpe-1);
      k = k+2;
      if index == 2*nedx*j - 1
        m = m+(nnpe-2)*(2*nedx+1)+1;
      end
    end
  end
  m = 4;
  for j = 1:nedy
    k = 2;
    for i = 2:2:2*nedx
      index = i + 2*(j-1)*nedx;
      intv = (j-1)*(xnodes+2*nmen);     # Intermediate variable
      IEN[:,index] = [m;m+xnodes+(nnpe-2)*nmen;m+xnodes+(nnpe-2)*nmen-3;xnodes+(k+1)+intv;xnodes+nmen+(k+1)+intv;m+xnodes+(nnpe-2)*nmen-1;m+xnodes+(nnpe-2)*nmen-2;xnodes+k+intv;xnodes+nmen+k+intv;0;0;0];
      m = m+(nnpe-1);
      k = k+2;
      if index == 2*nedx*j
        m = m+(nnpe-2)*(2*nedx+1)+1;
      end
    end
  end
  m = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 1;
  for j = 1:nedy
    for i = 1:2*nedx
      index = i + 2*(j-1)*nedx;
      IEN[10,index] = m;
      IEN[11,index] = m+1;
      IEN[12,index] = m+2;
      m = m+3;
    end
  end
  # Determining the location
  elem_edge_lengthx = lengthx/nedx;
  elem_edge_lengthy = lengthy/nedy;
  vtx_loc = zeros(2,nnp);
  m = 0;
  n = 0;
  for j = 1:(nedy+1)
    for i = 1:3:xnodes
      vtx_loc[1,i+(j-1)*(xnodes+2*nmen)] = m*elem_edge_lengthx;
      vtx_loc[2,i+(j-1)*(xnodes+2*nmen)] = n*elem_edge_lengthy;
      m = m+1;
    end
    n = n+1;
    m = 0;
  end 
  break

end

@doc """
### boundaryInfo

It generates information about nodes and elements on the boundary for
implementing boundary conditions.

**Inputs**

*  `xnodes` : Number of nodes along X-axis in one row
*  `ynodes` : Number of nodes along Y-axis in one column

**Outputs**
*  `NodeEdgex` : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  `NodeEdgey` : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.
*  `ElemEdgex` : A 2*(xnodes-1) matrix that stores boundary elements along the X-axis. The First row stores elements on the base of the geometry. The second row stores elements on the top of the geometry.
*  `ElemEdgey` : A 2*(ynodes-1) matrix that stores boundary elements along the Y-axis. The first row stores elements on the left side of the geopmetry while the second row stores elements on the right side of the geometry.

"""->

function boundaryInfo(xnodes,ynodes)
# Identifying which nodes are on the boundary
  nedx = xnodes - 1
	nedy = ynodes - 1
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

	# Identifying Elements that have boundary node(s)
	ElemEdgex = zeros(2,2*nedx)
	ElemEdgey = zeros(2,2*nedy)
    
    # Populating ElemEdgex
	for i = 1:(2*nedx)
	  ElemEdgex[1,i] = i
	  ElemEdgex[2,i] = 2*nedx*(ynodes-2) + i
	end

    # Populating ElemEdgey
	m = 0;
	for j = 1:2:2*nedy
	  ElemEdgey[1,j] = 1+2*nedx*m
	  ElemEdgey[1,j+1] = ElemEdgey[1,j]+1
	  ElemEdgey[2,j] = 2*nedx-1 + 2*nedx*m
	  ElemEdgey[2,j+1] = ElemEdgey[2,j] + 1
	  m = m+1
	end
	return NodeEdgex,NodeEdgey,ElemEdgex,ElemEdgey
end
end  # module
