module MeshType
# Function to create a simple Mesh type.

@doc """
### MeshType

It creates a simple mesh with triangle elements in a rectangular domain in 2D. In This code The horizontal axis is referred to as X-axis and the vertical axis is referred to as the Y-axis.

**Inputs**

*  'x' : Length on the X-axis
*  'y' : Length on the Y-axis
*  'elemSize' : length of the shortest edge of the element

**Intermediate Variables**

*  'xnodes' : Number of nodes along X-axis in one row
*  'ynodes' : Number of nodes along Y-axis in one column
*  'nnp'    : Number of nodal points
*  'nel'    : Number of elements mesh
*  'nedx'   : Number of element edges along x
*  'nedy'   : Number of element edges along y 
*  'nsd'    : Number of spatial dimensions
*  'IEN'    : Element information matrix
*  'NodeEdgex' : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  'NodeEdgey' : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.
*  'ElemEdgex' : A 2*(xnodes-1) matrix that stores boundary elements along the X-axis. The First row stores elements on the base of the geometry. The second row stores elements on the top of the geometry.
*  'ElemEdgey' : A 2*(ynodes-1) matrix that stores boundary elements along the Y-axis. The first row stores elements on the left side of the geopmetry while the second row stores elements on the right side of the geometry.

We assume a counter-clockwise numbering scheme for the the elements and edges with the numbering starting from the base of the geometry.
"""->

function createMesh(x,y,elemSize)
	c = (x/elemSize) + 1
	d = (y/elemSize) + 1
	c = round(c)
	d = round(d)
	xnodes = convert(Int,c)
	ynodes = convert(Int,d)
	nnp = xnodes*ynodes	# number of nodal points
	nel = 2*(xnodes - 1)*(ynodes - 1)
	nedx = xnodes - 1
	nedy = ynodes - 1
	nsd = 2
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
	return xnodes,ynodes,nnp,nel,IEN
end

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
	# println(NodeEdgex)
	# println(NodeEdgey)

	# Identifying Elements that have boundary node(s)
	ElemEdgex = zeros(2,2*nedx)
	ElemEdgey = zeros(2,2*nedy)

	for i = 1:(2*nedx)
	  ElemEdgex[1,i] = i
	  ElemEdgex[2,i] = 2*nedx*(ynodes-2) + i
	end
	# println(ElemEdgex)

	m = 0;
	for j = 1:2:2*nedy
	    ElemEdgey[1,j] = 1+2*nedx*m
	    ElemEdgey[1,j+1] = ElemEdgey[1,j]+1
	    ElemEdgey[2,j] = 2*nedx-1 + 2*nedx*m
	    ElemEdgey[2,j+1] = ElemEdgey[2,j] + 1
	    m = m+1
	end
	# println(ElemEdgey)
	return NodeEdgex,NodeEdgey,ElemEdgex,ElemEdgey
end
end  # module