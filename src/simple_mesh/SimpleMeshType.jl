module SimpleMeshType
# Creating the types.

export createMesh, nodeLocation, boundaryInfo

abstract SimpleMesh{T <: FloatingPoint}

@doc """
### createMesh

It creates a simple mesh with triangle elements in a rectangular domain in 2D. In This code The horizontal axis is referred to as X-axis and the vertical axis is referred to as the Y-axis.

**Fields**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nnpe`    : Number of nodes per element edge (Same for all edges)
*  `nnp`     : Number of nodal points
*  `nel`     : Number of elements mesh
*  `IEN`     : Element information matrix
*  `vtx_loc` : Location coordinates of the vertices

We assume a counter-clockwise numbering scheme for the the elements and edges with the numbering starting from the base of the geometry.
"""->
type createMesh{T} <: SimpleMesh{T}
	
	lengthx::FloatingPoint
	lengthy::FloatingPoint
	nedx::Int
	nedy::Int
	nnpe::Int
	nnp::Int
	nel::Int
	IEN::Array{Int}
	# vtx_loc::Array{T}
	function createMesh(lengthx,lengthy,nedx,nedy,nnpe)
		if nnpe == 2
			xnodes::Int = nedx + 1 # Total number of nodes on all adjacent element edges along X-axis
			ynodes::Int = nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
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
			#=
			vtx_loc = zeros(2,nnp);
		  elem_edge_lengthx = lengthx/nedx;
		  elem_edge_lengthy = lengthy/nedy;
		  for j = 1:(ynodes)
		    for i = 1:(xnodes)
		      vtx_loc[1,i+(j-1)*xnodes] = (i-1)*elem_edge_lengthx;
		      vtx_loc[2,i+(j-1)*xnodes] = (j-1)*elem_edge_lengthy;
		    end
		  end =#
		  return IEN

		elseif nnpe == 3
			xnodes = (nedx*2)+1; # Total number of nodes on all adjacent element edges along X-axis
		  ynodes = (nedy*2)+1; # Total number of nodes on all adjacent element edges along Y-axis
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
		  #=
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
		  end =#

		elseif nnpe==4
			xnodes = (nnpe-1)*nedx + 1; # Total number of nodes on all adjacent element edges along X-axis
		  ynodes = (nnpe-1)*nedy + 1; # Total number of nodes on all adjacent element edges along Y-axis
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
		  #=
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
		  end =#

		elseif nnpe==5
		  xnodes = (nnpe-1)*nedx + 1; # Total number of nodes on all adjacent element edges along X-axis
		  ynodes = (nnpe-1)*nedy + 1; # Total number of nodes on all adjacent element edges along Y-axis
		  nel = 2*nedx*nedy;
		  nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 6*nel;
		  IEN = zeros(18,nel);
		  m = 1;
		  nmen = 2*nedx +1; # Number of midpoint element nodes in a row
		  for j = 1:nedy
		    k = 2;
		    for i = 1:2:2*nedx
		      index = i + 2*(j-1)*nedx;
		      intv = (j-1)*(xnodes+(nnpe-2)*nmen);     # Intermediate variable
		      IEN[:,index] = [m;m+(nnpe-1);m+xnodes+(nnpe-2)*nmen;m+1;m+2;m+3;
		                      xnodes+k+intv;xnodes+nmen+k+intv;xnodes+(2*nmen)+k+intv;
		                      xnodes+(k-1)+intv;xnodes+nmen+(k-1)+intv;
		                      xnodes+(2*nmen)+(k-1)+intv;0;0;0;0;0;0];
		      m = m+(nnpe-1);
		      k = k+2;
		      if index == 2*nedx*j - 1
		        m = m+(nnpe-2)*(2*nedx+1)+1;
		      end
		    end
		  end
		  m = nnpe;
		  for j = 1:nedy
		    k = 2;
		    for i = 2:2:2*nedx
		      index = i + 2*(j-1)*nedx;
		      intv = (j-1)*(xnodes+(nnpe-2)*nmen);     # Intermediate variable
		      IEN[:,index] = [m;m+xnodes+(nnpe-2)*nmen;m+xnodes+(nnpe-2)*nmen-4;
		                      xnodes+(k+1)+intv;xnodes+nmen+(k+1)+intv;
		                      xnodes+(2*nmen)+(k+1)+intv;m+xnodes+(nnpe-2)*nmen-1;
		                      m+xnodes+(nnpe-2)*nmen-2;m+xnodes+(nnpe-2)*nmen-3;
		                      xnodes+k+intv;xnodes+nmen+k+intv;xnodes+(2*nmen)+k+intv;
		                      0;0;0;0;0;0];
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
		      IEN[13,index] = m;
		      IEN[14,index] = m+1;
		      IEN[15,index] = m+2;
		      IEN[16,index] = m+3;
		      IEN[17,index] = m+4;
		      IEN[18,index] = m+5;
		      m = m+6;
		    end
		  end

		  # Determining the location of the vertices and storing them in an array which
		  # is supposed to contain location of all the nodes 
		  #= 
		  elem_edge_lengthx = lengthx/nedx;
		  elem_edge_lengthy = lengthy/nedy;
		  vtx_loc = zeros(2,nnp);
		  m = 0;
		  n = 0;
		  for j = 1:(nedy+1)
		    for i = 1:4:xnodes
		      vtx_loc[1,i+(j-1)*(xnodes+(nnpe-2)*nmen)] = m*elem_edge_lengthx;
		      vtx_loc[2,i+(j-1)*(xnodes+(nnpe-2)*nmen)] = n*elem_edge_lengthy;
		      m = m+1;
		    end
		    n = n+1;
		    m = 0;
		  end
		end # Ends the if-else statement 
		=#
	end # Ends the function
end # Ends the type 

@doc """
### nodeLocation
It creates 2D array of length equal to the number of nodal points and stores the X & Y coordinates of all the nodes. 

**Fields**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nnpe`    : Number of nodes per element edge (same for all edges)
*  `elem_edge_lengthx`  : Length of element edge along the X-axis
*  `elem_edge_lengthy`  : Length of element edge along the Y-axis
*  `vtx_loc`  : 2D array which stores the X & Y coordinates of all nodes

NOTE: The function "nodeLocation" inside the type only calculates the coordinates of the nodes on the vertices. The other coordinate locations can be obtained by calling SummationByParts 

""" ->
type nodeLocation{T} <: SimpleMesh{T}
	
	lengthx::FloatingPoint
	lengthy::FloatingPoint
	nedx::Int
	nedy::Int
	nnpe::Int
	elem_edge_lengthx::FloatingPoint
	elem_edge_lengthy::FloatingPoint
	vtx_loc::Array{T}
	  
	function nodeLocation(lengthx,lengthy,nedx,nedy,nnpe)
		if nnpe == 2
			xnodes::Int = nedx + 1 # Total number of nodes on all adjacent element edges along X-axis 
			ynodes::Int = nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
			nnp = xnodes*ynodes	# number of nodal points
			nel = 2*(xnodes - 1)*(ynodes - 1)
			nsd = 2             # Number of spatial dimensions
		  elem_edge_lengthx = lengthx/nedx;
		  elem_edge_lengthy = lengthy/nedy;
		  vtx_loc = zeros(2,nnp);
		  for j = 1:(ynodes)
		    for i = 1:(xnodes)
		      vtx_loc[1,i+(j-1)*xnodes] = (i-1)*elem_edge_lengthx;
		      vtx_loc[2,i+(j-1)*xnodes] = (j-1)*elem_edge_lengthy;
		    end
		  end

		elseif nnpe == 3
			xnodes = (nedx*2)+1;
		  ynodes = (nedy*2)+1;
		  nel = 2*nedx*nedy;
		  nnp = xnodes*ynodes + nel;
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

		elseif nnpe==4
			xnodes = (nnpe-1)*nedx + 1;
		  ynodes = (nnpe-1)*nedy + 1;
		  nel = 2*nedx*nedy;
		  nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 3*nel;
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

		elseif nnpe==5
		  xnodes = (nnpe-1)*nedx + 1;
		  ynodes = (nnpe-1)*nedy + 1;
		  nel = 2*nedx*nedy;
		  nnp = xnodes*(nedy+1) + (nnpe-2)*nedy*(2*nedx+1) + 6*nel;
		  elem_edge_lengthx = lengthx/nedx;
		  elem_edge_lengthy = lengthy/nedy;
		  vtx_loc = zeros(2,nnp);
		  m = 0;
		  n = 0;
		  for j = 1:(nedy+1)
		    for i = 1:4:xnodes
		      vtx_loc[1,i+(j-1)*(xnodes+(nnpe-2)*nmen)] = m*elem_edge_lengthx;
		      vtx_loc[2,i+(j-1)*(xnodes+(nnpe-2)*nmen)] = n*elem_edge_lengthy;
		      m = m+1;
		    end
		    n = n+1;
		    m = 0;
		  end
		end # Ends if-else statement
	end # Ends the function
end # Ends the type

@doc """
### boundaryInfo

It generates information about nodes and elements on the boundary for implementing boundary conditions.

**Fields**

*  `nedx`   : Number of element edges along X-axis
*  `nedy`   : Number of element edges along Y-axis
*  `nnpe`   : Number of nodes per element edge (This value is same for all edges of an element)
*  `NodeEdgex` : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  `NodeEdgey` : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.
*  `HBedges` : A 2*nedx array which stores the edge numbers on the boundary along the X-axis. The first row stores the edges along the base of the rectangle. The second row stores edges on the top of the rectangle.
*  `VBedges` : A 2*nedy array which stores the edge numbers on the boundary along the Y-axis. The first row stores the vertical edges along the right side of the rectangle. The second row stores the edges on the left side of the rectangle.

The numbering of the edges is such that the horizontal edges are numbered first followed by the diagonal edges and the vertical edges respectively. The numbering is also done from left to right and bottom to top with the starting point always at the lower left corner.

"""->

type boundaryInfo{T} <: SimpleMesh{T}
	nedx::Int
	nedy::Int
	NodeEdgex::Array{Int}
	NodeEdgey::Array{Int}
	nedges::Int
	HBedges::Array{Int}
	VBedges::Array{Int}

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
		end # Ends the if statement
	end # Ends function boundaryNodeInfo

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

	end   # Ends function boundaryEdgeInfo
end # Ends the type

end # Ends Abstract type
end # Ends module