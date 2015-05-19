module Mesh2D
# implements a simple mesh for 2D domains with 
# !!!!!!!! This was hacked to make the vortex grid

using SummationByParts

export Mesh, makeperiodic!

type Mesh{T}
  numnodes::Int
  numelem::Int
  Nx::Int
  Ny::Int
  degree::Int
  len::Array{T} # length of sides
  x::Array{T} # node coordinates
  nodeidx::Array{Int} # node indices
  elemnodes::Array{Int} # element-local to global-node mapping
  bndryfaces::Array{Boundary} # boundary faces
  ifaces::Array{Interface} # element interfaces

  function Mesh(Nx::Int, Ny::Int, sbp::TriSBP{T}, lengthx::T, lengthy::T;
                periodic::Bool=false)
    len = zeros(T, 2)
    len[1] = lengthx; len[2] = lengthy
    numnodes, x, nodeidx, elemnodes, bndryfaces, ifaces = 
    Mesh2D.buildrectanglemesh(Nx, Ny, sbp, lengthx, lengthy, periodic=periodic)
    len = [lengthx; lengthy]
    numelem = 2*Nx*Ny
    new(numnodes, numelem, Nx, Ny, sbp.degree, len, x, nodeidx, elemnodes,
        bndryfaces, ifaces)
  end
end

@doc """
Builds a triangle-element mesh for a rectangular domain

**Inputs**

* `Nx`,`Ny`: number of intervals along each side
* `sbp`: triangle SBP operator of order p
* `lengthx`,`lengthy`: lenghts of the rectangular domain
* `periodic`: if True, makes the mesh doubly periodic

**Outputs**

* `numnodes`: number of nodes
* `xy`: node coordinates
* `elemnodes`: local to global mapping of nodes
* `bndryfaces`: boundary faces, if not periodic
* `ifaces`: element interfaces

"""->
function buildrectanglemesh{T}(Nx::Int, Ny::Int, sbp::TriSBP{T}, lengthx::T=1.0,
                               lengthy::T=1.0; periodic::Bool=false)
  @assert(lengthx > 0.0)
  @assert(lengthy > 0.0)
  @assert(Nx >= 1)
  @assert(Ny >= 1)
  #curvy = true  
  #curvy = false
  sector = true
  rin = 1.0
  rout = 3.0

  # get the SBP operator accuracy
  p = sbp.degree

  # find the number of nodes total
  numnodes = (Nx+1)*(Ny+1)
  if p == 2
    # add midedges and centroid
    numnodes += Nx*(Ny+1) + Ny*(Nx+1) + Nx*Ny
    numnodes += 2*Nx*Ny
  elseif p == 3
    # add 2 edge nodes and 3 bubble nodes
    numnodes += 2*(Nx*(Ny+1) + 2*Ny*(Nx+1) + Nx*Ny)
    numnodes += 6*Nx*Ny
  elseif p == 4
    # add 3 edge nodes and 6 bubble nodes
    numnodes += 3*(Nx*(Ny+1) + 2*Ny*(Nx+1) + Nx*Ny)
    numnodes += 12*Nx*Ny
  end

  # set the nodes and the element-node indexing
  xy = zeros(T, (2,numnodes))
  nodeidx = zeros(Int, (numnodes))
  elemnodes = zeros(Int, (sbp.numnodes, 2*Nx*Ny))

  # set the vertices
  ptr = 0
  idxptr = 0
  dx = lengthx/Nx
  dy = lengthy/Ny
  if sector
    function curvygrid(r,θ)
      return [r*cos(θ); r*sin(θ)]
    end
    for j = 1:(Ny+1)
      θ = (j-1)*dy
      for i = 1:(Nx+1)
        r = (i-1)*dx
        xy[:,(j-1)*(Nx+1)+i] = curvygrid(rin + (rout-rin)*r, 0.5*θ*pi)
        idxptr += 1
        nodeidx[(j-1)*(Nx+1)+i] = idxptr
      end
      if periodic # correct indices on right edge
        idxptr -= 1
        nodeidx[(j-1)*(Nx+1)+(Nx+1)] = nodeidx[(j-1)*(Nx+1)+1]
      end
    end
  else # not curvy
    for j = 1:(Ny+1)
      y = (j-1)*dy
      for i = 1:(Nx+1)
        x = (i-1)*dx
        xy[:,(j-1)*(Nx+1)+i] = [x; y]
        idxptr += 1
        nodeidx[(j-1)*(Nx+1)+i] = idxptr
      end
      if periodic # correct indices on right edge
        idxptr -= 1
        nodeidx[(j-1)*(Nx+1)+(Nx+1)] = nodeidx[(j-1)*(Nx+1)+1]
      end
    end
  end
  if periodic # correct indices on the upper edge
    for i = 1:Nx+1
      nodeidx[Ny*(Nx+1)+i] = nodeidx[i]
    end
    idxptr -= Nx
  end

  # set the vertices' local-to-global indexing
  for j = 1:Ny
    for i = 1:Nx
      eptr = (j-1)*Nx + i
      elemnodes[1:3, eptr] = [(j-1)*(Nx+1) + i;
                              (j-1)*(Nx+1) + i + 1;
                              j*(Nx+1)+i];
      eptr += Nx*Ny
      elemnodes[1:3, eptr] = [j*(Nx+1) + i + 1;
                              j*(Nx+1) + i;
                              (j-1)*(Nx+1) + i + 1];
    end
  end
  ptr += (Nx+1)*(Ny+1)
  
  # create edges
  for j = 1:Ny
    for i = 1:Nx
      eptr = (j-1)*Nx + i
      vtx = xy[:,elemnodes[1:3, eptr]].'
      xref = SummationByParts.calcnodes(sbp, vtx)

      # horizontal edge
      xy[:,ptr+1:ptr+p-1] = xref[:,3+1:3+p-1]
      nodeidx[ptr+1:ptr+p-1] = [idxptr+1:idxptr+p-1;]
      # add the nodes to the element
      elemnodes[3+1:3+p-1, eptr] = [ptr+1:ptr+p-1;]
      if j > 1 # this edge is shared with a flipped element
        elemnodes[3+1:3+p-1, Nx*(Ny-1)+eptr] = [ptr+p-1:-1:ptr+1;] # reverse
      end
      ptr += p-1
      idxptr += p-1

      # diagonal edge
      xy[:,ptr+1:ptr+p-1] = xref[:,3+(p-1)+1:3+2*(p-1)]
      nodeidx[ptr+1:ptr+p-1] = [idxptr+1:idxptr+p-1;]
      # add the nodes to the elements
      elemnodes[3+(p-1)+1:3+2*(p-1), eptr] = [ptr+1:ptr+p-1;]
      elemnodes[3+(p-1)+1:3+2*(p-1), Nx*Ny+eptr] = [ptr+p-1:-1:ptr+1;] # reverse
      ptr += p-1
      idxptr += p-1

      # vertical edge
      xy[:,ptr+1:ptr+p-1] = xref[:,3+2*(p-1)+1:3+3*(p-1)]
      nodeidx[ptr+1:ptr+p-1] = [idxptr+1:idxptr+p-1;]
      # add the nodes to the element
      elemnodes[3+2*(p-1)+1:3+3*(p-1), eptr] = [ptr+1:ptr+p-1;]
      if i > 1 # this edge is shared with a flipped element
        elemnodes[3+2*(p-1)+1:3+3*(p-1), Nx*Ny+eptr-1] = [ptr+p-1:-1:ptr+1;] # rev
      end
      ptr += p-1
      idxptr += p-1
    end
    # add vertical edge at the right side of domain
    eptr = Nx*Ny + j*Nx
    vtx = xy[:,elemnodes[1:3, eptr]].'
    xref = SummationByParts.calcnodes(sbp, vtx)
    xy[:,ptr+1:ptr+p-1] = xref[:,3+2*(p-1)+1:3+3*(p-1)]
    if periodic # match with nodes at left side of domain
      e2ptr = (j-1)*Nx + 1
      nodeidx[ptr+1:ptr+p-1] = nodeidx[elemnodes[3+3*(p-1):-1:3+2*(p-1)+1,e2ptr]]
    else # make new node indices
      nodeidx[ptr+1:ptr+p-1] = [idxptr+1:idxptr+p-1;]
      idxptr += p-1
    end
    # add the nodes to the element
    elemnodes[3+2*(p-1)+1:3+3*(p-1), eptr] = [ptr+1:ptr+p-1;]
    ptr += p-1
  end
  # add horizontal edges at the top of domain
  for i = 1:Nx
    eptr = 2*Nx*Ny - Nx + i
    vtx = xy[:,elemnodes[1:3, eptr]].'
    xref = SummationByParts.calcnodes(sbp, vtx)
    xy[:,ptr+1:ptr+p-1] = xref[:,3+1:3+(p-1)]
    if periodic # match with nodes at bottom of domain
      e2ptr = i
      nodeidx[ptr+1:ptr+p-1] = nodeidx[elemnodes[3+(p-1):-1:3+1,e2ptr]]
    else # make new node indices
      nodeidx[ptr+1:ptr+p-1] = [idxptr+1:idxptr+p-1;]
      idxptr += p-1
    end
    # add the nodes to the element
    elemnodes[3+1:3+(p-1), eptr] = [ptr+1:ptr+p-1;]
    ptr += p-1
  end

  # add the interior nodes
  if p > 1
    numintr = sbp.numnodes - sbp.numbndry
    for j = 1:Ny
      for i = 1:Nx
        # add interior nodes for unflipped element
        eptr = (j-1)*Nx + i
        vtx = xy[:,elemnodes[1:3, eptr]].'
        xref = SummationByParts.calcnodes(sbp, vtx)        
        xy[:,ptr+1:ptr+numintr] = xref[:,sbp.numbndry+1:end]
        nodeidx[ptr+1:ptr+numintr] = [idxptr+1:idxptr+numintr;]
        # add the nodes to the element
        elemnodes[sbp.numbndry+1:end, eptr] = [ptr+1:ptr+numintr;]
        ptr += numintr
        idxptr += numintr

        # add interior nodes for flipped element
        eptr += Nx*Ny
        vtx = xy[:,elemnodes[1:3, eptr]].'
        xref = SummationByParts.calcnodes(sbp, vtx)        
        xy[:,ptr+1:ptr+numintr] = xref[:,sbp.numbndry+1:end]
        nodeidx[ptr+1:ptr+numintr] = [idxptr+1:idxptr+numintr;]
        # add the nodes to the element
        elemnodes[sbp.numbndry+1:end, eptr] = [ptr+1:ptr+numintr;]
        ptr += numintr
        idxptr += numintr
      end
    end
  end
  
  # identify the boundary faces
  bndryfaces = Array(Boundary, 2*(Nx + Ny) )
  # boundary faces along the lower and upper edges
  ptr = 1
  for i = 1:Nx
    bndryfaces[ptr] = Boundary(i,1)
    ptr += 1
    bndryfaces[ptr] = Boundary(2*Nx*Ny - Nx + i, 1)
    ptr += 1
  end
  # boundary faces along the left and right edges
  for j = 1:Ny
    bndryfaces[ptr] = Boundary((j-1)*Nx + 1, 3)
    ptr += 1
    bndryfaces[ptr] = Boundary(Nx*Ny + j*Nx, 3)
    ptr += 1
  end
  
  # identify the interfaces
  if periodic 
    ifaces = Array(Interface, 3*Nx*Ny)
  else
    ifaces = Array(Interface, 3*Nx*Ny - Nx - Ny)
  end
  # add interfaces on diagonals
  ptr = 1
  for j = 1:Ny
    for i = 1:Nx
      eptr1 = (j-1)*Nx + i
      eptr2 = eptr1 + Nx*Ny
      ifaces[ptr] = Interface(eptr1, eptr2, 2, 2)
      ptr += 1
    end
  end
  # add horizontal interfaces   
  for j = 2:Ny
    for i = 1:Nx
      eptr1 = (j-1)*Nx + i
      eptr2 = eptr1 - Nx + Nx*Ny
      ifaces[ptr] = Interface(eptr1, eptr2, 1, 1)
      ptr += 1
    end
  end
  # add vertical interfaces
  for j = 1:Ny
    for i = 2:Nx
      eptr1 = (j-1)*Nx + i
      eptr2 = eptr1 - i + Nx*Ny
      ifaces[ptr] = Interface(eptr1, eptr2, 3, 3)
      ptr += 1
    end
  end
  if periodic
    # add horizontal periodic interfaces
    for i = 1:Nx
      eptr1 = i
      eptr2 = Nx*(Ny-1) + i + Nx*Ny
      ifaces[ptr] = Interface(eptr1, eptr2, 1, 1)
      ptr += 1
    end
    # add vertical periodic interfaces
    for j = 1:Ny
      eptr1 = (j-1)*Nx + 1
      eptr2 = j*Nx + Nx*Ny
      ifaces[ptr] = Interface(eptr1, eptr2, 3, 3)
      ptr += 1
    end
  end

  return numnodes, xy, nodeidx, elemnodes, bndryfaces, ifaces
end

@doc """
Writes the mesh to a file.  Triangles are written, as well as the individual
nodes.

**Inputs**

* `mesh`: the mesh object
* `filename`: the name of the file

"""->
function dumpmesh{T}(mesh::Mesh{T}, filename::String)
  fout = open(filename, "w")
  # header information
  println(fout, mesh.numnodes)
  println(fout, mesh.numelem)
  # output the nodes' x,y coordinates
  for i = 1:mesh.numnodes
    print(fout, mesh.x[1,i]," ")
  end
  print(fout, "\n")
  for i = 1:mesh.numnodes
    print(fout, mesh.x[2,i]," ")
  end
  print(fout, "\n")
  # output the triangle nodes indices
  for j = 1:3
    for k = 1:mesh.numelem
      print(fout, mesh.elemnodes[j,k]," ")
    end
    print(fout, "\n")
  end
  close(fout)
end

@doc """

Takes the given mesh and matches nodes on the left-right and top-bottom domains
to make the mesh doubly periodic. *Becareful*: after this is called, the x,y
coordinates are no longer useful.

**InOuts**

* `mesh`: A predefined 2D mesh

"""->
function makeperiodic!{T}(mesh::Mesh{T})
  p = mesh.degree
  # loop over the elmeents, and remap the indices
  for k = 1:size(mesh.elemnodes,2)
    for j = 1:size(mesh.elemnodes,1)
      mesh.elemnodes[j,k] = mesh.nodeidx[mesh.elemnodes[j,k]]
    end
    #println("element = ",k,": nodes = ",mesh.elemnodes[:,k])
  end
  mesh.numnodes -= (mesh.Nx+mesh.Ny+1 + mesh.Nx*(p-1) + mesh.Ny*(p-1))
end

end