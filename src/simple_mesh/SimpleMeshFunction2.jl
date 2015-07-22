# SimpleMeshFunction2.jl
# Created for the new simpleMesh type

function createMesh(lengthx,lengthy,nedx,nedy,nnpe,numDofPerNode)
  xnodes = (nnpe-1)*nedx + 1 # Total number of nodes on all adjacent element edges along X-axis
  ynodes = (nnpe-1)*nedy + 1 # Total number of nodes on all adjacent element edges along Y-axis
  nel = 2*nedx*nedy # number of elements
  nsd = 2  # number of spatial dimensions

  if nnpe == 2
  	nen = 3 # Number of element nodes
  	dofs = Array{Int}(numDofPerNode, nen, nel)
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

    println(dofs)



 #=
  elseif nnpe == 3
  	nen = 7 # Number of element nodes
    dofs = Array{Int}(numDofPerNode, nen, nel)

  elseif nnpe == 4
  	nen = 12 # Number of element nodes
  	dofs = Array{Int}(numDofPerNode, nen, nel)

  elseif nnpe == 5
  	nen = 18 # Number of element nodes
  	dofs = Array{Int}(numDofPerNode, nen, nel)
=#
  end
end