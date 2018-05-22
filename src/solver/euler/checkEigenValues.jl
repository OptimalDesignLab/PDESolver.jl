# Function to check the eigen values

function checkEigenValues(EigenValues, elementNo, elemcount)
  
  counter = 0
  for i = 1:length(EigenValues)
    if isreal(EigenValues[i]) ==  true
      if EigenValues[i] < 0 && abs(EigenValues[i]) > 1e-12
        counter += 1
        elemcount += 1
      end
    else
      # Check if the real and imaginary parts are positive
      println
      if real(EigenValues[i]) < 0 && abs(real(EigenValues[i])) > 1e-12
        counter += 1
        elemcount += 1
      elseif imag(EigenValues[i]) < 0 && abs(imag(EigenValues[i])) > 1e-12
        counter += 1
        elemcount += 1
      end
    end  # End if isreal(EigenValues[i]) ==  true 
  end # End for i = 1:length(EigenValues)
  if counter > 0
    println("Element ", elementNo, " has negative eigen values")
  end
  
  return nothing
end


function elementEigenValues{Tmsh,Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                       sbp::AbstractSBP, 
                                       eqn::EulerData{Tsol, Tres, Tdim})

  res_0 = zeros(eqn.res_vec)
  res_0_norm = physicsRhs(mesh, sbp, eqn, opts, res_0, (evalResidual,))
  pert = 1e-6
  numDofPerElement = mesh.numDofPerNode*mesh.numNodesPerElement
  elem_eigen_vals = zeros(numDofPerElement)
  elem_count = 0 # Count the number of elements with negative eigen values
  #=
  for i = 1:mesh.numEl
    elem_jac = zeros(numDofPerElement,numDofPerElement)
    elem_res_vec = zeros(numDofPerElement) # Element wise residual
    elem_q_vec = zeros(elem_res_vec)  # 1D array containing element wise conservative variables
    elemRes0 = zeros(elem_res_vec) # Collects the original residual
    entry_orig = zero(eltype(eqn.q_vec)) # stores one element of original res_vec) at a time
    orig_res_vec = copy(eqn.res_vec)

    orig_q_vec = copy(eqn.res_vec)
    for j = 1:mesh.numNodesPerElement # Populating elemRes0 & elem_q_vec
      for k = 1:mesh.numDofPerNode
        elemRes0[(j-1)*mesh.numDofPerNode + k] = eqn.res[k,j,i]
        elem_q_vec[(j-1)*mesh.numDofPerNode + k] = eqn.q[k,j,i]
      end
    end

    # println("elemRes0 - eqn.res_vec = \n", elemRes0 - eqn.res_vec)
    # println("elem_q_vec - eqn.q_vec = \n", elem_q_vec - eqn.q_vec)

    for j = 1:numDofPerElement # Perturb all dofs in the element  
      if j == 1
        entry_orig = elem_q_vec[j]
        elem_q_vec[j] += pert  # Perturb the jth element of q
      else
        elem_q_vec[j-1] = entry_orig
        entry_orig = elem_q_vec[j]
        elem_q_vec[j] += pert
      end
      #println("elem_q_vec - eqn.q_vec = \n", elem_q_vec - eqn.q_vec)
      # Disassembling q_vec to q at the element level
      for l = 1:mesh.numNodesPerElement
        for k = 1:mesh.numDofPerNode
          eqn.q[k,l,i] = elem_q_vec[(l-1)*mesh.numDofPerNode + k]
        end
      end
      evalResidual(mesh, sbp, eqn, opts) # Evaluate function with perturbed q
      
      # Populate elem_res_vec with perturbed residuals
      for l = 1:mesh.numNodesPerElement
        for k = 1:mesh.numDofPerNode
          elem_res_vec[(l-1)*mesh.numDofPerNode + k] = eqn.res[k,l,i]
        end
      end
      
      fill!(eqn.res_vec, 0.0)
      array3DTo1D(mesh, sbp, eqn, opts, eqn.res,  eqn.res_vec)
      # println("elem_res_vec - eqn.res_vec = \n", elem_res_vec - eqn.res_vec)
      # println("orig_res_vec - eqn.res_vec = \n", orig_res_vec - eqn.res_vec)

      nl_solvers.calcJacRow(sview(elem_jac, :, j), elemRes0, elem_res_vec, pert)
    end  # End for j = = 1:numDofPerElement
    
    # println("elem_res_vec = \n", elem_res_vec)
    eqn.q[mesh.numDofPerNode,mesh.numNodesPerElement,i] = entry_orig
    
    elem_eigen_vals = eigvals(elem_jac)
    checkEigenValues(elem_eigen_vals, i, elem_count)
    
    #= println("eigen values of element ", i)
    counter = 0
    for kappa = 1:length(elem_eigen_vals)
      if elem_eigen_vals[kappa] < 0 && abs(elem_eigen_vals[kappa]) > 1e-12
        counter += 1
      end
    end =#

    #println(elem_eigen_vals, '\n', '\n')
    #println("number of negative element eigen values = ", counter)
    
  end  

  println("Number of elements with negative eigen values = ", elem_count)
  =#
  res_0_norm = physicsRhs(mesh, sbp, eqn, opts, res_0, (evalResidual,))
  jac = zeros(Float64, mesh.numDof, mesh.numDof)
  res_0 = copy(eqn.res_vec)
  nl_solvers.calcJacFD(mesh, sbp, eqn, opts, evalResidual, res_0, pert, jac)
  println("Jacobian successfully computed") 


  eigen_values = eigvals(jac)
  counter = 0
  for i = 1:length(eigen_values)
    if eigen_values[i] < 0 && abs(eigen_values[i]) > 1e-12
      #println(eigen_values[i])
      counter += 1
    end
  end
  #println(eigen_values)
  println("number of negative eigen values for assembled residual = ", counter)
  println("Number of eigen values of assembled residual = ", length(eigen_values))
  printSolution("eigen_values.dat", eigen_values)  

  return nothing
end
