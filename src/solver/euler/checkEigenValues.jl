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


function elementEigenValues{Tmsh,Tsol,Tdim}(mesh::AbstractMesh{Tmsh}, 
                                       sbp::SBPOperator, 
                                       eqn::EulerData{Tsol,Tdim})

  res_0 = zeros(eqn.SL)
  res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)
  pert = 1e-6
  numDofPerElement = mesh.numDofPerNode*mesh.numNodesPerElement
  elem_eigen_vals = zeros(numDofPerElement)
  elem_count = 0 # Count the number of elements with negative eigen values
  #=
  for i = 1:mesh.numEl
    elem_jac = zeros(numDofPerElement,numDofPerElement)
    elemSL = zeros(numDofPerElement) # Element wise residual
    elemSL0 = zeros(elemSL)  # 1D array containing element wise conservative variables
    elemRes0 = zeros(elemSL) # Collects the original residual
    entry_orig = zero(eltype(eqn.SL0)) # stores one element of original SL) at a time
    orig_SL = copy(eqn.SL)

    orig_SL0 = copy(eqn.SL)
    for j = 1:mesh.numNodesPerElement # Populating elemRes0 & elemSL0
      for k = 1:mesh.numDofPerNode
        elemRes0[(j-1)*mesh.numDofPerNode + k] = eqn.res[k,j,i]
        elemSL0[(j-1)*mesh.numDofPerNode + k] = eqn.q[k,j,i]
      end
    end

    # println("elemRes0 - eqn.SL = \n", elemRes0 - eqn.SL)
    # println("elemSL0 - eqn.SL0 = \n", elemSL0 - eqn.SL0)

    for j = 1:numDofPerElement # Perturb all dofs in the element  
      if j == 1
        entry_orig = elemSL0[j]
        elemSL0[j] += pert  # Perturb the jth element of q
      else
        elemSL0[j-1] = entry_orig
        entry_orig = elemSL0[j]
        elemSL0[j] += pert
      end
      #println("elemSL0 - eqn.SL0 = \n", elemSL0 - eqn.SL0)
      # Disassembling SL0 to q at the element level
      for l = 1:mesh.numNodesPerElement
        for k = 1:mesh.numDofPerNode
          eqn.q[k,l,i] = elemSL0[(l-1)*mesh.numDofPerNode + k]
        end
      end
      evalEuler(mesh, sbp, eqn, opts) # Evaluate function with perturbed q
      
      # Populate elemSL with perturbed residuals
      for l = 1:mesh.numNodesPerElement
        for k = 1:mesh.numDofPerNode
          elemSL[(l-1)*mesh.numDofPerNode + k] = eqn.res[k,l,i]
        end
      end
      
      fill!(eqn.SL, 0.0)
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res,  eqn.SL)
      # println("elemSL - eqn.SL = \n", elemSL - eqn.SL)
      # println("orig_SL - eqn.SL = \n", orig_SL - eqn.SL)

      nl_solvers.calcJacRow(unsafe_view(elem_jac, :, j), elemRes0, elemSL, pert)
    end  # End for j = = 1:numDofPerElement
    
    # println("elemSL = \n", elemSL)
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
  res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)
  jac = zeros(Float64, mesh.numDof, mesh.numDof)
  res_0 = copy(eqn.SL)
  nl_solvers.calcJacFD(mesh, sbp, eqn, opts, evalEuler, res_0, pert, jac)
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