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

end