export eqn_deepcopy

# One reason for doing this is this case:
#   a = rand(2,2)
#   b = a
#   a[3] = 8
#   b[3] == 8
#   this is because 'a[3] =' is actually setindex!
"""
  EulerEquationMod.eqn_deepcopy

  This function performs a proper deepcopy (unlike julia's builtin deepcopy) 
    on an Euler equation object.
  It preserves reference topology (i.e. q & q_vec pointing to same array in DG schemes).

    Inputs:
      eqn
      mesh
      sbp
      opts

    Outputs:
      eqn_copy

"""
# TODO: check q/q_vec, res/res_vec
# TODO: tests
#  1. ensure copy matches
#  2. ensure changes to eqn don't affect eqn_copy
#  3. ensure changes to eqn_copy.q change eqn_copy.q_vec, same for res

function eqn_deepcopy{Tmsh, Tsol, Tres, Tdim}(eqn::EulerData_{Tsol, Tres, Tdim}, mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, opts::Dict)

  # over 100 fields, so it is necessary to write a better approach for copying than explicitly copying every named field

  # 1: call constructor on eqn_copy
  var_type = opts["variable_type"]

  eqn_copy = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)

  # 2: copy over fields

  for fdnm in fieldnames(eqn)    # loop over first level fieldnames in eqn

    fdnm_type = typeof(getfield(eqn, fdnm))    # get the super type of the current field

    # handle params
    if issubtype(fdnm_type, AbstractParamType)       # if this first level fieldname is of type ParamType; ex: eqn.params
      
      # loop over eqn.params, eqn.params_conservative, or eqn.params_entropy
      println(" is a subtype of AbstractParamType, fdnm: ", fdnm)

      for fdnm_lvl2 in fieldnames(getfield(eqn, fdnm))      # loop over 2nd level fieldnames
        if issubtype(fdnm_type, Array)              # if the 2nd level fieldname is of type Array; ex: eqn.params.q_vals

          # this does not work: setfield!(getfield(eqn_copy, a), b , getfield(getfield(eqn, a),b))
          #   must use copy, or else changing eqn's value changes eqn_copy
          setfield!(getfield(eqn_copy, fdnm), fdnm_lvl2, copy(getfield(getfield(eqn, fdnm), fdnm_lvl2)))

          # Note: this is assuming that there are no Arrays of Arrays inside an eqn.params (or params_entropy, etc)

        else            # if the 2nd level fieldname is not of type Array; ex: eqn.params.gamma

          # because copy is not defined for all non-array types, such as functions
          println("fdnm: ", fdnm)
          println("fdnm_lvl2: ", fdnm_lvl2)
          setfield!(getfield(eqn_copy, fdnm), fdnm_lvl2, getfield(getfield(eqn, fdnm), fdnm_lvl2))

        end

      end
      
    # handle arrays
    elseif issubtype(fdnm_type, Array)        # if this first level fieldname is of type Array; ex: eqn.q or eqn.q_face_send

      # handle array of arrays
      if issubtype(eltype(fdnm_type), Array)        # if this is an Array of Arrays; ex: eqn.q_face_send

        # first copy the outer array
        setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))      # copy is required here, as the innermost object is an array

        # then loop over array and copy all elements, which are each an array
        for i = 1:length(getfield(eqn, fdnm))

          # use getindex/setindex! ?
          setindex!(getfield(eqn_copy, fdnm), getindex(getfield(eqn, fdnm), i), i)

        end   # end of loop over elements (each an array) of the 1st level field, which is of type array

      else      # if this a simple Array; ex: eqn.q

        setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))      # copy is required here, as the innermost object is an array

      end

    else        # handle non-arrays

      setfield!(eqn_copy, fdnm, getfield(eqn, fdnm))      # copy is not defined for many of these non-array types: use assignment

    end

  end     # end of loop over first level fieldnames

      
  return eqn_copy

end
