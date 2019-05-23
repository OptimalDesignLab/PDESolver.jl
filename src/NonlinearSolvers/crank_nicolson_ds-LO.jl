
function getCNDSPCandLO(mesh, sbp, eqn, opts)
  
  jac_type = opts["jac_type"]

  if jac_type <= 2
    pc = PCNone(mesh, sbp, eqn, opts)
  elseif opts["use_volume_preconditioner"]
    pc = CNVolumePC(mesh, sbp, eqn, opts)
  else
    pc = CNMatPC(mesh, sbp, eqn, opts)
  end

  if jac_type == 1    # Dense jac
    lo = CNDSDenseLO(pc, mesh, sbp, eqn, otps)
  elseif jac_type == 2    # Sparse jac, julia matrix
    lo = CNDSSparseDirectLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 3
    lo = CNDSPetscMatLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 4
    lo = CNDSPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  end

  return pc, lo
end

#------------------------------------------------------------------------------
# Linear operators

mutable struct CNDSDenseLO <: AbstractDenseLO
  lo_inner::NewtonDenseLO
end

function CNDSDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonDenseLO(pc, mesh, sbp, eqn, opts)

  return CNDSDenseLO(lo_inner)
end

mutable struct CNDSSparseDirectLO <: AbstractSparseDirectLO
  lo_inner::NewtonSparseDirectLO
end

function CNDSSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)

  return CNDSSparseDirectLO(lo_inner)
end


mutable struct CNDSPetscMatLO <: AbstractPetscMatLO
  lo_inner::NewtonPetscMatLO
end

function CNDSPetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = NewtonPetscMatLO(pc, mesh, sbp, eqn, opts)

  return CNDSPetscMatLO(lo_inner)
end
