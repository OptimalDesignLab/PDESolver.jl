# list of all tags used by all physics modules

# header guard
if !isdefined(:tag_defs)
  tag_defs = false
end

if !tag_defs
  tag_defs = true
  # advection
  global const TAG_COMPLEX = "tag_complex"
  global const TAG_BC = "tag_bc"
  global const TAG_FLUX = "tag_flux"
  global const TAG_VOLUMEINTEGRALS = "tag_volumeintegral"
  global const TAG_CONVERGENCE = "tag_convergence"
  global const TAG_MMS = "tag_mms"
  global const TAG_FRONTEND = "tag_frontend"

  # euler
  global const TAG_ENTROPYVARS = "tag_entropyvars"
  global const TAG_NLSOLVERS = "tag_nlsolvers"
  global const TAG_ESS = "tag_ess"
  global const TAG_UTILS = "tag_utils"
  global const TAG_INPUT = "tag_input"
  global const TAG_MISC = "tag_misc"
  global const TAG_CURVILINEAR = "tag_curvilinear"

end

