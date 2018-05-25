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
  global const TAG_ADJOINT = "tag_adjoint"
  global const TAG_VISCOUS = "tag_viscous"
  global const TAG_CN = "tag_cn"

  # euler
  global const TAG_ENTROPYVARS = "tag_entropyvars"
  global const TAG_NLSOLVERS = "tag_nlsolvers"
  global const TAG_ESS = "tag_ess"
  global const TAG_UTILS = "tag_utils"
  global const TAG_INPUT = "tag_input"
  global const TAG_REVERSEMODE = "tag_reversemode"
  global const TAG_MISC = "tag_misc"
  global const TAG_CURVILINEAR = "tag_curvilinear"
  global const TAG_HOMOTOPY = "tag_homotopy"
  global const TAG_STAGGERED = "tag_staggered"
  global const TAG_PARALLEL_DERIVATIVES = "tag_parallel_derivatives"
  global const TAG_CHECKPOINT = "tag_checkpoint"
  global const TAG_JAC = "tag_jac"

  global const TAG_PARALLEL = "tag_parallel"

  # types of tests
  global const TAG_SHORTTEST = "tag_shorttest"
  global const TAG_LONGTEST = "tag_longtest"
  global const LengthTags = Set([TAG_SHORTTEST, TAG_LONGTEST])

  global const TAG_TMP = "tag_tmp"  # a temporary tag
end
