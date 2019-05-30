# list of all tags used by all physics modules

# header guard
if !isdefined(:tag_defs)
  tag_defs = false
end

if !tag_defs
  tag_defs = true
  # advection
  global const TAG_COMPLEX = "TAG_COMPLEX"
  global const TAG_BC = "TAG_BC"
  global const TAG_FLUX = "TAG_FLUX"
  global const TAG_VOLUMEINTEGRALS = "TAG_VOLUMEINTEGRAL"
  global const TAG_CONVERGENCE = "TAG_CONVERGENCE"
  global const TAG_MMS = "TAG_MMS"
  global const TAG_FRONTEND = "TAG_FRONTEND"
  global const TAG_ADJOINT = "TAG_ADJOINT"
  global const TAG_VISCOUS = "TAG_VISCOUS"
  global const TAG_CN = "TAG_CN"

  # euler
  global const TAG_ENTROPYVARS = "TAG_ENTROPYVARS"
  global const TAG_NLSOLVERS = "TAG_NLSOLVERS"
  global const TAG_ESS = "TAG_ESS"
  global const TAG_UTILS = "TAG_UTILS"
  global const TAG_INPUT = "TAG_INPUT"
  global const TAG_REVERSEMODE = "TAG_REVERSEMODE"
  global const TAG_MISC = "TAG_MISC"
  global const TAG_CURVILINEAR = "TAG_CURVILINEAR"
  global const TAG_HOMOTOPY = "TAG_HOMOTOPY"
  global const TAG_STAGGERED = "TAG_STAGGERED"
  global const TAG_PARALLEL_DERIVATIVES = "TAG_PARALLEL_DERIVATIVES"
  global const TAG_CHECKPOINT = "TAG_CHECKPOINT"
  global const TAG_JAC = "TAG_JAC"
  global const TAG_FUNCTIONAL = "TAG_FUNCTIONAL"

  global const TAG_PARALLEL = "TAG_PARALLEL"

  # types of tests
  global const TAG_SHORTTEST = "TAG_SHORTTEST"
  global const TAG_LONGTEST = "TAG_LONGTEST"
  global const LengthTags = Set([TAG_SHORTTEST, TAG_LONGTEST])

  global const TAG_TMP = "TAG_TMP"  # a temporary tag
end
