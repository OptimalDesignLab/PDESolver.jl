
abstract ExactSolutionType

export ExactSolutionType, calcErrorL2Norm

type ExactPolynomial <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactPolynomial, 
						  xy::AbstractArray{Tmsh}, 
						  params::ParamType{2},
						  qe::AbstractArray{Tsol, 1})
	sigma = 0.5	
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
	rhoInf = 1.0
	TInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	x = xy[1]
	y = xy[2]
	rho = (x-x*x)*(y-y*y) 
	u   = (x-x*x)*(y-y*y)
	v   = (x-x*x)*(y-y*y)
	T   = (x-x*x)*(y-y*y)
	rho = (sigma*rho + 1.0)*rhoInf
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf
	qe[1] = rho
	qe[2] = rho*u
	qe[3] = rho*v
	qe[4] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
	qe[4] = qe[4]*rho

	return nothing
end


type ExactPolynomial_1 <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactPolynomial_1, 
						  xy::AbstractArray{Tmsh}, 
						  params::ParamType{2},
						  qe::AbstractArray{Tsol, 1})

	sigma = 0.5
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = xy[1]
	y = xy[2]
	rho = (x-x*x)*(y-y*y) 
	u   = 0.0
	v   = 0.0 
	T   = 0.0 
	rho = (sigma*rho + 1.0)*rhoInf
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf
	qe[1] = rho
	qe[2] = rho*u
	qe[3] = rho*v
	qe[4] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
	qe[4] = qe[4]*rho

	return nothing
end


type ExactPolynomial_2 <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactPolynomial_2, 
						  xy::AbstractArray{Tmsh}, 
						  params::ParamType{2},
						  qe::AbstractArray{Tsol, 1})

	sigma = 0.5
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = xy[1]
	y = xy[2]
	rho = 0.0 
	u   = sqrt((x - x*x) * (y - y*y))
	v   = 0.0
	T   = (x - x*x) * (y - y*y)
	rho = (sigma*rho + 1.0)*rhoInf
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf
	qe[1] = rho
	qe[2] = rho*u
	qe[3] = rho*v
	qe[4] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
	qe[4] = qe[4]*rho

	return nothing
end



type ExactTrigonometric <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactTrigonometric, 
						  xy::AbstractArray{Tmsh}, 
						  params::ParamType{2},
						  qe::AbstractArray{Tsol, 1})
	sigma = 0.5
	pi = 3.14159265358979323846264338
	gamma = 1.4
	gamma_1 = gamma - 1.0
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x2 = 2*xy[1]*pi
	y2 = 2*xy[2]*pi
	x4 = 4*xy[1]*pi
	y4 = 4*xy[2]*pi
	sx2 = sin(x2)
	sy2 = sin(y2)
	sx4 = sin(x4)
	sy4 = sin(y4)
	cx2 = cos(x2)
	cx4 = cos(x4)
	cy2 = cos(y2)
	cy4 = cos(y4)
	#
	# Exact solution in form of primitive variables
	#
	rho = 0.25 * sx2 * sy2
	u   = 0.25 * sx4 * sy4
	v   = 0.25 * (cx4  + 1.0) * (cy4 + 1.0)
	T   = 0.25 * (1.0 - cx4) * (1.0 - cy4)
	rho = (sigma*rho + 1.0)*rhoInf 
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf

	qe[1] = rho
	qe[2] = rho*u
	qe[3] = rho*v
	qe[4] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
	qe[4] = qe[4]*rho

	return nothing
end


global const ExactDict = Dict{ASCIIString, ExactSolutionType}(
	"ExactTrigonometric" => ExactTrigonometric(),
	"ExactPolynomial"    => ExactPolynomial(),
	"ExactPolynomial_1"  => ExactPolynomial_1(),
	"ExactPolynomial_2"  => ExactPolynomial_2(),
)


function calcErrorL2Norm{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
												  sbp::AbstractSBP,
												  eqn::AbstractEulerData{Tsol, Tres},
												  opts)
	l2norm::Float64 = 0.
    lInfnorm::Float64 = 0.
    qe = Array(Tsol, mesh.numDofPerNode)
    exactFunc = ExactDict[opts["exactSolution"]]

	lInfnorm = -1.0
	elem_with_max_dq = 0
	numDofs_to_sum = mesh.numDofPerNode
    for e = 1:mesh.numEl
        for n = 1:mesh.numNodesPerElement
            xy = sview(mesh.coords, :, n, e)
            exactFunc(xy, eqn.params, qe)
            q = sview(eqn.q, :, n, e)
			jac = mesh.jac[n, e]

            for v = 1 : numDofs_to_sum
                # dq = Real(q[v] - qe[v])
                dq = Float64(q[v] - qe[v])
				l2norm += dq*dq*sbp.w[n]/jac

				if lInfnorm < abs(dq)
					lInfnorm = abs(dq)
				end

				# if abs(dq) > 1.0e-4
					# println(e, ", ", n, ", ", v, ", ", real(dq), ", ", real(q[v]), ", ", real(qe[v]))
				# end
			end
		end
    end
	
    return sqrt(l2norm), lInfnorm
end
