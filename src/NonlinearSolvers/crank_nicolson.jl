# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson

@doc """
crank_nicolson

"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        q_vec::AbstractVector, res_vec::AbstractVector, pre_func,
                        post_func, ctx, opts)

end
