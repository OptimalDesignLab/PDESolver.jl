# test boundary forces
function test_functional_integrate()
facts("--- Testing Boundary Functional Computation on CG Mesh ---") do
    clean_dict(arg_dict)
    ARGS[1] = "input_vals_functional_CG.jl"
    include("../../src/solver/advection/startup_advection.jl")  # initialization and construction
    println("use_DG = ", arg_dict["use_DG"])


    @fact mesh.isDG --> false
    @fact opts["functional_name1"] --> "qflux"
    @fact opts["analytical_functional_val"] --> roughly(2*(exp(1) - 1), atol=1e-12)
    @fact opts["geom_edges_functional1"] --> [2,3]
    
    fname = "./functional_error1.dat"
    error = readdlm(fname)

    @fact error[1] --> roughly(0.0060826244541961885, atol=1e-6)

    rm("./functional_error1.dat") # Delete the file
  end

  facts("--- Testing Boundary Functional Computation on DG Mesh ---") do
    ARGS[1] = "input_vals_functional_DG.jl"
    include("../../src/solver/advection/startup_advection.jl")  # initialization and construction


    @fact mesh.isDG --> true
    @fact opts["functional_name1"] --> "qflux"
    @fact opts["analytical_functional_val"] --> roughly(2*(exp(1) - 1), atol=1e-12)
    @fact opts["geom_edges_functional1"] --> [2,3]
    
    fname = "./functional_error1.dat"
    error = readdlm(fname)

    @fact error[1] --> roughly(1.657574600175115e-5, atol=1e-6)

    rm("./functional_error1.dat") # Delete the file
  end

  facts("--- Testing Adjoint Computation on DG Mesh ---") do
    
    @fact mesh.isDG --> true
    @fact opts["calc_adjoint"] --> true

    fname = "./adjoint_vector.dat"
    adjoint_vec = readdlm(fname)
    for i = 1:length(adjoint_vec)
      @fact adjoint_vec[i] --> roughly(-1.0 , atol=1e-10)
    end

    rm("./adjoint_vector.dat")
    
  end

  return nothing
end

#test_functional_integrate()
add_func1!(AdvectionTests, test_functional_integrate)
