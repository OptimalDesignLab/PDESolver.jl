push!(LOAD_PATH,"../src/simple_mesh/")

using PDESolver
using SimpleMesh

# Using fact check for the first time
facts("Testing SimpleMesh type") do
  # Enter your tests here
  context("Check if function createMesh() works or not") do
    # Build an element grid comprising of 2 triangular elements
    for p = 1:4
      if p == 1
        act_ien = [1 2;2 4;3 3]
        act_nnp = 4
        act_nel = 2
        act_vtx_loc = [0.0 1.0 0.0 1.0;
                       0.0 0.0 1.0 1.0]
        act_elem_edge_lengthx = 1.0
        act_elem_edge_lengthy = 1.0
        act_NodeEdgex = [1 2;3 4]
        act_NodeEdgey = [2 4;1 3]
        act_HBedges = [1;2]
        act_VBedges = [4;5]
        act_nedges = 5
        m = SimpleMesh{Float64}(1,1,1,1,p+1)
        res1 = act_ien - m.ien
        res2 = act_nnp - m.nnp
        res3 = act_nel - m.nel
        res4 = act_vtx_loc - m.vtx_loc
        res5 = act_elem_edge_lengthx - m.elem_edge_lengthx
        res6 = act_elem_edge_lengthy - m.elem_edge_lengthy
        res7 = act_NodeEdgex - m.NodeEdgex
        res8 = act_NodeEdgey - m.NodeEdgey
        res9 = act_HBedges - m.HBedges
        res10 = act_VBedges - m.VBedges
        res11 = act_nedges - m.nedges
        @fact res1 --> roughly(zeros(size(ien)), atol=1e-13)
        @fact res2 --> roughly(0, atol=1e-13)
        @fact res3 --> roughly(0, atol=1e-13)
        @fact res4 --> roughly(zeros(size(act_vtx_loc)), atol=1e-13)
        @fact res5 --> roughly(0, atol=1e-13)
        @fact res6 --> roughly(0, atol=1e-13)
        @fact res7 --> roughly(zeros(size(act_NodeEdgex)), atol=1e-13)
        @fact res8 --> roughly(zeros(size(act_NodeEdgey)), atol=1e-13)
        @fact res9 --> roughly(zeros(size(act_HBedges)), atol=1e-13)
        @fact res10 --> roughly(zeros(size(act_VBedges)), atol=1e-13)
        @fact res11 --> roughly(0, atol=1e-13)
      elseif p == 2
        act_ien = [1 3;3 9;7 7;2 6;5 8;4 5;10 11]
        act_nnp = 11
        act_nel  = 2
        act_vtx_loc = [0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                       0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0]
        act_elem_edge_lengthx = 1.0
        act_elem_edge_lengthy = 1.0
        act_NodeEdgex = [1 2 3;7 8 9]
        act_NodeEdgey = [3 6 9;1 4 7]
        act_HBedges = [1;2]
        act_VBedges = [4;5]
        act_nedges = 5
        m = SimpleMesh{Float64}(1,1,1,1,p+1)
        res1 = act_ien - m.ien
        res2 = act_nnp - m.nnp
        res3 = act_nel - m.nel
        res4 = act_vtx_loc - m.vtx_loc
        res5 = act_elem_edge_lengthx - m.elem_edge_lengthx
        res6 = act_elem_edge_lengthy - m.elem_edge_lengthy
        res7 = act_NodeEdgex - m.NodeEdgex
        res8 = act_NodeEdgey - m.NodeEdgey
        res9 = act_HBedges - m.HBedges
        res10 = act_VBedges - m.VBedges
        res11 = act_nedges - m.nedges
        @fact res1 --> roughly(zeros(size(ien)), atol=1e-13)
        @fact res2 --> roughly(0, atol=1e-13)
        @fact res3 --> roughly(0, atol=1e-13)
        @fact res4 --> roughly(zeros(size(act_vtx_loc)), atol=1e-13)
        @fact res5 --> roughly(0, atol=1e-13)
        @fact res6 --> roughly(0, atol=1e-13)
        @fact res7 --> roughly(zeros(size(act_NodeEdgex)), atol=1e-13)
        @fact res8 --> roughly(zeros(size(act_NodeEdgey)), atol=1e-13)
        @fact res9 --> roughly(zeros(size(act_HBedges)), atol=1e-13)
        @fact res10 --> roughly(zeros(size(act_VBedges)), atol=1e-13)
        @fact res11 --> roughly(0, atol=1e-13)
      elseif p == 3 
        act_ien = [1 4;4 14;11 11;2 7;3 10;6 13;9 12;5 6;8 9;15 18; 16 19;17 20]
        act_nnp = 20
        act_nel  = 2
        act_vtx_loc = [0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0]
        act_elem_edge_lengthx = 1.0
        act_elem_edge_lengthy = 1.0
        act_NodeEdgex = [1 2 3 4;11 12 13 14]
        act_NodeEdgey = [4 7 10 14;1 5 8 11]
        act_HBedges = [1;2]
        act_VBedges = [4;5]
        act_nedges = 5
        m = SimpleMesh{Float64}(1,1,1,1,p+1)
        res1 = act_ien - m.ien
        res2 = act_nnp - m.nnp
        res3 = act_nel - m.nel
        res4 = act_vtx_loc - m.vtx_loc
        res5 = act_elem_edge_lengthx - m.elem_edge_lengthx
        res6 = act_elem_edge_lengthy - m.elem_edge_lengthy
        res7 = act_NodeEdgex - m.NodeEdgex
        res8 = act_NodeEdgey - m.NodeEdgey
        res9 = act_HBedges - m.HBedges
        res10 = act_VBedges - m.VBedges
        res11 = act_nedges - m.nedges
        @fact res1 --> roughly(zeros(size(ien)), atol=1e-13)
        @fact res2 --> roughly(0, atol=1e-13)
        @fact res3 --> roughly(0, atol=1e-13)
        @fact res4 --> roughly(zeros(size(act_vtx_loc)), atol=1e-13)
        @fact res5 --> roughly(0, atol=1e-13)
        @fact res6 --> roughly(0, atol=1e-13)
        @fact res7 --> roughly(zeros(size(act_NodeEdgex)), atol=1e-13)
        @fact res8 --> roughly(zeros(size(act_NodeEdgey)), atol=1e-13)
        @fact res9 --> roughly(zeros(size(act_HBedges)), atol=1e-13)
        @fact res10 --> roughly(zeros(size(act_VBedges)), atol=1e-13)
        @fact res11 --> roughly(0, atol=1e-13)
      elseif p == 4
        act_ien = [1 5;5 19;15 15;2 8;3 11;4 14;7 18;10 17;13 16;6 7;9 10;12 13;20 26;21 27;
                   22 28;23 29;24 30;25 31]
        act_nnp = 31
        act_nel  = 2
        act_vtx_loc = [0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
        act_elem_edge_lengthx = 1.0
        act_elem_edge_lengthy = 1.0
        act_NodeEdgex = [1 2 3 4 5;15 16 17 18 19]
        act_NodeEdgey = [5 8 11 14 19;1 6 9 12 15]
        act_HBedges = [1;2]
        act_VBedges = [4;5]
        act_nedges = 5
        m = SimpleMesh{Float64}(1,1,1,1,p+1)
        res1 = act_ien - m.ien
        res2 = act_nnp - m.nnp
        res3 = act_nel - m.nel
        res4 = act_vtx_loc - m.vtx_loc
        res5 = act_elem_edge_lengthx - m.elem_edge_lengthx
        res6 = act_elem_edge_lengthy - m.elem_edge_lengthy
        res7 = act_NodeEdgex - m.NodeEdgex
        res8 = act_NodeEdgey - m.NodeEdgey
        res9 = act_HBedges - m.HBedges
        res10 = act_VBedges - m.VBedges
        res11 = act_nedges - m.nedges
        @fact res1 --> roughly(zeros(size(ien)), atol=1e-13)
        @fact res2 --> roughly(0, atol=1e-13)
        @fact res3 --> roughly(0, atol=1e-13)
        @fact res4 --> roughly(zeros(size(act_vtx_loc)), atol=1e-13)
        @fact res5 --> roughly(0, atol=1e-13)
        @fact res6 --> roughly(0, atol=1e-13)
        @fact res7 --> roughly(zeros(size(act_NodeEdgex)), atol=1e-13)
        @fact res8 --> roughly(zeros(size(act_NodeEdgey)), atol=1e-13)
        @fact res9 --> roughly(zeros(size(act_HBedges)), atol=1e-13)
        @fact res10 --> roughly(zeros(size(act_VBedges)), atol=1e-13)
        @fact res11 --> roughly(0, atol=1e-13)
      end # Ends if-else statement
    end   # Ends for loop
  end     # Ends context

end # Ends facts
