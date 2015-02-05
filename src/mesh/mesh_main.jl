# Main file contianing mesh related functions and types

# put module declaration here, if this is a module

# put include statements here





abstract AbstractMesh  # declare an abstract type called Mesh

# define concrete types here
type Mesh <: AbstractMesh
  # define some data  members here
  num_nodes::Int
  num_elements::Int
  num_edges::Int
  num_faces::Int

  # functions used to interact with the mesh data structure
  getnode_num::Function  # get information about node based on its number
  getnode_next::Function  # used to iterate over nodes
  getelement_num::Function # get information about element based on its number
  getelement_next::Function # used to iterate over elements
  getedge_num::Function
  getedge_next::Function
  getface_num::Function
  getface_next::Function

  # use default constructor, for now
end


#create some functions to interface with Simmetrix, or create a local mesh
# when constructing a Mesh object, the appropriate functions are passed to it
# this enforces a common interface for all mesh types
