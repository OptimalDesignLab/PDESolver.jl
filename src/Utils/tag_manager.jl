# managment of MPI tags

"""
  Manages MPI tags.  Specifically, keeps track of which tags are in use and
  provides for getting unused tags.  Currently it does not consider which tags
  are in use with which receiving processes.

  **Implementation Notes**

  The management scheme attempts to re-use tags as soon as they are freed.
  This requires an O(N) search in some cases to figure out the next tag.
  This could be optimized in various ways, but it doesn't seem worth it
  considering only a handful of tags are used at any one time.
"""
mutable struct MPITagManager
  used_tags::BitArray{1}  # true if tag in inuse, false otherwise
                          # if large numbers of tags are used, this could be
                          # a SparseVector{Bool} instead
  next_tag::Cint
  start_tag::Cint
end

"""
"""
function MPITagManager(starting_tag::Integer=1)

  used_tags = BitArray(0)
  next_tag = starting_tag

  return MPITagManager(used_tags, next_tag, starting_tag)
end

"""
  Returns the next free tag

  **Inputs**

   * mgr: an MPITagManager

  **Outputs**

   * a Cint containing the tag
"""
function getNextTag(mgr::MPITagManager)

  next_tag = mgr.next_tag
  next_tag_offset = next_tag - mgr.start_tag + 1

  resize_array = next_tag_offset > length(mgr.used_tags)
  if resize_array
    resize!(mgr.used_tags, next_tag_offset)
  end

  mgr.used_tags[next_tag_offset] = true
  
  # figure out what the next tag will be
  if resize_array
    new_next_tag = next_tag + 1
  else  
    # find next tag using O(N) search
    # hopefully this is the rare case because tags mostly get used and rarely
    # get freed
    new_next_tag = findfirst(mgr.used_tags, false) + mgr.start_tag - 1
    if new_next_tag == mgr.start_tag - 1
      new_next_tag = length(mgr.used_tags) + mgr.start_tag
    end
  end

  mgr.next_tag = new_next_tag

  return next_tag
end

"""
  Marks a tag as being in use.  This is useful if (for some reason) the user
  wants to use a specific tag number rather than calling [`getNextTag`](@ref)

  **Inputs**

   * mgr: [`MPITagManager`](@ref)
   * tag: an integer specifying the tag to be marked.  This tag must be
          greater than the starting tag passed to the `mgr` constructor
"""
function markTagUsed(mgr::MPITagManager, tag::Integer)
  
  if tag < mgr.start_tag
    error("cannot use tags less than start tag: $(mgr.start_tag)")
  end

  tag_offset = tag - mgr.start_tag + 1
  if tag_offset > length(mgr.used_tags)
    resize!(mgr.used_tags, tag_offset)
  end

  if mgr.used_tags[tag_offset]
    error("attempting to mark already-used tag $tag as used")
  end

  mgr.used_tags[tag_offset] = true

  return nothing
end

"""
  Marks a tag to be freed.  It is an error to free a tag that is not used.

  **Inputs**

   * tag: integer specifying the tag to be freed
"""
function freeTag(mgr::MPITagManager, tag::Integer)

  tag_offset = tag - mgr.start_tag + 1
  if tag < mgr.start_tag
    error("cannot use tags less than start tag: $(mgr.start_tag)")
  end

  if !mgr.used_tags[tag_offset]
    error("cannot free tag $tag, which is not in use")
  end

  mgr.used_tags[tag_offset] = false

  # figure out what tag to use next.  If possible, reuse this tag, but only
  # if the next tag is greater
  # Basically, try to use the first available tag
  if mgr.next_tag > tag
    mgr.next_tag = tag
  end

  return nothing
end



