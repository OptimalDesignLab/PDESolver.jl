# Parallel Constructs Documentations

```@meta
  CurrentModule = Utils
```

These function define the primative operations used by the physics modules to
exchange data in parallel.
When using these functions, the should not have to make any MPI calls directly,
they should all be encapsulated within the provided functions.

TODO: crossref to physics module documentation

The [Types and Basic API](@ref) section describes the [`SharedFaceData`](@ref)
datatype and the basic functions that operate on it.
The [Parallel Data Exchange](@ref) section describes the functions used
by the physics modules that start and finish parallel parallel communication.


## Types and Basic API

```@autodocs
  Modules = [Utils]
  Pages = ["Utils/parallel_types.jl"]
```

## Parallel Data Exchange

The functions in this section are used to start sending data in parallel and
finish receiving it.
All functions operate on a `Vector` of [`SharedFaceData`](@ref) that define
what data to send to which peer processes.  See [Parallel Overview](@ref) for
a high-level overview of how the code is parallelized.

Sending the data to the other processes is straight-forward.  Receiving it
(efficiently) is not.
In particular, [`finishExchangeData`] waits to receive data from one peer
process, calls a user supplied callback function to do calculations
involving the received data, and then waits for the next receive to finish.
This is significantly more efficient than waiting for all receives to finish
and then doing computations on all the data.

This section describes the API the physics modules use to do parallel 
communication.  The [Internals](@ref utils_parallel_internals) section
describes the helper functions used in the implementation.

```@docs
startSolutionExchange
exchangeData
finishExchangeData
TAG_DEFAULT
```


### [Internals](@id utils_parallel_internals)

These helper functions are used by the functions in
[Parallel Data Exchange](@ref).

```@docs
verifyReceiveCommunication
getSendDataFace
getSendDataElement
@mpi_master
@time_all
```
