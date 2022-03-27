using Gridap
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using FillArrays
using Test
using InteractiveUtils
using GridapDistributed
using PartitionedArrays

function mesh_channel(; D::Integer, N=32::Integer, parts=1, printmodel=false::Bool)

   """
   mesh_channel() generate a mesh for a channel; 
   Periodic boundaries in dimensions 1 and 3
   In dimensions 1 and 3 equally spaced
   In dimension 2 function distributed
   #Arguments
   - D::Integer number of dimensions (1 or 2)
   - N::Integer numer of cells in each dimension, deault value N=32
   - parts :: if distributed
   - printmodel::Boolean if true create vtk file pf the model
   """


   #D = 2 # Number of spatial dimensions

   Lx = 2 * pi
   Ly = 2
   Lz = 2 / 3 * pi
   nx = N
   ny = N
   nz = N


   #N = 32 # Partition (i.e., number of cells per space dimension)
   function stretching(x::Point)
      m = zeros(length(x))
      m[1] = x[1]

      gamma1 = 2.5
      m[2] = -tanh(gamma1 * (x[2])) / tanh(gamma1)
      if length(x) > 2
         m[3] = x[3]
      end
      Point(m)
   end

   if D > 2
      pmin = Point(0, -Ly / 2, -Lz / 2)
      pmax = Point(Lx, Ly / 2, Lz / 2)
      partition = (nx, ny, nz)
      periodic_tuple = (true, false, true)
      model_name = "model3d"

   else
      pmin = Point(0, -Ly / 2)
      pmax = Point(Lx, Ly / 2)
      partition = (nx, ny)
      periodic_tuple = (true, false)
      model_name = "model2d"

   end

   #partition = Tuple(Fill(N, D))
   if parts != 1
   model = CartesianDiscreteModel(parts, pmin, pmax, partition, map=stretching, isperiodic=periodic_tuple)
   else
      model = CartesianDiscreteModel(pmin, pmax, partition, map=stretching, isperiodic=periodic_tuple)
   end

   if printmodel
      writevtk(model, model_name)
   end
   return   model
end

