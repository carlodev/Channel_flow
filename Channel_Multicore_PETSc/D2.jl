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

L = 2 # Domain length in each space dimension
D = 3 # Number of spatial dimensions
n = 8 # Partition (i.e., number of cells per space dimension)
function main_ex1(parts)
    function stretching(x::Point)
        m = zeros(length(x))
        m[1] = x[1]^2
        for i in 2:D
            m[i] = x[i]
        end
        Point(m)
    end

    pmin = Point(Fill(0, D))
    pmax = Point(Fill(L, D))
    partition = Tuple(Fill(n, D))
    model = CartesianDiscreteModel(parts, pmin, pmax, partition, map=stretching,  isperiodic=(true,false,true))
    #model = CartesianDiscreteModel(pmin, pmax, partition, map=stretching,  isperiodic=(true,false,true))

    writevtk(model,"model")

    return model

end

partition = (1, 4, 1)
prun(main_ex1, mpi, partition)