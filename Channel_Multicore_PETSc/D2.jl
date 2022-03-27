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
D = 2 # Number of spatial dimensions
n = 4 # Partition (i.e., number of cells per space dimension)
function main_ex1(part)
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
    model = CartesianDiscreteModel(part, pmin, pmax, partition, map=stretching,  isperiodic=(true,false))
    writevtk(model,"model")

    return model

end

partition = (2, 2)
prun(main_ex1, mpi, partition)