using JLD2
using DataFrames, XLSX, CSV
using Gridap
using Plots
DNS = DataFrame(XLSX.readtable("Bazilev_Channel_DNS.xlsx", "Sheet1")...)
@load "example.jld2"

Ret = 395
ν = 0.0001472
fx = 0.00337204
ut = (ν * Ret * fx)^(1/3)

nn = 100
ev_nodes = LinRange(0.5, 0.99999, nn)
y =1
vt(i) = VectorValue(y,i,0)
ve(v) = v[1] 

vm = evaluate(uh, vt.(ev_nodes))
ev_nodes = (1 .-ev_nodes) .* ut./ ν
plot( ev_nodes, ve.(vm),  seriestype = :scatter)

plot!(DNS.y,DNS.U,  seriestype = :scatter)

