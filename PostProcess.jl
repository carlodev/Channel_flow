using JLD2
using DataFrames, XLSX, CSV
using Gridap
using Plots
using LaTeXStrings

DNS = DataFrame(XLSX.readtable("Bazilev_Channel_DNS.xlsx", "Sheet1")...)
@load "Channel_3d_station.jld2"




Ret = 395
ν = 0.0001472
fx = 0.00337204

ut = ν*Ret
nn = 1000

ev_nodes = ch_nodes[500:end]
log_nodes = LinRange(10, 500, 100)
visc_nodes = LinRange(0.1, 5, 100)

y =1
vt(i) = VectorValue(y,i,0)
ve(v) = v[1] 
χ=0.41
B=5.2

uplus(y) = 1 ./χ  .* log.(y) .+ B 
vm = ch_vm[500:end]
vm = vm ./ut /Um
ev_nodes = (1 .-ev_nodes) .* ut./ ν 
plot(ev_nodes, ve.(vm),  seriestype = :scatter, label="SUPG", xaxis=:log)
plot!(DNS.y, DNS.U,  seriestype = :scatter, label="DNS", xaxis=:log)
plot!(visc_nodes, visc_nodes,  linestyle=:dash, linewidth=3, label=L"U^+=y^+", xaxis=:log)
plot!(log_nodes, uplus(log_nodes),  linestyle=:dash, linewidth=3, label="LogLaw", xaxis=:log)
plot!(legend=:topleft)
savefig("Channel_pp.pdf")

