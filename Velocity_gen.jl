
using Plots
nu = 14.88 *10^(-6)
rho = 1.225
mu = nu*rho

u0 = 10

dpdx = -2000/20
h = 1

y = LinRange(0,1,100)
ux = 0.5.*mu.*dpdx.*(y.^2-y.*h)
y1 = y .- 0.5
plot(y1,ux)

