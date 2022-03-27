using Plots

N = 32;
gamma1=2.5;
y = zeros(N+1);
z = zeros(N+1)
y2 = zeros(N+1);

for i = 1:1:N+1
    y[i]= -tanh(gamma1*(1-2*i/(N)))/tanh(gamma1);

end

x = LinRange(-1, 1, N+1)

plot(z,y, marker=(:o))

for i = 1:1:N+1
    y2[i]= -tanh(gamma1*(x[i]))/tanh(gamma1);

end
plot!(z, y2, marker=(:ro))