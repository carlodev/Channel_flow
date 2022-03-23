using Plots

N = 32;
gamma1=2.5;
y = zeros(N);
z = zeros(N)
for i = 1:1:N
    y[i]= -tanh(gamma1*(1-2*i/(N-1)))/tanh(gamma1);

end

plot(z,y, marker=(:o))