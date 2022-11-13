% Plot the epsilon_r of water versus temperatures. It just a plof function
% do not effect the main function

f = @(x) 87.740-0.40008*x+9.398e-4*x.^2-1.410e-6*x.^3;
T = 0:1:100;
plot(T,f(T))