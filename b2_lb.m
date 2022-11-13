function length=b2_lb(T,sigma)
% Find Bjerrum length at given temperature
e = 1.602176634*1e-19;
Epsilon_0 = 8.854187817*1e-12;
k_B = 1.38064852*1e-23;
Epsilon_r = 80;

length = e.^2./(4.*pi.*Epsilon_0.*Epsilon_r.*k_B.*T)./(sigma.*1e-10);
end