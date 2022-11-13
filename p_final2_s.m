function p_result = p_final2_s(rhop,rhoa2,P_test)
% Calculate the osmotic pressure

PEsolution = INPUT;
lb =PEsolution.lb;
k_B = 1.38064852*1e-23;
sigma = PEsolution.sigma;
T = PEsolution.T;
transform = k_B*T/((sigma*1e-10)^3)*1e-6;
f=free_energy2(rhop,rhoa2,lb);

kk=amu2_s(rhop,rhoa2,'all');
amu_p=kk(:,1);
amu_a2=kk(:,2);
p_result = -(f'-amu_p.*rhop'-amu_a2.*rhoa2').*transform-P_test;
end