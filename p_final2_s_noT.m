function p_result = p_final2_s_noT(rhop,rhoa2,P_test)
% Calculate the osmotic pressure 
% No T means no need for the Bjerrum length input.
PEsolution = INPUT;
transform = PEsolution.transform;

f=free_energy2_noT(rhop,rhoa2);
kk=amu2_s(rhop,rhoa2,'all');
amu_p=kk(:,1);
amu_a2=kk(:,2);
p_result = -(f'-amu_p.*rhop'-amu_a2.*rhoa2').*transform-P_test;
end