function D=eq_b2(A)
PEsolution = INPUT;
zp_c = abs(PEsolution.zp_c);
za2 = PEsolution.za2;
rhop1 = A(1);
rhop2 = A(2);
lb=PEsolution.lb;
mu_p1_b=b2_amu(rhop1,rhop1*zp_c/za2,'all',lb);
mu_p2_b=b2_amu(rhop2,rhop2*zp_c/za2,'all',lb);
D(1) = (mu_p1_b(1)-mu_p2_b(1))*za2+(mu_p1_b(2)-mu_p2_b(2))*zp_c;
D(2) = p_final2_s(rhop1,rhop1*zp_c/za2,0)-p_final2_s(rhop2,rhop2*zp_c/za2,0);



