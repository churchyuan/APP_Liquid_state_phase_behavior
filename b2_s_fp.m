function Ans=b2_s_fp(rhop1,rhop2)
PEsolution = INPUT;
lb = PEsolution.lb;
zp_c = abs(PEsolution.zp_c);
za2 = PEsolution.za2;
rhoa1 = rhop1*zp_c/za2;
rhoa2 = rhop2*zp_c/za2;
mu_p1_b=amu2_s(rhop1,rhoa1,'all');
mu_p2_b=amu2_s(rhop2,rhoa2,'all');

mu_p1 = mu_p1_b(:,1)*za2+mu_p1_b(:,2)*zp_c;
mu_p2 = mu_p2_b(:,1)*za2+mu_p2_b(:,2)*zp_c;
p_p1 = p_final2_s_noT(rhop1,rhoa1,0);
p_p2 = p_final2_s_noT(rhop2,rhoa2,0); 
n=1;
p1=polyfit(mu_p1,p_p1,n);
p2=polyfit(mu_p2,p_p2,n);
syms x3;
f1=poly2sym(p1,x3);
f2=poly2sym(p2,x3);
x=vpasolve(f1==f2);
p1=abs(abs(mu_p1)-abs(x));
p2=abs(abs(mu_p2)-abs(x));
[~,b1]=min(p1);
[~,b2]=min(p2);
Ans=[rhop1(b1),rhop2(b2)];
Ans = [Ans,lb];

end