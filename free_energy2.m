function f = free_energy2(rhop,rhoa2,lb)
%---------------------INPUT----------------
PEsolution = INPUT;

zp_n = PEsolution.zp_n;
zp_e = PEsolution.zp_e;
za2 = PEsolution.za2;

Np = PEsolution.Np;
Na2 = PEsolution.Na2;

N1 = PEsolution.N1;
N2 = PEsolution.N2;

lamda_e = PEsolution.lamda_e; 
lamda_n = PEsolution.lamda_n;

mp = PEsolution.MPCE;
ma2 = PEsolution.MCa;
rhop_e = rhop.*lamda_e;
rhop_n = rhop.*lamda_n;
Bk = PEsolution.Bk;
Ak = PEsolution.Ak;

% In order to use PC-SAFT to describe the soft sphere model, two diameters
% should be used in different parts of free energy density
% Hard-sphere diameter should be used in PC-SAFT and soft-sphete diameter
% should be used in other models.
% sigma_i should be soft-sphere diameter and sigmai_v should be hard-sphere
% diameter.
% Repeat the same procedure in amu2_s, b2_amu and free_energy2_noT.

sigmap = PEsolution.sigmap_av;
sigmaa2 = PEsolution.sigmaa2;
sigmap_v = PEsolution.sigmap;
sigmaa2_v = PEsolution.sigmaa2;
sigmap_a2_v = (sigmap_v+sigmaa2_v)/2; 

T = PEsolution.T;
ep_a2 = PEsolution.epce_ca;
epM = PEsolution.eppM*T;
epA = PEsolution.eppA*T;
%% ------------- the ideal part of free energy density-------------
fid = rhop/Np.*(log(rhop/Np)-1)+rhoa2/Na2.*(log(rhoa2/Na2)-1);

%% ------------- the hard-sphere part of free energy density-------
n0p = rhop;
n0a2 = rhoa2;
n0p_e = rhop_e;
n0p_n = rhop_n;

n1p_e = rhop_e*sigmap/2;
n1p_n = rhop_n*sigmap/2;
n1a2 =rhoa2*sigmaa2/2;

n2p_e = pi*rhop_e*sigmap^2;
n2p_n = pi*rhop_n*sigmap^2;
n2a2 =pi*rhoa2*sigmaa2^2;


n3p_e = pi*rhop_e*sigmap^3/6;
n3p_n= pi*rhop_n*sigmap^3/6;
n3a2 =pi*rhoa2*sigmaa2^3/6;


n0 = n0p_e+n0p_n+n0a2;
n1 = n1p_e+n1p_n+n1a2;
n2 = n2p_e+n2p_n+n2a2;
n3 = n3p_e+n3p_n+n3a2;


fhs = -1.*n0.*log(1-n3)+(n1.*n2)./(1-n3)+(n2.^3./(36.*pi.*n3.^3)).*...
(n3.*log(1-n3)+n3.^2./((1-n3).^2));

% gamma (screening parameter) & Pn (size asymmetry factor)
Ga = 0;
Pn = 0;
for times = 1:20
gammat_a2 = (rhoa2./((1+Ga.*sigmaa2).^2)).*(za2-(pi.*Pn.*sigmaa2.^2)./(2.*(1-n3))).^2;
gammat_p_e = (rhop_e./((1+Ga.*sigmap).^2)).*(zp_e-(pi.*Pn.*sigmap.^2)./(2.*(1-n3))).^2;
gammat_p_n = (rhop_n./((1+Ga.*sigmap).^2)).*(zp_n-(pi.*Pn.*sigmap.^2)./(2.*(1-n3))).^2;

gamma = pi.*lb.*(gammat_p_e+gammat_p_n+gammat_a2);
Ga= sqrt(gamma);

pnt1p1_e = rhop_e.*sigmap.*zp_e./(1+Ga.*sigmap);
pnt1p1_n = rhop_n.*sigmap.*zp_n./(1+Ga.*sigmap);
pnt11_2 = rhoa2.*sigmaa2.*za2./(1+Ga.*sigmaa2);

pnt2 = rhop_e.*sigmap.^3./(1+Ga.*sigmap)+...
    rhoa2.*sigmaa2.^3./(1+Ga.*sigmaa2)+ ...
    rhop_n.*sigmap.^3./(1+Ga.*sigmap);
pnt3 = 1+pi./(2.*(1-n3)).*pnt2;
Pn = pnt1p1_e./pnt3+pnt1p1_n./pnt3+pnt11_2./pnt3;
end
a1 = 2.*pi.*lb.*(zp_e-(pi.*Pn.*sigmap.^2)./(2.*(1-n3)))./(Ga.*(1+Ga*sigmap));
%% --------- the chain connectivity part of free energy density----
ghs_pp = 1./(1-n3) +n2.*sigmap./(4.*(1-n3).^2) ;
y_pp_1 = (-Ga.^2.*a1.^2./(4.*pi.^2.*sigmap.*lb)+(lb.*zp_e.^2./sigmap));
fch = -1/Np.*n0p.*(N1.*y_pp_1+(N2+N1).*log(ghs_pp));
%% ---the electrostatic correlation part of free energy density ---
felt11_2 = (n0a2.*za2)./(1+Ga.*sigmaa2);
felt1p1_e = (n0p_e.*zp_e)./(1+Ga.*sigmap);
felt1p1_n = (n0p_n.*zp_n)./(1+Ga.*sigmap);


felt21_2 = Ga.*za2+(pi.*Pn.*sigmaa2)./(2.*(1-n3));
felt2p1_e = Ga.*zp_e+(pi.*Pn.*sigmap)./(2.*(1-n3));
felt2p1_n = Ga.*zp_n+(pi.*Pn.*sigmap)./(2.*(1-n3));

fel = -1.*lb.*(felt1p1_e.*felt2p1_e+...
         felt11_2.*felt21_2+felt1p1_n.*felt2p1_n)+Ga.^3./(3.*pi);
%% ---------the calcium-binding part of free energy density-----
n3 = n3p_e+n3a2+n3p_n;
n0p = n0p_e+n0p_n;
x=(n0a2/ma2)./(n0p/mp+n0a2/ma2);
m=ma2*x+mp*(1-x);
a0=Ak(1,1)+(m-1)./m*Ak(1,2)+(m-1).*(m-2)./m.^2*Ak(1,3);
a1=Ak(2,1)+(m-1)./m*Ak(2,2)+(m-1).*(m-2)./m.^2*Ak(2,3);
a2=Ak(3,1)+(m-1)./m*Ak(3,2)+(m-1).*(m-2)./m.^2*Ak(3,3);
a3=Ak(4,1)+(m-1)./m*Ak(4,2)+(m-1).*(m-2)./m.^2*Ak(4,3);
a4=Ak(5,1)+(m-1)./m*Ak(5,2)+(m-1).*(m-2)./m.^2*Ak(5,3);
a5=Ak(6,1)+(m-1)./m*Ak(6,2)+(m-1).*(m-2)./m.^2*Ak(6,3);
a6=Ak(7,1)+(m-1)./m*Ak(7,2)+(m-1).*(m-2)./m.^2*Ak(7,3);
I1=a0+a1.*n3+a2.*n3.^2+a3.*n3.^3+a4.*n3.^4+a5.*n3.^5+a6.*n3.^6;
fdsI=-2*pi*I1.*(2*n0p_e.*n0a2*ep_a2/T*sigmap_a2_v^3+n0p_e.^2*epA/T*sigmap_v^3+n0p_n.^2*epM/T*sigmap_v^3);

b0=Bk(1,1)+(m-1)./m*Bk(1,2)+(m-1).*(m-2)./m.^2*Bk(1,3);
b1=Bk(2,1)+(m-1)./m*Bk(2,2)+(m-1).*(m-2)./m.^2*Bk(2,3);
b2=Bk(3,1)+(m-1)./m*Bk(3,2)+(m-1).*(m-2)./m.^2*Bk(3,3);
b3=Bk(4,1)+(m-1)./m*Bk(4,2)+(m-1).*(m-2)./m.^2*Bk(4,3);
b4=Bk(5,1)+(m-1)./m*Bk(5,2)+(m-1).*(m-2)./m.^2*Bk(5,3);
b5=Bk(6,1)+(m-1)./m*Bk(6,2)+(m-1).*(m-2)./m.^2*Bk(6,3);
b6=Bk(7,1)+(m-1)./m*Bk(7,2)+(m-1).*(m-2)./m.^2*Bk(7,3);
M=1+m.*(8*n3-2*n3.^2)./(1-n3).^4+(1-m).*(20*n3-27*n3.^2+12*n3.^3-2*n3.^4)./(1-n3).^2./(2-n3).^2;
I2=b0+b1.*n3+b2.*n3.^2+b3.*n3.^3+b4.*n3.^4+b5.*n3.^5+b6.*n3.^6;
fdsII=-pi./M.*I2.*(2*n0p_e.*n0a2*(ep_a2/T)^2*sigmap_a2_v^3+n0p_e.^2*(epA/T)^2*sigmap_v^3+n0p_n.^2*(epM/T)^2*sigmap_v^3).*m;
     %% the sum of all parts
f = fid + fch + fel + fhs+ fdsI + fdsII;

end







