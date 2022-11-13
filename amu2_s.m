function amu = amu2_s(rhop,rhoa2,a)
% Function that calculates chemical potential of the salt-free solution
% rhop and rhoa2 denote the concentration of polymer and counterion,
% respectively.
PEsolution = INPUT;
lb = PEsolution.lb;
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
ep_a2 = PEsolution.epce_ca;

% In order to use PC-SAFT to describe the soft sphere model, two diameters
% should be used in different parts of free energy density.
% Hard-sphere diameter should be used in PC-SAFT and soft-sphete diameter
% should be used in other models.
% sigma_i should be soft-sphere diameter and sigmai_v should be hard-sphere
% diameter.
% Repeat the same procedure in b2_amu, free_energy and free_energy_noT.

sigmap = PEsolution.sigmap_av;
sigmaa2 = PEsolution.sigmaa2;
sigmap_v = PEsolution.sigmap;
sigmaa2_v = PEsolution.sigmaa2;
sigmap_a2_v = (sigmap_v+sigmaa2_v)/2; 
T = PEsolution.T;
epM = PEsolution.eppM*T;
epA = PEsolution.eppA*T;
%% ----- the ideal part of chemical potential-----------------------
d_fid_p = log(rhop./Np)./Np;
d_fid_a2 = log(rhoa2);

%% ----- the hard-sphere part of chemical potential-----------------
% Parameter of Rosenfeld
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

d_fhs_f = @(sigma) -log(1-n3)+sigma./2.*n2./((1-n3)) + pi.*sigma.^2.*(n1./(1-n3)+ ...
          n2.^2.*log(1-n3)./(12.*pi.*n3.^2)+n2.^2./(12.*pi.*n3.*(1-n3).^2))+ ...
          pi./6.*sigma.^3.*(n0./(1-n3)+n1.*n2./(1-n3).^2-n2.^3.*log(1-n3)./(18.*pi.*n3.^3) - ...
          n2.^3./(36.*pi.*n3.^2.*(1-n3))-n2.^3.*(1-3.*n3)./(36.*pi.*n3.^2.*(1-n3).^3));

d_fhs_p = d_fhs_f(sigmap);
d_fhs_a2 = d_fhs_f(sigmaa2);
%% --------- the chain connectivity part of chemical potential----

Ga = 0;
Pn = 0;
for times = 1:30
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

a1 = 2.*pi.*lb.*(zp_e-(pi.*Pn.*sigmap.^2)./(2.*(1-n3)))./(Ga.*(1+Ga.*sigmap));
ghs_pp = 1./(1-n3) +n2.*sigmap./(4.*(1-n3).^2) ;
y_pp = ghs_pp.*exp(-Ga.^2.*a1.^2./(4.*pi.^2.*sigmap.*lb)+(lb.*zp_e.^2./sigmap));

% Functions that simplify addition 
df1_dPn_f_sum =@(n0i,sigma,zi) -pi.*lb.*(n0i).*(zi-(pi.*Pn.*sigma^2)./(2.*(1-n3)))./(1+Ga.*sigma).* ...
    (pi.*sigma^2)./(1-n3)./(1+Ga.*sigma);
df1_dPn = df1_dPn_f_sum(n0p_e,sigmap,zp_e)+ ...
          df1_dPn_f_sum(n0a2,sigmaa2,za2)+ ...
          df1_dPn_f_sum(n0p_n,sigmap,zp_n);
      
df1_dGa_f_sum = @(n0i,sigma,zi) -pi.*lb.*2.*sigma.*(n0i).*(zi-(pi.*Pn.*sigma^2)./(2.*(1-n3))).^2./(1+Ga.*sigma).^3;
df1_dGa = df1_dGa_f_sum(n0p_e,sigmap,zp_e)+ ...
          df1_dGa_f_sum(n0a2,sigmaa2,za2)+ ...
          df1_dGa_f_sum(n0p_n,sigmap,zp_n);
      
df1_dn0i_f = @(sigma,zi) pi.*lb.*( (zi-(pi.*Pn.*sigma^2)./(2.*(1-n3)))./(1+Ga.*sigma) ).^2;


df1_dn3i_f_sum = @(n0i,sigma,zi) -pi.*lb.*(n0i).*(zi-(pi.*Pn.*sigma^2)./(2.*(1-n3)))./(1+Ga.*sigma).* ...
            (pi.*Pn.*sigma^2)./(1-n3).^2./(1+Ga.*sigma);
df1_dn3i = df1_dn3i_f_sum(n0p_e,sigmap,zp_e)+ ...
           df1_dn3i_f_sum(n0p_n,sigmap,zp_n)+ ...
           df1_dn3i_f_sum(n0a2,sigmaa2,za2);

df2_dGa_f_p1 = @(n1i,sigma,zi) 2.*sigma.*n1i.*zi./(1+Ga.*sigma).^2 ;
df2_dGa_f_p2 = @(n3i,sigma) n3i./(1+Ga.*sigma);
df2_dGa_f_p3 = @(n1i,sigma,zi) 2.*n1i.*zi./(1+Ga.*sigma) ;
df2_dGa_f_p4 = @(n3i,sigma,zi) n3i.*sigma./(1+Ga.*sigma).^2 ;

df2_dGa_p1 = (df2_dGa_f_p1(n1p_e,sigmap,zp_e)+...
               df2_dGa_f_p1(n1p_n,sigmap,zp_n)+ ...
               df2_dGa_f_p1(n1a2,sigmaa2,za2));      
df2_dGa_p2 = (df2_dGa_f_p2(n3p_e,sigmap)+ ...
                df2_dGa_f_p2(n3p_n,sigmap)+...
              df2_dGa_f_p2(n3a2,sigmaa2));
df2_dGa_p3 = (df2_dGa_f_p3(n1p_e,sigmap,zp_e)+...
              df2_dGa_f_p3(n1p_n,sigmap,zp_n)+...
               df2_dGa_f_p3(n1a2,sigmaa2,za2));     
df2_dGa_p4= (df2_dGa_f_p4(n3p_e,sigmap,zp_e)+...
             df2_dGa_f_p4(n3p_n,sigmap,zp_n)+...
               df2_dGa_f_p4(n3a2,sigmaa2,za2));   
           
df2_dGa = -df2_dGa_p1./(1+3./(1-n3).*df2_dGa_p2)+3./(1-n3).*df2_dGa_p3.* ...
    df2_dGa_p4./(1+3./(1-n3).*df2_dGa_p2).^2; 
         
df2_dn1i_f = @(zi,sigma) 2.*zi./(1+Ga.*sigma)./(1+3./(1-n3).*df2_dGa_p2 );

df2_dn3i_f_p1 = @(n1i,sigma,zi) 2.*n1i.*zi./(1+Ga.*sigma) ;
df2_dn3i_f_p2 = @(n3i,sigma) n3i./(1+Ga.*sigma);
df2_dn3i_p1 = (df2_dn3i_f_p1(n1p_e,sigmap,zp_e)+...
                df2_dn3i_f_p1(n1p_n,sigmap,zp_n)+...
               df2_dn3i_f_p1(n1a2,sigmaa2,za2));
df2_dn3i_p2 = (df2_dn3i_f_p2(n3p_e,sigmap)+...
                df2_dn3i_f_p2(n3p_n,sigmap)+...
               df2_dn3i_f_p2(n3a2,sigmaa2));
df2_dn3i_f = @(sigma) -df2_dn3i_p1./(1+3./(1-n3).*df2_dn3i_p2).^2.*(3./(1-n3).^2.*df2_dn3i_p2+ ...
             3./((1-n3).*(1+Ga.*sigma)));

df3_dPn = -pi^2.*lb.*sigmap^2./(Ga.*(1+Ga.*sigmap).*(1-n3)) ;

df3_dGa = -2.*pi.*lb.*(zp_e-pi.*Pn.*sigmap.^2./(2.*(1-n3))).*(1+2.*Ga.*sigmap)./ ...
          (Ga.^2.*(1+Ga.*sigmap).^2);
df3_dn3i =  -pi^2.*lb.*Pn.*sigmap.^2./(Ga.*(1+Ga.*sigmap).*(1-n3).^2); 


dGa_dn0i = @(sigma,zi) df1_dn0i_f(sigma,zi)./(2.*Ga - df1_dPn.*df2_dGa-df1_dGa);
dPn_dn0i = @(sigma,zi) df2_dGa.*dGa_dn0i(sigma,zi) ;
da1_dn0i = @(sigma,zi) df3_dPn.*dPn_dn0i(sigma,zi) + df3_dGa.*dGa_dn0i(sigma,zi) ;


dGa_dn1i = @(zi,sigma) df1_dPn.*df2_dn1i_f(zi,sigma)./(2.*Ga - df1_dPn.*df2_dGa-df1_dGa);
dPn_dn1i = @(zi,sigma) df2_dn1i_f(zi,sigma)+df2_dGa.*dGa_dn1i(zi,sigma) ;
da1_dn1i = @(zi,sigma) df3_dPn.*dPn_dn1i(zi,sigma) + df3_dGa.*dGa_dn1i(zi,sigma) ;


dGa_dn3i = @(sigma) (df1_dn3i+df1_dPn.*df2_dn3i_f(sigma))./(2.*Ga - df1_dPn.*df2_dGa-df1_dGa);
dPn_dn3i = @(sigma) df2_dn3i_f(sigma)+df2_dGa.*dGa_dn3i(sigma) ;
da1_dn3i = @(sigma) df3_dn3i+df3_dPn.*dPn_dn3i(sigma) + df3_dGa.*dGa_dn3i(sigma) ;


dy_samep1 = exp(lb.*zp_e.^2./sigmap-Ga.^2.*a1.^2./(4.*pi.^2.*sigmap.*lb));
dy_pp_dn0i = @(sigma,zi) -ghs_pp.*dy_samep1.*(Ga.^2.*a1./(2.*pi.^2.*sigmap.*lb).*da1_dn0i(sigma,zi)+ ...
             Ga.*a1.^2./(2.*pi.^2.*sigmap.*lb).*dGa_dn0i(sigma,zi));
dy_pp_dn1i = @(zi,sigma) -ghs_pp.*dy_samep1.*(Ga.^2.*a1./(2.*pi.^2.*sigmap.*lb).*da1_dn1i(zi,sigma)+ ...
             Ga.*a1.^2./(2.*pi.^2.*sigmap.*lb).*dGa_dn1i(zi,sigma));
dghs_pp_dn2i = sigmap./(4.*(1-n3).^2);
dghs_pp_dn3i = 1./(1-n3).^2+n2.*sigmap./(2.*(1-n3).^3);
dy_pp_dn2i = dghs_pp_dn2i.*dy_samep1;                  
dy_pp_dn3i = @(sigma) dghs_pp_dn3i.*dy_samep1-ghs_pp.*dy_samep1.*(Ga.^2.*a1./(2.*pi.^2.*sigmap.*lb).*da1_dn3i(sigma)+ ...
             Ga.*a1.^2./(2.*pi.^2.*sigmap.*lb).*dGa_dn3i(sigma));

% Chain connectivity for polymer 
dPhi_dn0i = -1./Np.*(N1.*log(y_pp)+N2.*log(ghs_pp))-N1./Np.*n0p./y_pp.*(dy_pp_dn0i(sigmap,zp_e).*lamda_e+dy_pp_dn0i(sigmap,zp_n).*lamda_n);
dPhi_dn1i = -N1./Np.*n0p./y_pp.*(dy_pp_dn1i(zp_e,sigmap).*lamda_e+dy_pp_dn1i(zp_n,sigmap).*lamda_n);
dPhi_dn2i = -N1./Np.*n0p./y_pp.*dy_pp_dn2i-N2./Np.*n0p./ghs_pp.*dghs_pp_dn2i;
dPhi_dn3i = -N1./Np.*n0p./y_pp.*dy_pp_dn3i(sigmap)-N2./Np.*n0p./ghs_pp.*dghs_pp_dn3i;

d_fch_p = dPhi_dn0i + 1./2.*sigmap.*dPhi_dn1i + pi.*sigmap.^2.*dPhi_dn2i + pi./6.*sigmap.^3.*dPhi_dn3i;

% Chain connectivity for ions
dPhi_dn02 = -N1./Np.*n0p./y_pp.*dy_pp_dn0i(sigmaa2,za2);
dPhi_dn12 = -N1./Np.*n0p./y_pp.*dy_pp_dn1i(za2,sigmaa2);
dPhi_dn22 = -N1./Np.*n0p./y_pp.*dy_pp_dn2i-N2./Np.*n0p./ghs_pp.*dghs_pp_dn2i;
dPhi_dn32 = -N1./Np.*n0p./y_pp.*dy_pp_dn3i(sigmaa2)-N2./Np.*n0p./ghs_pp.*dghs_pp_dn3i;
d_fch_a2 = dPhi_dn02 + 1./2.*sigmaa2.*dPhi_dn12 + pi.*sigmaa2.^2.*dPhi_dn22 + pi./6.*sigmaa2.^3.*dPhi_dn32;



%% ---the electrostatic correlation part of chemical potential ---

fel_p1_f = @(n0i,z,sigma) n0i.*z./(1+Ga.*sigma) ;
fel_p2_f = @(z,sigma) Ga.*z+pi.*Pn.*sigma./2./(1-n3);

dfel_p1_f_dn0i =@(n0i,z,sigma,Ga_sigma,Ga_z) z./(1+Ga.*sigma)-n0i.*z.*sigma.*dGa_dn0i(Ga_sigma,Ga_z)./(1+Ga.*sigma).^2 ;
dfel_p1_f_dn0i_p =@(n0i,z,sigma,Ga_sigma,Ga_z) z./(1+Ga.*sigma)-n0i.*z.*sigma.*dGa_dn0i(Ga_sigma,Ga_z)./(1+Ga.*sigma).^2;
dfel_p1_f_dn0i_p_s =@(n0i,z,sigma,Ga_sigma,Ga_z) -n0i.*z.*sigma.*dGa_dn0i(Ga_sigma,Ga_z)./(1+Ga.*sigma).^2;
dfel_p1_f_dn0i_0inot0i =@(n0i,z,sigma,Ga_sigma,Ga_z) -n0i.*z.*sigma.*dGa_dn0i(Ga_sigma,Ga_z)./(1+Ga.*sigma).^2 ;
dfel_p1_f_dn0i_0inot0i_p =@(n0i,z,sigma,Ga_sigma,Ga_z) -n0i.*z.*sigma.*dGa_dn0i(Ga_sigma,Ga_z)./(1+Ga.*sigma).^2 ;

dfel_p1_f_dn1i_p1_f =@(n0i,z,sigma) -n0i.*z.*sigma./(1+Ga.*sigma).^2;

dfel_p1_f_dn1i =@(n0i,z,sigma,dsigma,dz) dfel_p1_f_dn1i_p1_f(n0i,z,sigma).*dGa_dn1i(dz,dsigma);
dfel_p1_f_dn1i_p =@(n0i,z,sigma,dsigma,dz) dfel_p1_f_dn1i_p1_f(n0i,z,sigma).*dGa_dn1i(dz,dsigma);
dfel_p1_f_dn3i =@(n0i,z,sigma,dsigma) dfel_p1_f_dn1i_p1_f(n0i,z,sigma).*dGa_dn3i(dsigma);
dfel_p1_f_dn3i_p =@(n0i,z,sigma,dsigma) dfel_p1_f_dn1i_p1_f(n0i,z,sigma).*dGa_dn3i(dsigma);
   
dfel_p2_f_dn0i =@(z,sigma,dsigma,dz) z.*dGa_dn0i(dsigma,dz)+pi.*sigma./2./(1-n3).*dPn_dn0i(dsigma,dz);
dfel_p2_f_dn0i_p =@(z,sigma,dsigma,dz) z.*dGa_dn0i(dsigma,dz)+pi.*sigma./2./(1-n3).*dPn_dn0i(dsigma,dz);
dfel_p2_f_dn1i =@(z,sigma,dsigma,dz) z.*dGa_dn1i(dz,dsigma)+pi.*sigma./2./(1-n3).*dPn_dn1i(dz,dsigma);
dfel_p2_f_dn1i_p =@(z,sigma,dsigma,dz) z.*dGa_dn1i(dz,dsigma)+pi.*sigma./2./(1-n3).*dPn_dn1i(dz,dsigma);
dfel_p2_f_dn3i =@(z,sigma,dsigma,dz) z.*dGa_dn3i(dsigma)+pi.*Pn.*sigma./(2.*(1-n3).^2)+pi.*sigma./2./(1-n3).*dPn_dn3i(dsigma);
dfel_p2_f_dn3i_p =@(z,sigma,dsigma,dz) z.*dGa_dn3i(dsigma)+pi.*Pn.*sigma./(2.*(1-n3).^2)+pi.*sigma./2./(1-n3).*dPn_dn3i(dsigma);
          


lamda_n1=lamda_n;

dfel_dn0i_p_1_1 = (-lb.*((dfel_p1_f_dn0i_p(n0p_e,zp_e,sigmap,sigmap,zp_e).*lamda_e+dfel_p1_f_dn0i_p_s(n0p_e,zp_e,sigmap,sigmap,zp_n).*lamda_n).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*(dfel_p2_f_dn0i_p(zp_e,sigmap,sigmap,zp_e).*lamda_e+dfel_p2_f_dn0i_p(zp_e,sigmap,sigmap,zp_n).*lamda_n))) ;%Ga.^2./pi.*dGa_dn0i(sigmap,zp)
       
dfel_dn0i_p_1_3 = (-lb.*((dfel_p1_f_dn0i_0inot0i_p(n0a2,za2,sigmaa2,sigmap,zp_e).*lamda_e+dfel_p1_f_dn0i_0inot0i_p(n0a2,za2,sigmaa2,sigmap,zp_n).*lamda_n).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*(dfel_p2_f_dn0i_p(za2,sigmaa2,sigmap,zp_e).*lamda_e+dfel_p2_f_dn0i_p(za2,sigmaa2,sigmap,zp_n).*lamda_n))) ;   
         

dfel_dn0i_p =  dfel_dn0i_p_1_1  +dfel_dn0i_p_1_3 +Ga.^2./pi.*(dGa_dn0i(sigmap,zp_e).*lamda_e+dGa_dn0i(sigmap,zp_n).*lamda_n1);        

         
         
         

dfel_dn0i_a2_1_2 = (-lb.*(dfel_p1_f_dn0i_0inot0i(n0p_e,zp_e,sigmap,sigmaa2,za2).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*dfel_p2_f_dn0i(zp_e,sigmap,sigmaa2,za2))) ;     
dfel_dn0i_a2_1_3 = (-lb.*(dfel_p1_f_dn0i(n0a2,za2,sigmaa2,sigmaa2,za2).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*dfel_p2_f_dn0i(za2,sigmaa2,sigmaa2,za2))) ;   

dfel_dn0i_a2 =  dfel_dn0i_a2_1_2 +dfel_dn0i_a2_1_3 +Ga.^2./pi.*dGa_dn0i(sigmaa2,za2);   




% dfel_dn1i------         
dfel_dn1i_p_1 = (-lb.*((dfel_p1_f_dn1i_p(n0p_e,zp_e,sigmap,sigmap,zp_e).*lamda_e+dfel_p1_f_dn1i_p(n0p_e,zp_e,sigmap,sigmap,zp_n).*lamda_n).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*(dfel_p2_f_dn1i_p(zp_e,sigmap,sigmap,zp_e).*lamda_e+dfel_p2_f_dn1i_p(zp_e,sigmap,sigmap,zp_n).*lamda_n))) ;%+Ga.^2./pi.*dGa_dn1i_f(zp,sigmap)
         

dfel_dn1i_p_3 = (-lb.*((dfel_p1_f_dn1i_p(n0a2,za2,sigmaa2,sigmap,zp_e).*lamda_e+dfel_p1_f_dn1i_p(n0a2,za2,sigmaa2,sigmap,zp_n).*lamda_n).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*(dfel_p2_f_dn1i_p(za2,sigmaa2,sigmap,zp_e).*lamda_e+dfel_p2_f_dn1i_p(za2,sigmaa2,sigmap,zp_n).*lamda_n)));          
         
    
dfel_dn1i_p =  dfel_dn1i_p_1+dfel_dn1i_p_3+Ga.^2./pi.*(dGa_dn1i(zp_e,sigmap).*lamda_e+dGa_dn1i(zp_n,sigmap).*lamda_n1); 

dfel_dn1i_a2_1 = (-lb.*(dfel_p1_f_dn1i(n0p_e,zp_e,sigmap,sigmaa2,za2).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*dfel_p2_f_dn1i(zp_e,sigmap,sigmaa2,za2))) ;
       
dfel_dn1i_a2_3 = (-lb.*(dfel_p1_f_dn1i(n0a2,za2,sigmaa2,sigmaa2,za2).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*dfel_p2_f_dn1i(za2,sigmaa2,sigmaa2,za2)));          

dfel_dn1i_a2 =  dfel_dn1i_a2_1+dfel_dn1i_a2_3+Ga.^2./pi.*dGa_dn1i(za2,sigmaa2) ; 


% dfel_dn3i------
dfel_dn3i_p_1 = (-lb.*(dfel_p1_f_dn3i_p(n0p_e,zp_e,sigmap,sigmap).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*(dfel_p2_f_dn3i_p(zp_e,sigmap,sigmap,zp_e).*lamda_e+dfel_p2_f_dn3i_p(zp_e,sigmap,sigmap,zp_n).*lamda_n)));     
         

dfel_dn3i_p_3 = (-lb.*(dfel_p1_f_dn3i_p(n0a2,za2,sigmaa2,sigmap).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*(dfel_p2_f_dn3i_p(za2,sigmaa2,sigmap,zp_e).*lamda_e+dfel_p2_f_dn3i_p(za2,sigmaa2,sigmap,zp_n).*lamda_n)));    
   
dfel_dn3i_p=dfel_dn3i_p_1+dfel_dn3i_p_3+Ga.^2./pi.*dGa_dn3i(sigmap) ;     


dfel_dn3i_a2_1 = (-lb.*(dfel_p1_f_dn3i(n0p_e,zp_e,sigmap,sigmaa2).*fel_p2_f(zp_e,sigmap)+ ...
             fel_p1_f(n0p_e,zp_e,sigmap).*dfel_p2_f_dn3i(zp_e,sigmap,sigmaa2,za2)));     
 
dfel_dn3i_a2_3 = (-lb.*(dfel_p1_f_dn3i(n0a2,za2,sigmaa2,sigmaa2).*fel_p2_f(za2,sigmaa2)+ ...
             fel_p1_f(n0a2,za2,sigmaa2).*dfel_p2_f_dn3i(za2,sigmaa2,sigmaa2,za2)));           
  
dfel_dn3i_a2=dfel_dn3i_a2_1+dfel_dn3i_a2_3+Ga.^2./pi.*dGa_dn3i(sigmaa2)  ;    

   

% summary
      
dfel_dn1i_p = 1./2.*sigmap.*dfel_dn1i_p;
dfel_dn2i_p = 0;
dfel_dn3i_p = pi./6.*sigmap.^3.*dfel_dn3i_p ;

      

dfel_dn1i_a2 = 1./2.*sigmaa2.*dfel_dn1i_a2 ;
dfel_dn2i_a2 = 0;
dfel_dn3i_a2 = pi./6.*sigmaa2.^3.*dfel_dn3i_a2 ;

    

d_fel_p = dfel_dn0i_p+dfel_dn1i_p+dfel_dn2i_p+dfel_dn3i_p;

d_fel_a2 = dfel_dn0i_a2+dfel_dn1i_a2+dfel_dn2i_a2+dfel_dn3i_a2;

%% -------the calcium-binding part of chemical potential-----------

n0p = n0p_e+n0p_n;
n3 = n3p_e+n3a2+n3p_n;

x=(n0a2/ma2)./(n0p/mp+n0a2/ma2);
m=ma2*x+mp*(1-x);

dmn02=(mp-ma2)*mp*ma2*n0a2./(n0p*ma2+n0a2*mp).^2;
dmn01=(ma2-mp)*mp*ma2*n0p./(n0p*ma2+n0a2*mp).^2;

a0=Ak(1,1)+(m-1)./m*Ak(1,2)+(m-1).*(m-2)./m.^2*Ak(1,3);
a1=Ak(2,1)+(m-1)./m*Ak(2,2)+(m-1).*(m-2)./m.^2*Ak(2,3);
a2=Ak(3,1)+(m-1)./m*Ak(3,2)+(m-1).*(m-2)./m.^2*Ak(3,3);
a3=Ak(4,1)+(m-1)./m*Ak(4,2)+(m-1).*(m-2)./m.^2*Ak(4,3);
a4=Ak(5,1)+(m-1)./m*Ak(5,2)+(m-1).*(m-2)./m.^2*Ak(5,3);
a5=Ak(6,1)+(m-1)./m*Ak(6,2)+(m-1).*(m-2)./m.^2*Ak(6,3);
a6=Ak(7,1)+(m-1)./m*Ak(7,2)+(m-1).*(m-2)./m.^2*Ak(7,3);

da0=1./m.^2*Ak(1,2)+(3./m.^2-4./m.^3)*Ak(1,3);
da1=1./m.^2*Ak(2,2)+(3./m.^2-4./m.^3)*Ak(2,3);
da2=1./m.^2*Ak(3,2)+(3./m.^2-4./m.^3)*Ak(3,3);
da3=1./m.^2*Ak(4,2)+(3./m.^2-4./m.^3)*Ak(4,3);
da4=1./m.^2*Ak(5,2)+(3./m.^2-4./m.^3)*Ak(5,3);
da5=1./m.^2*Ak(6,2)+(3./m.^2-4./m.^3)*Ak(6,3);
da6=1./m.^2*Ak(7,2)+(3./m.^2-4./m.^3)*Ak(7,3);


I1=a0+a1.*n3+a2.*n3.^2+a3.*n3.^3+a4.*n3.^4+a5.*n3.^5+a6.*n3.^6;

temp=da0+da1.*n3+da2.*n3.^2+da3.*n3.^3+da4.*n3.^4+da5.*n3.^5+da6.*n3.^6;

dI1n01=temp.*dmn01;
dI1n02=temp.*dmn02;

dI1n3=a1+2*a2.*n3+3*a3.*n3.^2+4*a4.*n3.^3+5*a5.*n3.^4+6*a6.*n3.^5;

temp=2*n0p_e.*n0a2*ep_a2/T*sigmap_a2_v^3+n0p_e.^2*epA/T*sigmap_v^3+n0p_n.^2*epM/T*sigmap_v^3;



dfdsI_n02=-4*pi*I1.*(+n0a2*ep_a2/T*sigmap_a2_v^3*lamda_e+n0p_e*(epA/T)*sigmap_v^3*lamda_e+n0p_n*(epM/T)*sigmap_v^3*lamda_n)-2*pi*dI1n02.*temp;
dfdsI_n01=-4*pi*I1.*(n0p_e*ep_a2/T*sigmap_a2_v^3)-2*pi*dI1n01.*temp;

dfdsI_n3=-2*pi*dI1n3.*temp;

dfdsI2=dfdsI_n02+pi/6*sigmap^3*dfdsI_n3;
dfdsI1=dfdsI_n01+pi/6*sigmaa2^3*dfdsI_n3;

b0=Bk(1,1)+(m-1)./m*Bk(1,2)+(m-1).*(m-2)./m.^2*Bk(1,3);
b1=Bk(2,1)+(m-1)./m*Bk(2,2)+(m-1).*(m-2)./m.^2*Bk(2,3);
b2=Bk(3,1)+(m-1)./m*Bk(3,2)+(m-1).*(m-2)./m.^2*Bk(3,3);
b3=Bk(4,1)+(m-1)./m*Bk(4,2)+(m-1).*(m-2)./m.^2*Bk(4,3);
b4=Bk(5,1)+(m-1)./m*Bk(5,2)+(m-1).*(m-2)./m.^2*Bk(5,3);
b5=Bk(6,1)+(m-1)./m*Bk(6,2)+(m-1).*(m-2)./m.^2*Bk(6,3);
b6=Bk(7,1)+(m-1)./m*Bk(7,2)+(m-1).*(m-2)./m.^2*Bk(7,3);

db0=1./m.^2*Bk(1,2)+(3./m.^2-4./m.^3)*Bk(1,3);
db1=1./m.^2*Bk(2,2)+(3./m.^2-4./m.^3)*Bk(2,3);
db2=1./m.^2*Bk(3,2)+(3./m.^2-4./m.^3)*Bk(3,3);
db3=1./m.^2*Bk(4,2)+(3./m.^2-4./m.^3)*Bk(4,3);
db4=1./m.^2*Bk(5,2)+(3./m.^2-4./m.^3)*Bk(5,3);
db5=1./m.^2*Bk(6,2)+(3./m.^2-4./m.^3)*Bk(6,3);
db6=1./m.^2*Bk(7,2)+(3./m.^2-4./m.^3)*Bk(7,3);


M=1+m.*(8*n3-2*n3.^2)./(1-n3).^4+(1-m).*(20*n3-27*n3.^2+12*n3.^3-2*n3.^4)./(1-n3).^2./(2-n3).^2;


temp=(8*n3-2*n3.^2)./(1-n3).^4-(20*n3-27*n3.^2+12*n3.^3-2*n3.^4)./(1-n3).^2./(2-n3).^2;

dMn01=temp.*dmn01;
dMn02=temp.*dmn02;

dMn3=m.*(8+20*n3-4*n3.^2)./(1-n3).^5+(1-m).*(2*n3.^3+12*n3.^2-48*n3+40)./(1-n3).^3./(2-n3).^3;

I2=b0+b1.*n3+b2.*n3.^2+b3.*n3.^3+b4.*n3.^4+b5.*n3.^5+b6.*n3.^6;


temp=db0+db1.*n3+db2.*n3.^2+db3.*n3.^3+db4.*n3.^4+db5.*n3.^5+db6.*n3.^6;

dI2n01=temp.*dmn01;
dI2n02=temp.*dmn02;

dI2n3=b1+2*b2.*n3+3*b3.*n3.^2+4*b4.*n3.^3+5*b5.*n3.^4+6*b6.*n3.^5;

temp=2*n0p_e.*n0a2*(ep_a2/T)^2*sigmap_a2_v^3+n0p_e.^2*(epA/T)^2*sigmap_v^3+n0p_n.^2*(epM/T)^2*sigmap_v^3;

dfdsII_n02=-2*pi./M.*I2.*m.*(n0a2.*(ep_a2/T)^2*sigmap_a2_v^3*lamda_e+n0p_e*(epA/T)^2*sigmap_v^3*lamda_e+n0p_n*(epM/T)^2*sigmap_v^3*lamda_n)...
    +pi./M.^2.*dMn02.*I2.*m.*temp...
    -pi./M.*dI2n02.*m.*temp...
    -pi./M.*I2.*dmn02.*temp;

dfdsII_n01=-2*pi./M.*I2.*m.*(n0p_e*(ep_a2/T)^2*sigmap_a2_v^3)...
    +pi./M.^2.*dMn01.*I2.*m.*temp...
    -pi./M.*dI2n01.*m.*temp...
    -pi./M.*I2.*dmn01.*temp;

dfdsII_n3=pi./M.^2.*dMn3.*I2.*m.*temp...
    -pi./M.*dI2n3.*m.*temp;


dfdsII2=dfdsII_n02+pi/6*sigmap^3*dfdsII_n3;
dfdsII1=dfdsII_n01+pi/6*sigmaa2^3*dfdsII_n3;

%% ------------------------final part----------------------------------------
amu_p = d_fid_p+d_fhs_p+d_fch_p+d_fel_p+(dfdsI2+dfdsII2);
amu_a2 = d_fid_a2+d_fhs_a2+d_fch_a2+d_fel_a2+dfdsI1+dfdsII1;


switch a
    case 'a2'
        amu=amu_a2;
    case 'p'
        amu=amu_p;
    case 'all'
        amu=[amu_p',amu_a2'];
end

end