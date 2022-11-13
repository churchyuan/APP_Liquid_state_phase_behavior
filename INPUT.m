function PEsolution = INPUT
% Input data of the system
% sigma is the size (diameter) of a atom and the unit is Å(0.1nm).
% z means the electricity and p is polymer, a is positive atom and b is negative atom.
% There are two cations: H+ (a1) and Ca2+ (a2).
global  lamda_e N1 Np ep epM epA lb za2 zp_e
k_B = 1.38064852*1e-23;
% %-----------ion radii
PEsolution.sigma = 4.0; % Reference size
PEsolution.sigmap = 4.0/PEsolution.sigma; 
PEsolution.sigmaa2 = 4.0/PEsolution.sigma;

% Free energy of electrostatic correlation 
PEsolution.zp_n = 0;
PEsolution.zp_e = zp_e;
PEsolution.za2 = za2;

% Temperature
PEsolution.T = 298.15; 

% Bjerrum length
PEsolution.lb=lb;

PEsolution.lamda_e =lamda_e ;
PEsolution.lamda_n = 1-PEsolution.lamda_e;
PEsolution.Np = Np;
PEsolution.N1= N1;
PEsolution.N2 = Np-N1-1;

%-----------------------------------------
PEsolution.zp_c = -PEsolution.zp_e*PEsolution.lamda_e ;
% Number of the polymer's connection
PEsolution.Na = 1;
PEsolution.Na2 = 1;
PEsolution.Nb = 1;

PEsolution.transform = k_B*PEsolution.T/((PEsolution.sigma*1e-10)^3)*1e-6;

% Parametes for CBI
PEsolution.Ak = [0.9105631445 -0.3084016918 -0.0906148351 0.7240946941 -0.5755498075 0.0976883116
    0.6361281449 0.1860531159 0.4527842806 2.2382791861 0.6995095521 -0.2557574982
    2.6861347891 -2.5030047259 0.5962700728 -4.0025849485 3.8925673390 -9.1558561530
    -26.547362491 21.419793629 -1.7241829131 -21.003576815 -17.215471648 20.642075974
    97.759208784 -65.255885330 -4.1302112531 26.855641363 192.67226447 -38.804430052
    -159.59154087 83.318680481 13.776631870 206.55133841 -161.82646165 93.626774077
    91.297774084 -33.746922930 -8.6728470368 -355.60235612 -165.20769346 -29.666905585];
PEsolution.Bk = PEsolution.Ak(:,4:6);
PEsolution.Ak = PEsolution.Ak(:,1:3);
PEsolution.ep=ep;
 PEsolution.eppM = epM;
PEsolution.eppA = epA ;
e_t=PEsolution.ep*PEsolution.T ;
PEsolution.eee =PEsolution.ep;
PEsolution.eca = e_t;
PEsolution.epce_ca = PEsolution.eca;
M=1;
PEsolution.M=M;
PEsolution.MCa = M;
PEsolution.MPCE =M*PEsolution.Np;
% Functions for invoking the soft-sphere model
PEsolution.sigmap_v=PEsolution.sigmap;
PEsolution.sigmaa2_v=PEsolution.sigmaa2;
PEsolution.sigmapce_ca_v = (PEsolution.sigmaa2_v+PEsolution.sigmap_v)/2;
PEsolution.sigmap_av=round(PEsolution.sigmap_v,1);

end