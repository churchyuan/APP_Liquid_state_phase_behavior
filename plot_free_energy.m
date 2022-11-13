function www=plot_free_energy(lb1,lb2,deltalb)
global lb;

    PEsolution = INPUT;
    lamda_e=PEsolution.lamda_e ;
    N1=PEsolution.N1 ;
    Np=PEsolution.Np ;
    e_t=PEsolution.eee;
    epM = PEsolution.eppM;
    epA = PEsolution.eppA;
    zp_e = PEsolution.zp_e;
    ep = PEsolution.ep;
    zp_c = abs(PEsolution.zp_c);
    za2 = PEsolution.za2;

seewrong=0;

% If seewrong is equal to 1, solution figure is obtained. It shows what
% the error is. Also, the region of rhop1 can be extended to make sure that
% result exists.

www=zeros(round(lb2-lb1)+1,4);
zz=1;
for lb=lb1:-deltalb:lb2
%   rhop1 = logspace(-20,0,10000);
    rhop1=[1.00577595930163e-15	0.259108548154065]; % use paramter show the free energy line
  rhoa1 = rhop1*zp_c/za2;
    if lb ==0
        lb=1e-6;
    end
    f = free_energy2_noT(rhop1,rhoa1);
semilogx(rhop1,f)
hold on

%     eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);

end

xlabel('\rho_p \sigma^3')
ylabel('free energy/ k_BT')
legend('lb =7\sigma','lb =6','lb =5','lb =4','lb =3','lb =2','lb =1','lb =0')