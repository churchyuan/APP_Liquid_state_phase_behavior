function www=balance2_s(lb1,lb2,deltalb)
global lb;

seewrong=0;

% If seewrong is equal to 1, solution figure is obtained. It shows what
% the error is. Also, the region of rhop1 can be extended to make sure that
% result exists.

www=zeros(round(lb2-lb1)+1,4);
zz=1;
for lb=lb1:-deltalb:lb2
    if lb ==0
        lb=1e-6;
    end
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
    step = 2000;

    rhop1 = logspace(-200,-15,step);
    rhop2 = logspace(-15,-10,step);
    rhop3 = logspace(-10,-4,step);
    rhop4 = logspace(-4,-2.55,step);
    rhop5 = logspace(-2.55,-0.05,step);

    rhoa1 = rhop1*zp_c/za2;
    rhoa2 = rhop2*zp_c/za2;
    rhoa3 = rhop3*zp_c/za2;
    rhoa4 = rhop4*zp_c/za2;
    rhoa5 = rhop5*zp_c/za2;
    rhop5(rhoa5>0.87)=[];
    rhoa5(rhoa5>0.87)=[];




    %% drawing the solution diagram
    if(seewrong==1)
        plot(amu2_s(rhop1,rhoa1,'p')*za2+amu2_s(rhop1,rhoa1,'a2')*zp_c,p_final2_s_noT(rhop1,rhoa1,0),'y.')
        hold on
        plot(amu2_s(rhop2,rhoa2,'p')*za2+amu2_s(rhop2,rhoa2,'a2')*zp_c,p_final2_s_noT(rhop2,rhoa2,0),'g.')
        plot(amu2_s(rhop3,rhoa3,'p')*za2+amu2_s(rhop3,rhoa3,'a2')*zp_c,p_final2_s_noT(rhop3,rhoa3,0),'b.')
        plot(amu2_s(rhop4,rhoa4,'p')*za2+amu2_s(rhop4,rhoa4,'a2')*zp_c,p_final2_s_noT(rhop4,rhoa4,0),'k.')
        plot(amu2_s(rhop5,rhoa5,'p')*za2+amu2_s(rhop5,rhoa5,'a2')*zp_c,p_final2_s_noT(rhop5,rhoa5,0),'r.')
    end
    %% main calculation
    rhop=[rhop1,rhop2,rhop3,rhop4,rhop5];
    rhoa=[rhoa1,rhoa2,rhoa3,rhoa4,rhoa5];
    x=amu2_s(rhop,rhoa,'p')*za2+amu2_s(rhop,rhoa,'a2')*zp_c;
    y=p_final2_s_noT(rhop,rhoa,0);
    A=calculation_c(x,y,rhop);
    if(lb<1e-5 )
        disp('the lb is reach its mimimum (1e-4)')
        if (length(A)==3)
            options=optimset('Display','none','TolFun',.0000001,'MaxIter',5000,'MaxFunEvals',5000);
            [xsolu,err]=fsolve(@eq_b2,A(1:2),options);
            www(zz,1:4)=[xsolu,A(3),sum(abs(err))];
            eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);
            return
        else
        eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);
            return
        end
    end
    if(A==-2)
        error('There is no result, please extend the region of rhop1');
    end
    if (A==-1 & isempty(www))
        return;
    end
    if( A==-1 & abs(www(zz-1,1)-www(zz-1,2))<www(zz-1,2)*0.01 )
        disp(['The reached maximum: lb = ',num2str(lb),'σ.']);
        disp('The code is finished.');
        www=www(1:zz-1,:);
        eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);
        return
    elseif(A==-1 & abs(www(zz-1,1)-www(zz-1,2))>1e-5 )
        www2=balance2_s(lb+deltalb,lb2,deltalb/2);
        www=[www(1:zz-1,:);www2];
        eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);
        return
    end

    options=optimset('Display','none','TolFun',.0000001,'MaxIter',5000,'MaxFunEvals',5000);
    [xsolu,err]=fsolve(@eq_b2,A(1:2),options);
    www(zz,1:4)=[xsolu,A(3),sum(abs(err))];

    zz=zz+1;
end
eval(['save(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',',''www''',');']);
end




    






