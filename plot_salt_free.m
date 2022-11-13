% Function that plots a phase diagram of salt-free solution
% cp is the critical point (as lamda_e epM rhopc and lbc).
% cp(j1,:) = [lamda_e,epM,c]; can be rewritten for different values.
% cp2 is the data of rhop and lb in phase diagram.
clear cp cp2;
j1=1;
dbstop if error
% lamda_ei=linspace(0.00001,1,10);
for zp_e =[-1] %  valence for polymer charged part
    for za2 =[1:1:10] %  valence for counterion
        for Np=1e4 % total length of the chain
            for lamda_e=[1]% charged fraction[0.1:0.05:0.5] [0.6:0.1:1]
                for N1=ceil(Np*lamda_e*1)-1 % number of the charged chemical bonding
                    for ep=[0.3] % epsilon_AC - strength of polymer charged-counterion interaction
                        for epM=0.2 % epsilon_B - strength of polymer neutral-neutral interaction
                            for epA = 0.2 % epsilon_A - strength of polymer charged-charged interaction
                                [c,c2] = draw(Np,N1,lamda_e,ep,epM,epA,zp_e,za2);
                                cp(j1,:) = [za2,c];
                                cp2(1:size(c2,1),j1*2-1:j1*2)=c2;
                                j1 = j1+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

xlabel('\rho_{polymer}\sigma^3');
ylabel('lb/\sigma');
title('influence of N');
function [c,c2] = draw(Np,N1,lamda_e,ep,epM,epA,zp_e,za2)
if exist(['Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat'],'file')
eval(['load(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',');']);

data=www;
p1=data(1:end,1);
p3=data(1:end,2);
lb=data(1:end,3);

hold on
plot([p1;flipud(p3)],[lb;flipud(lb)],'-','Linewidth',2)
plot((p1(end)+p3(end))/2,lb(end),'.','MarkerSize',18,'color',[1 0.5 0]);
c = [(p1(end)+p3(end))/2,lb(end)];
c2 = [[p1;flipud(p3)],[lb;flipud(lb)]];
xlabel('\rho_{polymer}\sigma^3');
ylabel('lb/\sigma');
title('Phase Diagram');
else
    c=[0 8];
    c2=[0 0];
end
end