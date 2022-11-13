function draw(Np,N1,lamda_e,ep,epM,epA,zp_e,za2)

 eval(['load(''', 'Np',num2str(Np),'N1_',num2str(N1),'lamda_e',num2str(lamda_e),'_epsilon_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'zp_',num2str(zp_e),'zs1_',num2str(za2),'.mat''',');']);
 
data=www;
p1=data(:,1);
p3=data(:,2);
lb=data(:,3);
hold on
plot((p1(end)+p3(end))/2,lb(end),'.','MarkerSize',18,'color',[1 0.5 0]);
plot(p1,lb);
plot(p3,lb);

end