colors
ifig1=7;
ifig2=8;
ifig3=9;
ifig4=1;
ifig5=11;
ifig6=12;
figure(ifig1),clf
figure(ifig2),clf
figure(ifig3),clf
figure(ifig4),clf
figure(ifig5),clf
figure(ifig6),clf

dataxisq100
h=pi/100;
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/200;
dataxisq200
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/400;
dataxisq400
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/800;
dataxisq800
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/1600;
dataxisq1600
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

vex=v; udex=ud; vddex=vdd;
dataxisq100
n=100;
h=pi/n;
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5); 
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on
%disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud),max(errvdd)))
disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud)*3*h,max(errvdd)*9*h^2))

n=200;
h=pi/n;
dataxisq200
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on
%disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud),max(errvdd)))
disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud)*3*h,max(errvdd)*9*h^2))

n=400;
h=pi/n;
dataxisq400
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on
disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud)*3*h,max(errvdd)*9*h^2))

n=800;
h=pi/n;
dataxisq800
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on
%disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud),max(errvdd)))
disp(sprintf('%d %0.3g %0.3g %0.3g',n,max(errv),max(errud)*3*h,max(errvdd)*9*h^2))


isize=17;
figure(ifig1), xlabel('z'), ylabel('v_Q'), set(gca,'fontsize',isize)
%print -dpdf fig7
pause
figure(ifig2), xlabel('z'), ylabel('ud_Q'), set(gca,'fontsize',isize)
%print -dpdf fig8
pause
figure(ifig3), xlabel('z'), ylabel('vdd_Q'), set(gca,'fontsize',isize)
%print -dpdf fig9
pause
figure(ifig4), xlabel('z'), ylabel('E[v_Q]'), set(gca,'fontsize',isize)
%print -dpdf fig10
pause
figure(ifig5), xlabel('z'), ylabel('E[ud_Q]'), set(gca,'fontsize',isize)
%print -dpdf fig11
pause
figure(ifig6), xlabel('z'), ylabel('E[vdd_Q]'), set(gca,'fontsize',isize)
%print -dpdf fig12
