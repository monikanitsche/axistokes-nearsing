colors
ifig1=5;
ifig2=6;
ifig3=7;
ifig4=8;
ifig5=9;
ifig6=10;
figure(ifig1),clf
figure(ifig2),clf
figure(ifig3),clf
figure(ifig4),clf
figure(ifig5),clf
figure(ifig6),clf

dataxism100v
h=pi/100;
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/200;
dataxism200v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/400;
dataxism400v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/800;
dataxism800v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

h=pi/1600;
dataxism1600v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
figure(ifig1), plot(z,v,'linewidth',2),hold on
figure(ifig2), plot(z,ud,'linewidth',2),hold on
figure(ifig3), plot(z,vdd,'linewidth',2),hold on

vex=v; udex=ud; vddex=vdd;
dataxism100v
h=pi/100;
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5); 
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on

h=pi/200;
dataxism200v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on

h=pi/400;
dataxism400v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on

h=pi/800;
dataxism800v
z=a(:,2); v=a(:,3); ud=a(:,4); vdd=a(:,5);
errv=abs(v-vex); errud=abs(ud-udex); errvdd=abs(vdd-vddex);
figure(ifig4), plot(z,errv,'linewidth',2),hold on
figure(ifig5), plot(z,errud,'linewidth',2),hold on
figure(ifig6), plot(z,errvdd,'linewidth',2),hold on

isize=17;
figure(ifig1), xlabel('z'), ylabel('v_M'), set(gca,'fontsize',isize)
%print -dpdf fig1, pause
figure(ifig2), xlabel('z'), ylabel('ud_M'), set(gca,'fontsize',isize)
%print -dpdf fig2, pause
figure(ifig3), xlabel('z'), ylabel('vdd_M'), set(gca,'fontsize',isize)
%print -dpdf fig3, pause
figure(ifig4), xlabel('z'), ylabel('E[v_M]'), set(gca,'fontsize',isize)
%print -dpdf fig4, pause
figure(ifig5), xlabel('z'), ylabel('E[ud_M]'), set(gca,'fontsize',isize)
%print -dpdf fig5, pause
figure(ifig6), xlabel('z'), ylabel('E[vdd_M]'), set(gca,'fontsize',isize)
%print -dpdf fig6
