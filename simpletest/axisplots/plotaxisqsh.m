clear
colors
ifig1=5;
ifig2=6;
figure(ifig1),clf
figure(ifig2),clf

dataxisq100sh
h=pi/100;
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/200;
dataxisq200sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/400;
dataxisq400sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/800;
dataxisq800sh
z=a(:,2); v=a(:,3);
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/1600;
dataxisq1600sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

vex=v; 
dataxisq100sh
n=100;
h=pi/n;
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))
%[emax,imax]=max(errv)

n=200;
h=pi/n;
dataxisq200sh
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

n=400;
h=pi/n;
dataxisq400sh
z=a(:,2); v=a(:,3);
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

n=800;
h=pi/n;
dataxisq800sh
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

isize=17;
figure(ifig1), xlabel('z'), ylabel('v_Q'), set(gca,'fontsize',isize)
figure(ifig2), xlabel('z'), ylabel('E[v_Q]'), set(gca,'fontsize',isize)
