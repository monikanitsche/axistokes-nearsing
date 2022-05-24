clear
colors
ifig1=5;
ifig2=6;
figure(ifig1),clf
figure(ifig2),clf

dataxism100sh
h=pi/100;
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/200;
dataxism200sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/400;
dataxism400sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/800;
dataxism800sh
z=a(:,2); v=a(:,3);
figure(ifig1), plot(z,v,'linewidth',2),hold on

h=pi/1600;
dataxism1600sh
z=a(:,2); v=a(:,3); 
figure(ifig1), plot(z,v,'linewidth',2),hold on

vex=v; 
dataxism100sh
n=100;
h=pi/n;
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))
%[emax,imax]=max(errv)

n=200;
h=pi/n;
dataxism200sh
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

n=400;
h=pi/n;
dataxism400sh
z=a(:,2); v=a(:,3);
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

n=800;
h=pi/n;
dataxism800sh
z=a(:,2); v=a(:,3); 
errv=abs(v-vex); 
figure(ifig2), plot(z,errv,'linewidth',2),hold on
disp(sprintf('%d %0.3g ',n,max(errv) ))

isize=17;
figure(ifig1), xlabel('z'), ylabel('v_Q'), set(gca,'fontsize',isize)
figure(ifig2), xlabel('z'), ylabel('E[v_Q]'), set(gca,'fontsize',isize)
