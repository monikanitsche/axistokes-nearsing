clear
ifig1= 4;
ifig2= 5;

isize=18;
offblue=[0 0.4470 0.7410]; %	'#0072BD'	 off blue
offred=[0.8500 0.3250 0.0980]; %	'#D95319'     off red	
orangy=[0.9290 0.6940 0.1250]; %	'#EDB120'	orangy
purple=[0.4940 0.1840 0.5560]; %	'#7E2F8E'	ipurple
offgreen=[0.4660 0.6740 0.1880]; %	'#77AC30' green	
dkcyan=[0.3010 0.7450 0.9330]; %	'#4DBEEE'dark cyan	
dkbrown=[0.6350 0.0780 0.1840]; %	'#A2142F'	dark brown/red


%----------------Emax[Ms] O(h^3)------------------------------
% DATA in dat..3r.m (tst.m)

% n   best         h          2h         3h        GLOBAL
%     x>.6        x>h        x>2h       x>3h       all x
a=[
 100  1.96e-07   5.08e-07   5.08e-07   5.08e-07   5.08e-07
 200  2.43e-08   6.29e-08   6.29e-08   6.29e-08   6.29e-08
 400  2.97e-09   7.60e-09   7.60e-09   7.60e-09   7.60e-09
 800  3.51e-10   8.58e-10   8.58e-10   8.58e-10   8.58e-10
];
n=a(:,1); 
h=1./n;
e3rbest=a(:,2);
e3r=a(:,6);
%----------------Emax[Ms] O(h^4)------------------------------
% DATA in dat..4r.m (tst.m)

% n   best         h          2h         3h        GLOBAL
%     x>.6        x>h        x>2h       x>3h       all x
a=[
100  2.80e-09   1.10e-08   1.10e-08   1.10e-08   3.56e-05
200  1.69e-10   1.90e-09   7.03e-10   7.03e-10   6.54e-07
400   1.02e-11   9.55e-10   1.42e-10   4.29e-11   1.06e-08
800   6.61e-13   8.50e-11   8.50e-11   3.15e-12   3.36e-10
];
e4rbest=a(:,2);
e4r3h=a(:,5);

% n   best       top0       top3      top4     top5    top6
%     x>.6       x>.25      x>.05     x>.04    x>.03     x>0.02
a=[
100   2.80e-09   9.49e-09   1.10e-08   1.10e-08   1.10e-08   1.58e-08
200   1.69e-10   6.00e-10   7.03e-10   7.03e-10   7.03e-10   7.03e-10
%400   1.02e-11   3.66e-11   4.29e-11   4.29e-11   4.29e-11   4.29e-11
 400   1.04e-11   3.77e-11   4.46e-11   4.46e-11   4.46e-11   4.46e-11
800   6.61e-13   2.29e-12   2.67e-12   2.67e-12   2.67e-12   2.67e-12
];
e4top5=a(:,6);

%----------------Emax[Ms] O(h^3)------------------------------
figure(ifig1),clf
loglog(h,e3r,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
%loglog(h,corpl,'.-.','linewidth',1,'markersize',20,'color',offblue)
hh=[.0015,0.009]; loglog(hh,0.27*hh.^3,'--','linewidth',2,'markersize',20,'color',offred)
%loglog(h,e4rh,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
loglog(h,e4top5,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
%loglog(h,e4rbest,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
hh=[.0015,0.009]; loglog(hh,0.6*hh.^4,'--','linewidth',2,'markersize',20,'color',offred)
text(0.005,2.3e-8,'$\sim h^3$','fontsize',isize+2,'Interpreter','latex')
text(0.005,2.3e-10,'$\sim h^4$','fontsize',isize+2,'Interpreter','latex')
%legend('$||\mathbf{x}_0-\mathbf{x}_{ax}||\ge 3h$','all gridpoints $\mathbf{x}_0$','location','southeast','Interpreter','latex','fontsize',isize+1)
set(gca,'fontsize',18) %, xlabel('h')
ylabel('$\max||\mathbf{u^s}-\mathbf{u^s}_{ex}||$','Interpreter','latex','fontsize',isize+2)
xlabel('$1/n$','Interpreter','latex','fontsize',isize+2)
ylim([2.e-12 2.e-6])
set(gca,'YTick',10.^[-11:-6])
text(0.00135,0.9e-6,'(a)','fontsize',isize+2)
axis square

%----------------Emax[Qs] O(h^3)------------------------------
% DATA in dat..3r.m
figure(ifig2),clf

h=1./[100 200 400 800];
e3rbest =[7.8e-7 9.7e-8 1.2e-8 1.5e-9];  % 7.8e-9
e3qbest =[7.8e-7 9.7e-8 1.2e-8 1.5e-9];  

e3r3h   =[1.2e-6 1.5e-7 1.9e-8 2.1e-9];  % 7.8e-9
e3q3h   =[1.2e-6 1.5e-7 1.9e-8 2.1e-9];  

%dble: using 255,233,1 in r3.m no fix
%first 5 columns, before extrapolation from axis
%Global errors confirmed for all "YES"
% n   best         h          2h         3h        GLOBAL
%     x>.6        x>h        x>2h       x>3h       all x
a=[ 
100  7.85e-07   1.21e-06   1.21e-06   1.21e-06   1.21e-06  %YES
200  9.67e-08   4.53e-07   1.55e-07   1.55e-07   1.55e-07  %YES
400  1.20e-08   4.17e-07   2.77e-08   1.90e-08   1.90e-08  %YES
800  1.48e-09   3.98e-08   3.98e-08   2.09e-09   2.09e-09  %YES
];
e3rbest=a(:,2);
e3r=a(:,6);


%----------------Emax[Qs] O(h^4)------------------------------
%dble: using 366,344,2 in r4.m no fix
% n   best       top0       top3      top4     top5
%     x>.6       x>.25      x>.05     x>.04    x>.03
a=[
100 1.05e-08   5.30e-08   6.58e-08   1.28e-07   5.23e-07 %same as q
200 6.28e-10   2.99e-09   3.40e-09   7.92e-09   2.74e-08 %same as q
400 4.01e-11   1.81e-10   2.05e-10   2.05e-10   6.84e-10 %same as q
%800 1.24e-11   1.24e-11   1.24e-11   1.24e-11   1.24e-11 %slightly larger
 800 9.41e-12   1.02e-11   1.10e-11   1.10e-11  1.14e-11
];
er4best=a(:,2);
er4top0=a(:,3);
er4top3=a(:,4);

%quad: using 366,344,2 in q4.m no fix
%     x>.6       x>.25      x>.05     x>.04    x>.03
a=[
100  1.05e-08   5.30e-08   6.58e-08   1.28e-07   5.23e-07  
200  6.28e-10   2.99e-09   3.40e-09   7.92e-09   2.74e-08 
400  4.02e-11   1.81e-10   2.05e-10   2.05e-10   6.84e-10
800  2.43e-12   1.03e-11   1.10e-11   1.10e-11   1.14e-11
];
eq4best=a(:,2);
eq4top0=a(:,3);
eq4top3=a(:,4);

%---------------- Q plots ------------------------------
loglog(h,e3r,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
%loglog(h,corpl,'.-.','linewidth',1,'markersize',20,'color',offblue)
hh=[.0015,0.009]; loglog(hh,0.7*hh.^3,'--','linewidth',2,'markersize',20,'color',offred)
loglog(h,er4top3,'.-','linewidth',3,'markersize',20,'color',offblue),hold on
hh=[.0015,0.009]; loglog(hh,3*hh.^4,'--','linewidth',2,'markersize',20,'color',offred)
text(0.005,5.5e-8,'$\sim h^3$','fontsize',isize+2,'Interpreter','latex')
text(0.005,11.e-10,'$\sim h^4$','fontsize',isize+2,'Interpreter','latex')
%legend('$||\mathbf{x}_0-\mathbf{x}_{ax}||\ge 3h$','all gridpoints $\mathbf{x}_0$','location','southeast','Interpreter','latex','fontsize',isize+1)
set(gca,'fontsize',18) %, xlabel('h')
ylabel('$\max||\mathbf{u^d}-\mathbf{u^d}_{ex}||$','Interpreter','latex','fontsize',isize+2)
xlabel('$1/n$','Interpreter','latex','fontsize',isize+2)
ylim([2.e-12 2.e-6])
set(gca,'YTick',10.^[-11:-6])
text(0.00135,0.9e-6,'(b)','fontsize',isize+2)
axis square

figure(ifig1),print -depsc Fig07a-maxerrsgle
pause
figure(ifig2),print -depsc Fig07b-maxerrdble


