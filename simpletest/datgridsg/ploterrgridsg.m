%plots absolute errors in single layer velocity computed 
% 1) without any corrections:             T_6[G] 
% 2) with correction near sheet:          T_6[G] + E[H]
% 3) with correction near sheet and axis: T_6[G] + E[H] + E[B]
%Figs for paper:
%  errMs4th, errMs4thcorr, errMs6thcorr 

clear
isize=19;
ifig2=1; ifig3=2; ifig4=3; 
axisvec=[-1.2,0,0,1.1];

colvecs=-10:0.1:-2; colvecd=-14:0.1:1;
colormap('default')
% default=parula. also jet, pink, jet, hot, cool, spring, summer, autumn, winter


%figure(ifig1),clf 
figure(ifig2),clf 
figure(ifig3),clf 
figure(ifig4),clf 
colvec=colvecd;

%"Exact" data
load dats1600r4;

z1ex=z1; z2ex=z2; z3ex=z3; z4ex=z4; z5ex=z5; z6ex=z6; 
errmax=0;
errmaxpre=0;
errmaxtop=0;
errmaxnocor=0;
errmaxtop2=0;
errmaxtop3=0;

%FOURTH ORDER DATA %Real, No Fix 
cax=[-14,-6];
minlev=1.e-14;
%np=100; load dats100r4;
%np=200; load dats200r4;
%np=400; load dats400r4;
np=800; load dats800r4;

%THIRD ORDER DATA %Real, no Fix
%cax=[-10,-2];
%minlev=1.e-11;
%np=100; load dats100r3;
%np=200; load dats200r3;
%np=400; load dats400r3;
%np=800; load dats800r3;

a=z1; b=z1ex;  ng=71;   mg=length(a)/ng; ploterrgridsuball
a=z2; b=z2ex;  ng=151;  mg=length(a)/ng; ploterrgridsuball
a=z3; b=z3ex;  ng=501;  mg=length(a)/ng; ploterrgridsuball
a=z4; b=z4ex;  ng=1001; mg=length(a)/ng; ploterrgridsuball
a=z4; b=z4ex;  ng=1001; mg=length(a)/ng; ploterrgridsuball3
a=z5; b=z5ex;  ng=501;  mg=length(a)/ng; ploterrgridsuball
a=z6; b=z6ex;  ng=151;  mg=length(a)/ng; ploterrgridsuball
%a=z7; b=z7ex;  ng=11;   mg=length(a)/ng; ploterrgridsuball


figure(ifig4)
np
alf=linspace(0,pi,np+1); x=0.95*sin(alf); y=-0.95*cos(alf);
%plot(y,x,'.k','markersize',15)
%figure(ifig1), axis square, axis equal, axis(axisvec)
%set(gca,'fontsize',isize)
figure(ifig2), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)
figure(ifig3), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)
figure(ifig4), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)

%plotfin

%ploterrfinishm

%figure(1), print -depsc errqn100nocorr  %no correction
%figure(2), print -depsc errqn100onecorr  %strip
%figure(3), print -depsc errqn100twocorr  %strip+axis
%figure(4), print -depsc errqn100threecorr  %strip+axis+interp
%figure(5), print -depsc errqn200nocorr  %no correction
%figure(6), print -depsc errqn200onecorr  %strip
%figure(7), print -depsc errqn200twocorr  %strip+axis
%figure(8), print -depsc errqn200threecorr  %strip+axis+interp
%figure(9),  print -depsc errqn400nocorr  %no correction
%figure(10), print -depsc errqn400onecorr  %strip
%figure(11), print -depsc errqn400twocorr  %strip+axis
%figure(12), print -depsc errqn400threecorr  %strip+axis+interp
%figure(13), print -depsc errqn800nocorr  %no correction
%figure(14), print -depsc errqn800onecorr  %strip
%figure(15), print -depsc errqn800twocorr  %strip+axis
%figure(16), print -depsc errqn800threecorr  %strip+axis+interp
