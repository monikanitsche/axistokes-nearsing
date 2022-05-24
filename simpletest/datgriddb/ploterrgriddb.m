%plots absolute errors in dble layer

clear
isize=19;
ifig2=1; ifig3=2; ifig4=3; 
axisvec=[-1.2,0,0,1.1];

colvecs=-10:0.1:-2; colvecd=-14:0.1:1;
colormap('default')
% default=parula. also jet, pink, jet, hot, cool, spring, summer, autumn, winter

figure(ifig2),clf 
figure(ifig3),clf 
figure(ifig4),clf 
colvec=colvecd;

load datd1600r4;
z1ex=z1; z2ex=z2; z3ex=z3; z4ex=z4; z5ex=z5; z6ex=z6; 
%z4ex=z4;

%FOURTH ORDER DATA
cax=[-14,-6];
minlev=1.e-14;

%Real, No Fix 
%np=100; load datd100r4;
%np=200; load datd200r4;
%np=400; load datd400r4;
np=800; load datd800r4;

%Quad, No Fix 
%np=100; load datd100q4;
%np=200; load datd200q4;
%np=400; load datd400q4;
%np=800; load datd800q4;

%THIRD ORDER DATA
%cax=[-10,1];
%minlev=1.e-11;
%Real, no Fix
%np=100; load datd100r3;
%np=200; load datd200r3;
%np=400; load datd400r3;
%np=800; load datd800r3;

%Quad, no Fix
%np=100; load datd100q3;
%np=200; load datd200q3;
%np=400; load datd400q3;
%np=800; load datd800q3;

a=z1; b=z1ex;  ng=71;   mg=length(a)/ng; ploterrgridsuball
a=z2; b=z2ex;  ng=151;  mg=length(a)/ng; ploterrgridsuball
a=z3; b=z3ex;  ng=501;  mg=length(a)/ng; ploterrgridsuball
a=z4; b=z4ex;  ng=1001; mg=length(a)/ng; ploterrgridsuball2
a=z5; b=z5ex;  ng=501;  mg=length(a)/ng; ploterrgridsuball
a=z6; b=z6ex;  ng=151;  mg=length(a)/ng; ploterrgridsuball

figure(ifig2), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)
figure(ifig3), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)
figure(ifig4), axis square, axis equal, axis(axisvec)
set(gca,'fontsize',isize)

%ploterrfinishq
