figure(ifig2) ,colorbar off %, set(gca,'XTickLabel',''), 
text(-1.19,1.155,'(a)','fontsize',isize+2) %,set(gca,'XTickLabel','')
%xlabel('z')
ylabel('r')
%print -depsc testerrgrid1
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,-0.0]; y=[0.,.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
%xlabel('z')
set(gca,'XTick',-1.2:.4:0)

figure(ifig3),colorbar off , set(gca,'YTickLabel','') %, set(gca,'XTickLabel','')
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,-0.0]; y=[0.,.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
text(-1.19,1.155,'(b)','fontsize',isize+2)
%xlabel('z')
set(gca,'XTick',-1.2:.4:0)
figure(ifig4), set(gca,'YTickLabel','') %, set(gca,'XTickLabel','')
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,-0.0]; y=[0.,.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
%xlabel('z')
set(gca,'XTick',-1.2:.4:0)
text(-1.19,1.155,'(c)','fontsize',isize+2)
text(0.082,-0.002,'<','fontsize',isize-1)
text(0.155,1.105,'+','fontsize',isize+2)
cbh=colorbar('YTick',[-10:1]);
%xlabel('z')
%print -depsc testerrgrid3

pause
%figure(ifig2), print -depsc Fig06a-sgerr100nocorr 
pause
%figure(ifig3), print -depsc Fig06b-sgerr100corr1
pause
%figure(ifig4), print -depsc Fig06c-sgerr100corr2

%figure(5), print -depsc figerrgrid1 
%figure(6), print -depsc figerrgrid2 
%figure(8), print -depsc figerrgrid3
%figure(2), print -depsc errqn100onecorrnowedge  %strip
%figure(3), print -depsc errqn100twocorrnowedge  %strip+axis
%figure(4), print -depsc errqn100threecorrnowedge  %strip+axis+interp

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
