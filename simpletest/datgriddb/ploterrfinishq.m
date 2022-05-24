figure(ifig2) ,colorbar off %, set(gca,'XTickLabel',''), 
text(-1.19,1.155,'(d)','fontsize',isize+2) %,set(gca,'XTickLabel','')
%xlabel('z')
ylabel('r')
%print -depsc testerrgrid1
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,0.]; y=[0.,0.01]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
xlabel('z')
set(gca,'XTick',-1.2:.4:0)

figure(ifig3),colorbar off , set(gca,'YTickLabel','') %, set(gca,'XTickLabel','')
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,0.]; y=[0.,0.01]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
text(-1.19,1.155,'(e)','fontsize',isize+2)
xlabel('z')
set(gca,'XTick',-1.2:.4:0)
%print -depsc testerrgrid2
figure(ifig4), set(gca,'YTickLabel','') %, set(gca,'XTickLabel','')
x=[-1.4,-1.1]; y=[0.8,1.1]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
x=[-0.01,0.]; y=[0.,0.01]; [xx,yy]=meshgrid(x,y); zz=-10*ones(size(xx));
[c,h]=contourf(xx,yy,zz,colvec); set(h,'edgecolor','none'); 
xlabel('z')
set(gca,'XTick',-1.2:.4:0)
text(-1.19,1.155,'(f)','fontsize',isize+2)
text(0.082,-0.002,'<','fontsize',isize-1)
text(0.135,1.105,'+','fontsize',isize+2)
cbh=colorbar('YTick',[-10:1]);
%xlabel('z')
%print -depsc testerrgrid3

%figure(1), print -depsc figerrgrid4 
%figure(2), print -depsc figerrgrid5 
%figure(4), print -depsc figerrgrid6
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
