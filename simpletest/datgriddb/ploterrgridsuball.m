x=a(:,1); y=a(:,2)+4; 
u=a(:,3); v=a(:,4); corrqu=a(:,5); corrqv=a(:,6); corrbu=a(:,7); corrbv=a(:,8); unew=a(:,9); vnew=a(:,10);
uex=b(:,3); vex=b(:,4); corrquex=b(:,5); corrqvex=b(:,6); corrbuex=b(:,7); corrbvex=b(:,8); unewex=b(:,9); vnewex=b(:,10);

uex=uex+corrquex+corrbuex;
vex=vex+corrqvex+corrbvex;

uc=u+corrqu;
vc=v+corrqv;
ucb=u+corrqu+corrbu;
vcb=v+corrqv+corrbv;

ind=find(unew~=0);
length(ind)
if (length(ind)>0) 
   ucb(ind)=unew(ind);
   vcb(ind)=vnew(ind);
end

x=reshape(x,[ng,mg]); y=reshape(y,[ng,mg]); 
uex=reshape(uex,[ng,mg]); vex=reshape(vex,[ng,mg]); 

u=reshape(u,[ng,mg]); v=reshape(v,[ng,mg]); 
uc=reshape(uc,[ng,mg]); vc=reshape(vc,[ng,mg]); 
ucb=reshape(ucb,[ng,mg]); vcb=reshape(vcb,[ng,mg]); 

%FIGURE 2
figure(ifig2)
erru=abs(u-uex);
errv=abs(v-vex);
err=sqrt(erru.^2+errv.^2);
err=max(err,minlev); origerr=max(max(err));
z=log10(err);
[c,h]=contourf( y,x, z,colvec); set(h,'edgecolor','none'); hold on
caxis(cax); colorbar

%FIGURE 3
figure(ifig3)
erru=abs(uc-uex);
errv=abs(vc-vex);
err=sqrt(erru.^2+errv.^2);
err=max(err,minlev); maxerrcorr=max(max(err));
errcorr=err;
z=log10(err);
[c,h]=contourf( y,x, z,colvec); set(h,'edgecolor','none'); hold on
caxis(cax); colorbar

%FIGURE 4
figure(ifig4)
erru=abs(ucb-uex);
errv=abs(vcb-vex);
err=sqrt(erru.^2+errv.^2);
%err=erru;
err=max(err,minlev); 
z=log10(err);
[c,h]=contourf( y,x, z,colvec); set(h,'edgecolor','none'); hold on
caxis(cax); colorbar
