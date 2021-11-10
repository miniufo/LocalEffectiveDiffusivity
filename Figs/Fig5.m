clear;clc;

xgrid=1;
ygrid=399;
tgrid=1440;
res=5.5; % km
tvect=(1:1:tgrid);
xvect=(0:1:xgrid-1)*5.5;
yvect=(0:1:ygrid-1)*5.5;
smoothKNum=31;
smoothKeff=1;
load mycolor;
colormap(myjet);
fontsize=8;

[tcoord, ycoord]=meshgrid(tvect,yvect);

path='D:/Data/MITgcm/barotropicDG/BetaCartRL/cartRL_advSchemes/Leith1_k0/';

keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff1','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr1','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr1v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr1g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,1)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('centered 2nd order');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff2','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr2','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr2v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr2g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,6)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('3rd upwind');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff4','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr4','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr4v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr4g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,11)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('2nd order DST (Lax-Wendroff)');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff7','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr7','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr7v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr7g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,16)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('3rd order DST with limiter');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff8','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr8','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr8v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr8g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,21)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2); %#ok<*NBRAK>
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('7th order one step');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff9','t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr9','t',[1,tgrid]);
trv2=read_grads([path,'squaredGrdN.ctl'],'tr9v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr9g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum=df./trg2/2/86400;
knum=smooth(knum,smoothKNum);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:tgrid; keff(:,i)=keff(:,i)*knum(i);end
for i=1:ygrid; keff(i,:)=smooth(keff(i,:),smoothKeff);end

subplot(6,5,26)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3,2.34],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
scatter(930,470,50,'k','filled');
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model year)');
% ylabel('Y-coordinate (km)');
% title('2nd order Prather');
disp(num2str(max(max(keff))))
hold off;




disp('   ')





xgrid=1;
ygrid=399;
tgrid=1439;
res=5.5; % km
tvect=(1:1:tgrid);
xvect=(0:1:xgrid-1)*5.5;
yvect=(0:1:ygrid-1)*5.5;
load mycolor;
colormap(myjet);

var='KyPos1';
smoothKeff=17;


[tcoord, ycoord]=meshgrid(tvect,yvect);

path='D:/Data/MITgcm/barotropicDG/BetaCartRL/cartRL_advSchemes/Leith1_k0/';

keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr1.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr1.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,2)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,(1:0.2:3),'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('centered 2nd order');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr2.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr2.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,7)
pcolor(tcoord,ycoord,keff/max(max(keff))*1.3);
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,(1:0.2:3),'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('3rd upwind');
disp(num2str(max(max(keff))/1.3))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr4.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr4.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,12)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,(1:0.2:3),'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('2nd order DST (Lax-Wendroff)');
disp(num2str(max(max(keff))))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr7.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr7.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,17)
pcolor(tcoord,ycoord,keff/max(max(keff))*1.9);
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,(1:0.2:3),'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('3rd order DST with limiter');
disp(num2str(max(max(keff))/1.9))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr8.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr8.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,22)
pcolor(tcoord,ycoord,keff/max(max(keff))*1.6);
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('7th order one step');
disp(num2str(max(max(keff))/1.6))
hold off;



keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr9.ctl',var,'t',[1,tgrid]);
tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/diffAlongContourtr9.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr=reshape(tr,[ygrid,tgrid]);
keff(keff==-9.99e8)=nan;
tr(tr==-9.99e8)=nan;
for i=1:ygrid;
    keff(i,:)=smooth(keff(i,:),smoothKeff);
end

subplot(6,5,27)
pcolor(tcoord,ycoord,keff/max(max(keff))*1.1);
shading flat;
%colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,[1:0.2:3],'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','','fontsize',fontsize);
% xlabel('time (model day)');
% ylabel('Y-coordinate (km)');
% title('2nd order Prather');
disp(num2str(max(max(keff))/1.1))
hold off;



subplot(4,3,12)
pcolor(tcoord,ycoord,keff/max(max(keff)));
shading flat;
colorbar southoutside;
caxis([0,1]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr,(1:0.2:3),'LineColor','k','LineWidth',1.2);
clabel(c,h,'fontsize',fontsize);
set(gca,'xlim',[1 tgrid],'xtick',(0:180:tgrid),'xticklabel',(0:180:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'fontsize',fontsize);
xlabel('time (model day)');
ylabel('Y-coordinate (km)');
title('2nd order Prather');
hold off;


