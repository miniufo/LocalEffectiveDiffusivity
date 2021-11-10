clear;clc;

xgrid=559;
ygrid=399;
tgrid=3600;
res=5.5; % km
tvect=(1:1:tgrid);
xvect=(0:1:xgrid-1)*5.5;
yvect=(0:1:ygrid-1)*5.5;
load mycolor;
colormap(myjet);

[xcoord, ycoord]=meshgrid(xvect,yvect);

tstep=1441;

tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr1','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,1)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('centered 2nd-order');


tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr2','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,6)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('3rd-order upwind');


tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr4','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,11)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('2nd-order DST');


tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr7','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,16)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('3rd-order DST limiter');


tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr8','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,21)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('7th-order one step');


tr=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr9','t',[tstep,tstep]);
tr=reshape(tr,[xgrid,ygrid])';
tr(tr==-9.99e8)=nan;

subplot(6,5,26)
pcolor(xcoord,ycoord,tr);
shading flat;
%colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
set(gca,'xticklabel','');
% xlabel('x-coordinate (km)');
% ylabel('y-coordinate (km)');
title('2nd-order Prather');


subplot(6,5,30)
pcolor(xcoord,ycoord,tr);
shading flat;
colorbar southoutside;
caxis([1.4,2.6]);
box on;
hold on
%contour(xcoord,ycoord,tr,1:0.2:3,'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
xlabel('x-coordinate (km)');
ylabel('y-coordinate (km)');
title('2nd-order Prather');




path='D:/Data/MITgcm/barotropicDG/BetaCartRL/cartRL_advSchemes/Leith1_k0/';

tgrid=3600;
xgrid=559;
ygrid=399;

trmax=read_grads([path,'aveVarExN.ctl'],'tr1max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr1ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr1min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,2)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
% set(gca,'xticklabel',0:1:10);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(1700, 4.5, 'centered 2nd order')
title('centered 2nd order');


trmax=read_grads([path,'aveVarExN.ctl'],'tr2max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr2ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr2min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,7)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(1800, 4.5, '3rd order upwind')
title('3rd order upwind');


trmax=read_grads([path,'aveVarExN.ctl'],'tr4max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr4ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr4min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,12)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(2000, 4.5, '2nd order DST')
title('2nd order DST (Lax-Wendroff)');


trmax=read_grads([path,'aveVarExN.ctl'],'tr7max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr7ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr7min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,17)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(1200, 4.5, '3rd order DST with limiter')
title('3rd order DST with limiter');


trmax=read_grads([path,'aveVarExN.ctl'],'tr8max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr8ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr8min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,22)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(1700, 4.5, '7th-order one step')
title('7th-order one step');


trmax=read_grads([path,'aveVarExN.ctl'],'tr9max','t',[1,tgrid]);
trave=read_grads([path,'aveVarExN.ctl'],'tr9ave','t',[1,tgrid]);
trmin=read_grads([path,'aveVarExN.ctl'],'tr9min','t',[1,tgrid]);
trmax=reshape(trmax,[1,tgrid]);
trave=reshape(trave,[1,tgrid]);
trmin=reshape(trmin,[1,tgrid]);

subplot(6,5,27)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
% xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel','');
ylim([-1 5]);
xlim([0 tgrid]);
% legend('max','ave','min','location','southeast');
text(1150, 4.5, '2nd order moment Prather')
title('2nd order moment Prather');

subplot(6,5,29)
plot(1:tgrid,trmax,'red',1:tgrid,trave,'green',1:tgrid,trmin,'blue','linewidth',1.2);
xlabel('time (model year)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel',0:1:10);
ylim([-1 5]);
xlim([0 tgrid]);
legend('max','ave','min','location','southeast');
title('2nd order moment Prather');
