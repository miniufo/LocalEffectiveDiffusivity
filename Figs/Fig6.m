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

[tcoord, ycoord]=meshgrid(tvect,yvect);

keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','nkeff9','t',[1,tgrid]);
tr9=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNGrd.ctl','tr9','t',[1,tgrid]);

keff=reshape(keff,[ygrid,tgrid]);
tr9=reshape(tr9,[ygrid,tgrid]);

keff(keff==-9.99e8)=nan;
tr9(tr9==-9.99e8)=nan;

subplot(2,1,1)
pcolor(tcoord,ycoord,keff);
shading flat;
colorbar;
caxis([0,1400]);
box on;
hold on;
[c, h]=contour(tcoord,ycoord,tr9,[1:0.1:3, 2.34],'LineColor','k','LineWidth',1.2);
clabel(c,h);
set(gca,'xlim',[1 tgrid],'xtick',(0:360:tgrid),'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
scatter(930,470,100,'k','filled');
xlabel('time (day)');
ylabel('Y-coordinate (km)');
title('evolution of tracer and Keff in transformed coord.');
hold off;



[xcoord, ycoord]=meshgrid(xvect,yvect);

tstep=930;

keff=read_grads('I:/cartRL_advSchemes/Leith1_k0/KeffNXY.ctl','nkeff9','t',[tstep,tstep]);
tr9=read_grads('I:/cartRL_advSchemes/Leith1_k0/StatN.ctl','tr9','t',[tstep,tstep]);

keff=reshape(keff,[xgrid,ygrid])';
tr9=reshape(tr9,[xgrid,ygrid])';

keff(keff==-9.99e8)=nan;
tr9(tr9==-9.99e8)=nan;

subplot(2,2,3)
pcolor(xcoord,ycoord,tr9);
shading flat;
colorbar southoutside;
caxis([1,3]);
box on;
hold on
contour(xcoord,ycoord,tr9,[2.34 2.34],'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
xlabel('x-coordinate (km)');
ylabel('y-coordinate (km)');
title('tracer in physical coord. (930 day)');

subplot(2,2,4)
pcolor(xcoord,ycoord,keff);
shading flat;
colorbar southoutside;
caxis([0,1400]);
box on;
hold on
contour(xcoord,ycoord,tr9,[2.34 2.34],'LineColor','k','LineWidth',0.8);
set(gca,'xlim',[0 xgrid-1]*res,'xtick',0:500:(xgrid-1)*res,'ylim',[0 ygrid-1]*res,'ytick',0:500:(ygrid-1)*res);
xlabel('x-coordinate (km)');
ylabel('y-coordinate (km)');
title('Keff in physical coord. (930 day)');

