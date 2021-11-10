clear;clc;

path='D:/Data/MITgcm/barotropicDG/BetaCartRL/cartRL_advSchemes/Leith1_k0/';

tgrid=3600;
xgrid=559;
ygrid=399;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr1v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr1g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum1=df./trg2/2/86400;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr2v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr2g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum2=df./trg2/2/86400;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr4v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr4g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum3=df./trg2/2/86400;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr7v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr7g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum4=df./trg2/2/86400;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr8v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr8g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum5=df./trg2/2/86400;

trv2=read_grads([path,'squaredGrdN.ctl'],'tr9v2','t',[1,tgrid]);
trg2=read_grads([path,'squaredGrdN.ctl'],'tr9g2','t',[1,tgrid]);
trv2=reshape(trv2,[1,tgrid]);
trg2=reshape(trg2,[1,tgrid]);
df=[0 -diff(trv2)];
knum6=df./trg2/2/86400;


smoothKNum=31;
knum1=smooth(knum1,smoothKNum);
knum2=smooth(knum2,smoothKNum);
knum3=smooth(knum3,smoothKNum);
knum4=smooth(knum4,smoothKNum);
knum5=smooth(knum5,smoothKNum);
knum6=smooth(knum6,smoothKNum);
%smoothed(1:60)=nan;
%smoothed(end-60:end)=nan;

subplot(2,3,1)
plot(1:tgrid,knum1,'g',1:tgrid,knum2,'r',1:tgrid,knum3,'b',1:tgrid,knum4,'c',1:tgrid,knum5,'m',1:tgrid,knum6,'k','linewidth',1.2);
xlim([1 tgrid]);
ylim([0 35]);
xlabel('time (model year)');
ylabel('diffusivity (m^2 s^-^1)');
set(gca,'xtick',0:360:tgrid);
set(gca,'xticklabel',0:1:10);
legend('centered 2nd-order','3rd-order upwind','2nd-order DST','3rd-order DST limiter','7th-order one step','2nd-order Prather','location','northeast');
title('spurious diffusivity (Knum)');
