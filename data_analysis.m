clear;clc;close all;
display(' ')
warning off
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20161227_02/';

folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Matlab/Analysis/Irina/Tested/4.7T/BSA_15/cll_47_20150630_01/'
%% T1
scan='CLL_T1_morepts_01';
name=[folder,scan];

procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
%figure();plot(abs(fid));figure();
fid=squeeze(mean(fid(6:10,:)));
fid=reshape(fid(:),length(procpar.c_phase),length(procpar.ti));
fid=squeeze(mean(fid));
sig=abs(fid).*sign(real(fid));
if sig(end)<sig(1)
    sig=abs(fid).*-sign(real(fid));
end
TI=procpar.ti;
fxn=inline('params(1) - params(2)*exp(-TI*params(3)) - params(4)*exp(-TI*params(5))','params','TI');
options=optimset('Display','off');
params = lsqnonlin(@(params) sig'-fxn(params,TI),[max(sig),(max(sig)-min(sig))/2,1,(max(sig)-min(sig))/2,10],[0 0 0 0 0],[],options);
avg_T1=(params(3)^-1*params(2)+params(5)^-1*params(4))/(params(2)+params(4));
figure(1);semilogy(TI,params(1)-fxn(params,TI),'b-',TI,params(1)-sig,'r*',TI,abs(fxn(params,TI)-sig'),'b--')
legend('Fitted M_0 - Fit','Fitted M_0 - Data','Absolute Residuals','Location','NorthEast')
title('Hard inversion recovery')
xlabel('Time (s)')
ylabel('Signal (a.u.)')
frac1=params(2)/(params(2)+params(4));
display('Hard inversion pulse:')
display(['Sample is ',num2str(round(frac1*100)),'% T1 = ',num2str(round(1000*1/params(3))),' ms']);
display(['Sample is ',num2str(100-round(frac1*100)),'% T1 = ',num2str(round(1000*1/params(5))),' ms']);
%% MT Henk
[params CI]=script04132012(folder);
PSR=params(4)/params(1)/avg_T1;
display(['Henkelman PSR using average hard-inverted T1: ',num2str(round(PSR*100)),'%']);
legend('w_1/(2*pi) = .127 kHz','w_1/(2*pi) = .256 kHz','w_1/(2*pi) = .383 kHz','w_1/(2*pi) = .639 kHz','w_1/(2*pi) = 1.277 kHz')
macroT2=params(3);
display(['Henkelman exchange rate: ',num2str(params(1)),' per sec']);
display(['Henkelman macromolecular T2: ',num2str(params(3)*1e6),' us']);
%% FID
scan='spuls_01';
name=[folder,scan];

procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
step=1/procpar.sw;
[dummy, ind]=max(abs(fid));
figure();semilogy(1e6*(step*ind:step:step*length(fid)),abs(fid(ind:end)),'r-^',1e6*(step*ind:step:1e-4),dummy*exp(macroT2^-1*step*ind)*exp(-macroT2^-1*(step*ind:step:1e-4)),'b-')
xlabel('Time (us)')
ylabel('Signal')
title('FID')
display(['Initial slope of FID corresponds to a T2 of ',num2str(1e6*(-(log(abs(fid(ind+1)))-log(abs(fid(ind))))/step)^-1),' us'])
legend('FID','Macromolecular component predicted by Henkelman fit')
% MERA_1D(abs(fid(ind:end)),step*ind:step:step*length(fid),1,'mc');

%% MT Inversion recovery
scan='CLL_T1_morepts_02';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
fid=squeeze(mean(fid(6:10,:)));
fid=reshape(fid(:),length(procpar.c_phase),length(procpar.ti));
fid=squeeze(mean(fid));
sig=abs(fid).*sign(real(fid));
if sig(end)<sig(1)
    sig=abs(fid).*-sign(real(fid));
end
TI=procpar.ti;
fxn=inline('params(1) - params(2)*exp(-TI*params(3)) - params(4)*exp(-TI*params(5))','params','TI');
options=optimset('Display','off');
params = lsqnonlin(@(params) sig'-fxn(params,TI),[max(sig),(max(sig)-min(sig))/2,1,(max(sig)-min(sig))/2,10],[0 0 0 0 0],[],options);
figure();plot(TI,params(1)-fxn(params,TI),'b-',TI,params(1)-sig,'r*',TI,abs(fxn(params,TI)-sig'),'b--')
legend('Fitted M_0 - Fit','Fitted M_0 - Data','Absolute Residuals','Location','NorthEast')
title('Soft inversion recovery')
xlabel('Time (s)')
ylabel('Signal (a.u.)')

params(2)=params(2)/params(1);
params(4)=params(4)/params(1);
params(1)=1;
if params(5)>params(3)
    temp=params(4);
    params(4)=params(2);
    params(2)=temp;
    
    temp=params(5);
    params(5)=params(3);
    params(3)=temp;
end
bfast=params(2);
bslow=params(4);
Rfast=params(3);
Rslow=params(5);
Sm=.82; %assumed
PSR=bfast/(bfast+bslow+1-Sm);
fracfast=bfast/(bfast+bslow);
fracslow=bslow/(bfast+bslow);
display('Soft inversion pulse:')
display(['Sample is ',num2str(round(fracfast*100)),'% T1 = ',num2str(round(1000*1/Rfast)),' ms']);
display(['Sample is ',num2str(round(fracslow*100)),'% T1 = ',num2str(round(1000*1/Rslow)),' ms']);
display(['Goch/Gore PSR = ',num2str(round(PSR*100)),'%']);
display(['Goch/Gore Exchange rate = ',num2str(Rfast),' per sec'])

%% T2
scan='CLL_T2_01';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
fid=sum(fid,3);
fid=reshape(fid,procpar.ppe,procpar.nume);
fid=sum(fid,1);
TE=[1:procpar.nume]*procpar.te;
% MERA_1D(abs(fid),TE,1,'mc');
data.D = abs(fid); data.t = TE'; fitting = []; analysis.graph = 'y';
[~] = MERA(data,fitting,analysis);
title('100 us echo time')
%% T2 longerecho
scan='CLL_T2_longerecho_01';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
fid=sum(fid,3);
fid=reshape(fid,procpar.ppe,procpar.nume);
fid=sum(fid,1);
TE=[1:procpar.nume]*procpar.te;
% MERA_1D(abs(fid),TE,1,'mc');
data.D = abs(fid); data.t = TE'; fitting = []; analysis.graph = 'y';
[~] = MERA(data,fitting,analysis)
title('1 ms echo time')
%% ADCSE
scan='CLL_ADCSE_01';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
% figure();plot(abs(fid));figure();
fid=mean(fid(1:end,:),1);
bvalues=procpar.bvalue;
fxn=inline('params(1)*exp(-bvalues*params(2))','params','bvalues');

bvalues=bvalues(1:end-0);
fid=fid(1:end-0);

options=optimset('Display','off');
params=lsqnonlin(@(params) abs(fid)'-fxn(params,bvalues),[1.1*max(abs(fid)),1/max(bvalues)],[0 0],[],options);
figure();semilogy(bvalues,abs(fid),'r*',bvalues,fxn(params,bvalues),'b-')
ADC=params(2);
title('Diffusion, short diffusion time')
xlabel('Time (s)')
ylabel('Signal (a.u.)')
display(['ADC, diffusion time = ',num2str(procpar.tDELTA*1000),' ms: ',num2str(ADC*1000),' um^2/ms'])

%% ADCSTE
scan='CLL_ADCSTE_01';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
% figure();plot(abs(fid));figure();
fid=mean(fid(1:end,:),1);
bvalues=procpar.bvalue;

fxn=inline('params(1)*exp(-bvalues*params(2))','params','bvalues');
options=optimset('Display','off');

bvalues=bvalues(1:end-0);
fid=fid(1:end-0);

params=lsqnonlin(@(params) abs(fid)'-fxn(params,bvalues),[max(abs(fid)),1/max(bvalues)],[0 0],[],options);
figure();semilogy(bvalues,abs(fid),'r*',bvalues,fxn(params,bvalues),'b-')
ADC=params(2);
display(['ADC, diffusion time = ',num2str(procpar.tDELTA*1000),' ms: ',num2str(ADC*1000),' um^2/ms'])
title('Diffusion, long diffusion time')
xlabel('Time (s)')
ylabel('Signal (a.u.)')
