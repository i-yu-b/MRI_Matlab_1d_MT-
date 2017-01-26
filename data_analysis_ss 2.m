clear;close all;clc
folder=['/Users/doesmd/Dropbox/projects/Myelin Goo Swago/data/sample1/'];
format short eng
%% T1
% hard inversion pulse

scan='CLL_T2_IR_02';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
Ni = length(procpar.ti);
Ne = procpar.nume;


% procpar.pw_inv is inversion pulse duration
display(procpar.pw_inv) % ~10 µs

[RE,IM]=load_fid(name);
s=RE+1i*IM;
% 4 points/echo, 1000 echoes, 16 TI, 4step phase cycle
s = reshape(s,4,1000,Ni,4); 
s = squeeze(mean(mean(s,1),4));
    
sig=abs(s).*sign(real(s)); % phase the signal
if sig(1,end)<sig(1,1)
    sig=abs(s).*-sign(real(s));
end

for k = 1:Ni
    sig(:,k) = sig(:,end)-sig(:,k);
end

%% 2D MERA analysis
data.D = sig;
data.t = (1:Ne)'*procpar.te;
data.t2 = procpar.ti;
analysis.graph = 'y';
fitting.twoD = 'y';
output = MERA(data,fitting,analysis)

%% 1D analysis by spin pool
fitting.twoD = 'n';
output = MERA(data,[],analysis)


TI=procpar.ti; % inversion times
N = length(TI);

% monoexponential fit
fx1 = @(params,y,ti) (y/max(y)- ...
   params(1)*(1 - params(2)*exp(-ti*params(3))));
p0 = [1,2,4];
lb = [0 0 0];
ub = [];
[params1,rnorm1,res1,exf1,otp1,~,jc1] =  ...
   lsqnonlin(fx1,p0,lb,ub,[],sig,TI);
P1 = params1(1)*max(sig)*(1-params1(2)*exp(-TI*params1(3)));



% biexponential fit
fx2 = @(params,y,ti) (y/max(y)- ...
   params(1)*(1 - params(2)*exp(-ti*params(3)) ...
   - params(4)*exp(-ti*params(5))));
p0 = [1,1,4,1,1];
lb = [0 0 0 0 0];
ub = [];
[params2,rnorm2,res2,exf2,otp2,~,jc2] =  ...
   lsqnonlin(fx2,p0,lb,ub,[],sig,TI);
P2 = params2(1)*max(sig)*(1-params2(2)*exp(-TI*params2(3))- ...
    params2(4)*exp(-TI*params2(5)));
semilogx(TI,sig,'d',TI,P1,TI,P2)


f12 = (rnorm1-rnorm2)/(5-3)/(rnorm2/(N-5))
fc12 = finv(0.95,2,N-5)

% since f12 > fc12, the higher order model is justified;

figure(1);
semilogx(TI,sig/max(sig),'bd-','MarkerSize',10,'LineWidth',2)
hold on
semilogx(TI,sig/max(sig)-fx2(params2,sig,TI),'r','LineWidth',2)
hold off
axis([min(TI)/2,max(TI)*2,-1.1,1.1])

title('Hard inversion recovery','FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Normalized Signal (a.u.)','FontSize',14)
frac1=params2(2)/(params2(2)+params2(4));
display('Hard inversion pulse:')
display(['Sample is ',num2str(round(frac1*100)), ...
   '% T1 = ',num2str(round(1000*1/params2(3))),' ms']);
display(['Sample is ',num2str(100-round(frac1*100)), ...
   '% T1 = ',num2str(round(1000*1/params2(5))),' ms']);
sig1 = sig;
Tau(1) = procpar.pw_inv;

% from this measurement we can see that the inversion-recovery signal is claerly
% bi-exponential, with component amplitudes around 80/20 for T1s around
% 300/100ms


%% T1
% soft inversion pulse

scan='';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);

% procpar.pw_inv is inversion pulse duration
display(procpar.pw_inv) % ~1 ms

[RE,IM]=load_fid(name);
fid=RE+1i*IM;
%figure();plot(abs(fid));figure();
fid=squeeze(mean(fid(6:10,:)));
fid=reshape(fid(:),length(procpar.c_phase),length(procpar.ti));
fid=squeeze(mean(fid))';
sig=abs(fid).*sign(real(fid)); % phase the signal
if sig(end)<sig(1)
    sig=abs(fid).*-sign(real(fid));
end

TI=procpar.ti; % inversion times
N = length(TI);

% monoexponential fit
fx1 = @(params,y,ti) (y/max(y)- ...
   params(1)*(1 - params(2)*exp(-ti*params(3))));
p0 = [1,2,4];
lb = [0 0 0];
ub = [];
[params1,rnorm1,res1,exf1,otp1,~,jc1] =  ...
   lsqnonlin(fx1,p0,lb,ub,[],sig,TI);



% biexponential fit
fx2 = @(params,y,ti) (y/max(y)- ...
   params(1)*(1 - params(2)*exp(-ti*params(3)) ...
   - params(4)*exp(-ti*params(5))));
p0 = [1,1,4,1,1];
lb = [0 0 0 0 0];
ub = [];
[params2,rnorm2,res2,exf2,otp2,~,jc2] =  ...
   lsqnonlin(fx2,p0,lb,ub,[],sig,TI);

f12 = (rnorm1-rnorm2)/(5-3)/(rnorm2/(N-5));
fc12 = finv(0.95,2,N-5);

% since f12 > fc12, the higher order model is justified;

figure(1);
semilogx(TI,sig/max(sig),'bd-','MarkerSize',10,'LineWidth',2)
hold on
semilogx(TI,sig/max(sig)-fx2(params2,sig,TI),'r','LineWidth',2)
hold off
axis([min(TI)/2,max(TI)*2,-1.1,1.1])

title('Hard inversion recovery','FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Normalized Signal (a.u.)','FontSize',14)
frac1=params2(2)/(params2(2)+params2(4));
display('Soft inversion pulse:')
display(['Sample is ',num2str(round(frac1*100)), ...
   '% T1 = ',num2str(round(1000*1/params2(3))),' ms']);
display(['Sample is ',num2str(100-round(frac1*100)), ...
   '% T1 = ',num2str(round(1000*1/params2(5))),' ms']);

sig2 = sig;
Tau(2) = procpar.pw_inv;


% from this measurement we can see that the inversion-recovery signal is claerly
% bi-exponential, with component amplitudes around 70/30 for T1s around
% 300/100ms. Thus, the shift to lower power inversion results in a larger fast
% component, which is consistent with an MT effect

%% fit hard and soft IR to two pool model w/exchange
% this fit assumes that we're only measuring signal directly from one component
% which is consistent with TE = 1 ms.
clc
% Moa = X(1);
% Mob = X(2);
% kba = X(3);
% R1a = X(4);
% alpha_a = X(5);
% R1b = 1s^-1
% alpha_b =  0.83
s = sig1; %hard ir
s = sig2; %soft ir

X0 = [max(s)/2, max(s)/3,10,2,-0.9]; % initial guess for soft IR, alpha = 0.8
lb = [0 0 0 0 -1];

[Xfit,rnorm,res,exf,out,lam,jac] = ...
    lsqnonlin(@sirfit,X0,[lb],[],[],TI,s);
% [Xfit,R,J,COVB,MSE] = nlinfit(TI,sig,@sirnlin,X0(1:5));
disp(Xfit')
disp(Xfit(2)/Xfit(1))
ci = nlparci(Xfit,res,'jacobian',jac)

% alpha_a fits to ~ -0.65, which corresponds to about 130° flip -- why so low??
% This was probably due to no recalibrating the inversion pulse

% For the soft IR case, using Sb = 0.83, PSR = 1.17 when R1a = R1b, but PSR
% = 0.24 when R1b = 1;

% another demonstration that the 2D analysis is necessary

%% 

% now try fitting both at once

% Moa = X(1);
% Mob = X(2);
% kba = X(3);
% R1a = X(4);
% alpha_a(1) = X(5);
% alpha_a(2) = X(6);
% R2b = X(7);


X0 = [max(sig1)/2, max(sig1)/3,10,2,-0.6,-0.9,1/100e-6]; 
lb = [0 0 0 0 -1 -1 0];

[Xfit,rnorm,res,exf,out,lam,jac] = ...
    lsqnonlin(@sirfit2,X0,[lb],[],[],TI,[sig1;sig2],Tau);

ci = nlparci(Xfit,res,'jacobian',jac)

PSR = Xfit(2)/Xfit(1)

% BUT THERE's A HUGE dependence on R1b
% setting R1b = R1a gives R2b ~ 2microsec!
% setting R1b = 1 gives non-unique solution?!


%% MT Henk
[params CI]=func04132012_07012015(folder);
avg_T1 = 1/Xfit(4);
PSR=params(4)/params(1)/avg_T1;
display(['Henkelman PSR using average hard-inverted T1: ', ...
    num2str(round(PSR*100)),'%']);
legend('w_1/(2*pi) = .127 kHz','w_1/(2*pi) = .256 kHz','w_1/(2*pi) = .383 kHz','w_1/(2*pi) = .639 kHz','w_1/(2*pi) = 1.277 kHz')
macroT2=params(3);
display(['Henkelman exchange rate: ', ...
    num2str(params(1)),' per sec']);
display(['Henkelman macromolecular T2: ', ...
    num2str(params(3)*1e6),' us']);
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
%% MT - Goch/Gore
scan='';
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
figure();semilogy(TI,params(1)-fxn(params,TI),'b-',TI,params(1)-sig,'r*',TI,abs(fxn(params,TI)-sig'),'b--')
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
scan='CLL_T2_IR_01';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
fid=sum(fid,2);
fid=reshape(fid,procpar.ppe,procpar.nume);
fid=sum(fid,1);
TE=[1:procpar.nume]*procpar.te;
MERA_1D(abs(fid),TE,1,'mc');
title('100 us echo time')
%% T2 longerecho
scan='CLL_T2_IR_02';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
fid=sum(fid,2);
fid=reshape(fid,procpar.ppe,procpar.nume);
fid=sum(fid,1);
TE=[1:procpar.nume]*procpar.te;
MERA_1D(abs(fid),TE,1,'mc');
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
figure();semilogy(bvalues*1000,abs(fid),'r*',bvalues*1000,fxn(params,bvalues),'b-')
ADC=params(2);
title('Diffusion, short diffusion time')
xlabel('b (ms/um^2)')
ylabel('Signal (a.u.)')
display(['ADC, diffusion time = ',num2str(procpar.tDELTA*1000),' ms: ',num2str(ADC*1000),' um^2/ms'])

%% ADCSTE
scan='CLL_ADCSTE_03';
name=[folder,scan];
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
% figure();plot(abs(fid));figure();
fid=mean(fid(2:end,:),1);
bvalues=procpar.bvalue;

fxn=inline('params(1)*exp(-bvalues*params(2))','params','bvalues');
options=optimset('Display','off');

bvalues=bvalues(1:end-0);
fid=fid(1:end-0);

params=lsqnonlin(@(params) abs(fid)'-fxn(params,bvalues),[max(abs(fid)),1/max(bvalues)],[0 0],[],options);
figure();semilogy(bvalues*1000,abs(fid),'r*',bvalues*1000,fxn(params,bvalues),'b-')
ADC=params(2);
display(['ADC, diffusion time = ',num2str(procpar.tDELTA*1000),' ms: ',num2str(ADC*1000),' um^2/ms'])
title('Diffusion, long diffusion time')
xlabel('b (ms/um^2)')
ylabel('Signal (a.u.)')



%% IR-CPMG
scan='';
name=[folder,scan];

procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
%figure();plot(abs(fid));figure();
fid=reshape(fid(:),procpar.ppe,procpar.nume,length(procpar.c_phase),length(procpar.ti));
fid = sum(fid,3);
fid=squeeze(mean(fid(6:10,:,:,:)));
te = procpar.te*[1:25];
out = MERA_1D(abs(fid),te,1,'mc');

sig=abs(fid).*sign(real(fid));
if sig(end)<sig(1)
    sig=abs(fid).*-sign(real(fid));
end
TI=procpar.ti;
fxn=inline('params(1) - params(2)*exp(-TI*params(3)) - params(4)*exp(-TI*params(5))','params','TI');
options=optimset('Display','off');
params = lsqnonlin(@(params) sig-fxn(params,TI),[max(sig),(max(sig)-min(sig))/2,1,(max(sig)-min(sig))/2,10],[0 0 0 0 0],[],options);
avg_T1=(params(3)^-1*params(2)+params(5)^-1*params(4))/(params(2)+params(4));
figure(1);semilogy(TI,params(1)-fxn(params,TI),'b-',TI,params(1)-sig,'r*',TI,abs(fxn(params,TI)-sig),'b--')
legend('Fitted M_0 - Fit','Fitted M_0 - Data','Absolute Residuals','Location','NorthEast')
title('Hard inversion recovery')
xlabel('Time (s)')
ylabel('Signal (a.u.)')
frac1=params(2)/(params(2)+params(4));
display('Hard inversion pulse:')
display(['Sample is ',num2str(round(frac1*100)),'% T1 = ',num2str(round(1000*1/params(3))),' ms']);
display(['Sample is ',num2str(100-round(frac1*100)),'% T1 = ',num2str(round(1000*1/params(5))),' ms']);


%% 2D IR-CPMG analysis


Moa = X(1);
Mob = X(2);
kba = X(3);
R1a = X(4);
R2a = X(5);
R2b = X(6);
R1b = 1;
alpha_a = -0.991
alpha_b = 0.83;

scan='';
name=[folder,scan];
display(procpar.pw_inv) % hard inversion
procpar=parsepp([name,'.fid/procpar']);
[RE,IM]=load_fid(name);
fid=RE+1i*IM;
%figure();plot(abs(fid));figure();
fid=reshape(fid(:),procpar.ppe,procpar.nume,length(procpar.c_phase),length(procpar.ti));
fid = sum(fid,3);
fid=squeeze(mean(fid(6:10,:,:,:)));
te = procpar.te*[1:25]';
ti = procpar.ti;


sig3=abs(fid).*sign(real(fid));
if sig3(end)<sig3(1)
    sig3=abs(fid).*-sign(real(fid));
end
% sig3 = zeros(size(fid));
% for k = 1:length(ti)
%    sig3(:,k) = fid(:,end)-fid(:,k);
% end
% sig3 = abs(sig3(:,1:end-1));


% out = MERA_1D(abs(sig3),te,1,'mc');
fop = optimset('Display','iter')

ze = [1:25];
g1 = @(A)(A(1,1)*exp(-te(ze)*A(1,2))+A(2,1)*exp(-te(ze)*A(2,2)))
g = @(A)(abs(sig3(ze,1))-g1(A));
A0 = [max(sig3(:,end)) 40; max(sig3(:,end))/2,12e3];
lb = [0 0; 0 0];
[Afit,rnorm,res,exf,out,lam,jac] = lsqnonlin(g,A0,lb,[],fop);
disp(Afit)



fop = optimset('Display','iter')
X0 = [max(sig3(:)), max(sig3(:))*0.25,10,2,50,2e3]; % initial guess for soft IR, alpha = 0.8
lb = [0 0 0 0 0 1e3];
[Xfit,rnorm,res,exf,out,lam,jac] = lsqnonlin(@ircpmgfit,X0,[lb],[],fop,ti,te(ze),sig3(ze,:));
% [Xfit,R,J,COVB,MSE] = nlinfit(TI,sig,@sirnlin,X0(1:5));
disp(Xfit)
disp(Xfit(2)/Xfit(1))
ci = nlparci(Xfit,res,'jacobian',jac)









