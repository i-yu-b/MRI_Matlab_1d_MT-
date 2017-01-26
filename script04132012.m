function [params CI]=script04132012(folder)
% folder='c:/users/chris/desktop/for mark/myelin experiments/cll_47_20120517_01/';
gammabar=42.57e6;
figure();

for ii=1:5
    scan=['CLL_henk_MT_0',num2str(ii)];
    [RE,IM]=load_fid([folder, scan]);
    fid=RE+1i*IM;
    procpar=parsepp([folder,scan,'.fid/procpar']);
    coarse(ii)=procpar.sat_rf_coarse_DG;
    fine(ii)=procpar.sat_rf_fine_DG;
    data(:,ii)=squeeze(abs(fid(1,:)));
end
offsets=procpar.cest_off_DG;
relativew1=10.^((coarse+20*log10(min(fine,4095)/4095))/20);
relativew1=relativew1/relativew1(1);
data=data/max(data(:));

startpoint=11;
data=data(startpoint:end,:);
offsets=offsets(startpoint:end);

w1=2*pi*gammabar*(3e-6*relativew1);

%superLorentzian parameter order: R, Rb, T2b, RM0b/Ra, T1a/T2a
options=optimset('Display','off','TolFun',eps,'TolX',eps);

% params=lsqnonlin(@superLorentziancostfxn,[1.5,1,1e-4,.5,40],[0 0 0 0 0],[],options,w1,offsets,data)
params=fminsearch(@superLorentziancostfxn,[1 1 1e-5 1 10],options,w1,offsets,data,1);
[params, ~, residuals, ~, ~, ~, jacobian]=lsqnonlin(@superLorentziancostfxn,params,[0 0 0 0 0],[],options,w1,offsets,data);
xlabel('Saturation offset frequency (Hz)')  
ylabel('Normalized signal intensity (unitless)')
title('Henkelman-style MT measurement')
% params(6)=sqrt(sum(residuals.^2))/length(data(:));
CI=nlparci(params(1:5),residuals, 'jacobian',jacobian);


function g=superlorentzian(w,T2)
step=.001;
theta=0:step:(pi/2);
integrand=sin(theta).*sqrt(2/pi)*T2./abs(3*cos(theta).^2-1).*exp(-2*(w*T2./(3*cos(theta).^2-1)).^2);
g=step*trapz(integrand);


function cost=superLorentziancostfxn(params, w1, offsets, data, varargin)
R=params(1);
Rb=params(2);
Rb=1;
T2b=params(3);
p4=params(4);
p5=params(5);

w1=w1(:)';%row
offsets=offsets(:);%column

Rrfb=zeros([length(offsets), length(w1)]);
for ii=1:length(offsets)
    Rrfb(ii,:)=pi*superlorentzian(2*pi*offsets(ii),T2b) * w1.^2;
end

num=Rb*p4 + Rrfb + Rb + R;
den=p4*(Rb+Rrfb) + (1  +  (1./(2*pi*offsets)).^2*w1.^2*p5).*(Rb+Rrfb+R);
sig=num./den;


semilogx(offsets,sig,'-',offsets,data,'*')
drawnow
cost=sig-data;
cost=cost(:);
if nargin>4
    cost=sum(cost.^2);
end