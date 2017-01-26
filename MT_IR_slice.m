function output = MT_IR_slice(folder,scan)
% function fitting T1 data soft pulse Inv-Recovery echo pulse sequence
% input: folder - folder name of current experiment
%        scan - number of the experiment
% output: [PSR kfr], where PSR - pool-size ratio, kfr - exchange rate from
%                    free to restricted pool
%%
scanname = ['ir_1d_slice_select_0',num2str(scan)];
procpar = parsepp([folder,scanname,'.fid/procpar']);
[RE,IM] = load_fid([folder,scanname]);
fid = RE+1i*IM;
TI = procpar.ti;
% average across mid points of echo
fid = squeeze(mean(fid(6:10,:,:)));
data = abs(fid).*sign(real(fid));
if data(end) < data(1)
    data = abs(fid).*-sign(real(fid));
end

%plot data in semilogscale and save
p2 = figure(2);
plot(TI, data,'-o');
title('Inversion-Recovery data')
xlabel('TI, s')
ylabel('M_{free}')

savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
saveas(p2,[folder,'/Results/MT_IR_soft_pulse.fig'])

%% fit data
parameters_fit = lsqnonlin(@(parameters) ...
                          T1_IR_cost(TI,data,parameters,2), ...
                           [max(data),(max(data)-min(data))/2,0.1, ...
                           (max(data)-min(data))/2,1],[0 0 0 0 0],[]);

                       % arrange parameters in right order: 
% Mzf/Mzf_eq=bf_plus*exp(-R1_plus*t)+bf_minus*exp(-R1_minus*t)+1
% R1_plus >> R1_minus
if parameters_fit(5)<parameters_fit(3)
    temp = parameters_fit(2);
    parameters_fit(2) = parameters_fit(4);
    parameters_fit(4) = temp;
            
    temp = parameters_fit(3);
    parameters_fit(3) = parameters_fit(5);
    parameters_fit(5) = temp;
end

Mf_eq = 1; %equlibrium value for free pool
bf_plus = -parameters_fit(2)/parameters_fit(1); %normalize to Meq
bf_minus = -parameters_fit(4)/parameters_fit(1); %normalize to Meq
R1_plus = 1/parameters_fit(3); 
R1_minus = 1/parameters_fit(5); 

fit = parameters_fit(1)*(1+bf_plus*exp(-TI*R1_plus)+ ...
                           bf_minus *exp(-TI*R1_minus));

   
% plot and save data
p3 = figure(3);
plot(TI,data,'ok',TI,fit,'-r',TI,abs(data-fit),'b--');   
title('Inversion-Recovery soft pulse data, biexp fitting')
xlabel('TI, s')
ylabel('log(M_{free}-M_{free_{0}})')
legend('Experimental data','Fitted data','Absolute Residuals', ...
       'Location','NorthEast')
save([folder,'/Results/T1_IR_hard_pulse_fit'])
saveas(p3,[folder,'/Results/T1_IR_hard_pulse_fit.fig'])

% calculate parameters
krf = R1_plus; % approximation
Sm = 0.82; %assumed
PSR = bf_plus/(bf_plus+bf_minus+1-Sm);

output = [PSR krf*PSR];



