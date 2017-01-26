function output = T1_IR_hard(folder,firstscan,flag)
% function fitting T1 data in hard pulse Inv-Recovery echo pulse sequence
% input: folder - folder name of current experiment
%        scan - number of the experiment
%        flag - if flag=1 monoexp fit, if flag =2 biexp fit
% output: if flag ==1 T1f - T1 of free pool
%         if flag ==2 [Mo bf T1f br T1r], where T1f - T1 of free pool
%                                               T1r - T1 of restricted pool

scanname = ['CLL_T1_morepts_0',num2str(firstscan)];
procpar = parsepp([folder,scanname,'.fid/procpar']);
[RE,IM] = load_fid([folder,scanname]);
fid = RE+1i*IM;
TI = procpar.ti;
% average across mid points of echo
fid = squeeze(mean(fid(6:10,:,:)));

% sum up phase cycling data
fid=reshape(fid,length(TI),length(procpar.c_phase));
fid = mean(fid,2);

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
saveas(p2,[folder,'/Results/T1_IR_hard_pulse.fig'])

% fit data
parameters_fit = lsqnonlin(@(parameters) ...
                           T1_IR_cost(TI,data,parameters,flag), ...
                           [max(data),(max(data)-min(data))/2,1, ...
                           (max(data)-min(data))/2,10],[0 0 0 0 0],[]);
p2 = figure(2);
M0 = parameters_fit(1); % M0 free
if flag == 1 
        b = parameters_fit(2);
        T1f = parameters_fit(3); %T1 observed of free pool
        fit = M0-b*exp(-TI/T1f);
        output = T1f;
end
    
if flag == 2
        if parameters_fit(5)>parameters_fit(3)
            temp = parameters_fit(2);
            parameters_fit(2) = parameters_fit(4);
            parameters_fit(4) = temp;
            
            temp = parameters_fit(3);
            parameters_fit(3) = parameters_fit(5);
            parameters_fit(5) = temp;
        end
                
        bf = parameters_fit(2);
        T1f = parameters_fit(3); 
        br = parameters_fit(4);
        T1r = parameters_fit(5); 
        output = parameters_fit;
        fit = M0-bf*exp(-TI/T1f)-br*exp(-TI/T1r);
end

% plot and save data
p3 = figure(3);
semilogy(TI,abs(M0-data),'ok',TI,abs(M0-fit),'-r',TI,abs(data-fit),'b--');
if flag == 1
    title('Inversion-Recovery hard pulse data, monoexp fitting')
end
if flag == 2
    title('Inversion-Recovery hard pulse data, biexp fitting')
end
xlabel('TI, s')
ylabel('log(M_{free}-M_{free_{0}})')
legend('Experimental data','Fitted data','Absolute Residuals', ...
       'Location','NorthEast')
save([folder,'/Results/T1_IR_hard_pulse_fit'])
saveas(p3,[folder,'/Results/T1_IR_hard_pulse_fit.fig'])