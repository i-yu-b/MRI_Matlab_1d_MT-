function parameters = MT_Henckelman_new(folder,omega1,scan1,num_power,scan2)
% function fitting MT offset resonance experiments with different B1 power
% input: folder - folder name of current experiment
%        omega1 - B1 amplitude in Hz
%        numscans - number of experiment series with different B1 power
%        firstscan - number of the first experiment
% output: parameters - fitting parameters,
%         where parameters(1) - R, Henckelman exchange rate, s-1
%               parameters(2) - R1 rate of restricted pool, s-1
%               parameters(3) - T2 of restricted pool, s
%               parameters(4) - R*Mor/R1f, where Mor - Mo value 
%                                                      of restricted pool, 
%                                                R1f - R1 rate of free pool
%               parameters(5) - 1/(R1f*T2f), where T1f - T1 value of free 
%                                                        pool

    if scan1<10
        scanname = ['mt_henk_1d_slice_select_0',num2str(scan1)];
    else
        scanname = ['mt_henk_1d_slice_select_',num2str(scan1)];
    end
    [RE,IM] = load_fid([folder, scanname]);
    fid = RE+1i*IM;
    procpar = parsepp([folder,scanname,'.fid/procpar']);
    coarse(1:num_power) = procpar.sat_rf_coarse_DG(1:num_power); % attenuation of MT pulse in db
    fine(1:num_power) = ones(num_power,1)*procpar.sat_rf_fine_DG; 
    % attenuation of MT pulse in db
    offsets = procpar.cest_off_DG; %Hz
    fid=reshape(fid,32,length(procpar.sat_rf_coarse_DG),length(offsets));
    fid=permute(fid,[1 3 2]);
    for i=1:num_power
        data(:,i) = squeeze(abs(fid(1,:,i)))/max(squeeze(abs(fid(1,:,i))));
    end
    
    if scan2<10
        scanname = ['mt_henk_1d_slice_select_0',num2str(scan2)];
    else
        scanname = ['mt_henk_1d_slice_select_',num2str(scan2)];
    end
    [RE,IM] = load_fid([folder, scanname]);
    fid = RE+1i*IM;
    procpar = parsepp([folder,scanname,'.fid/procpar']);
    coarse(num_power+1) = procpar.sat_rf_coarse_DG;
    fine(num_power+1) = procpar.sat_rf_fine_DG; 
    data(:,num_power+1) = squeeze(abs(fid(1,1,:)))/...
                               max(squeeze(abs(fid(1,1,:))));
                           size(data)
 
attenuation_db = coarse(:)+20*log10(min(fine(:),4095)/4095);
attenuation_mag = db2mag(attenuation_db); %in magnitude values

savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
save([folder,'/Results/MT_Henck_raw'])

% plot data
p1=figure(1);
semilogx(offsets/1000, data(:,:),'-o');
title ('Magnetization transfer data. Henkelman approach')
xlabel('\Delta\omega, kHz')
ylabel('M_{free}/M_{free_{0}}')

lgd = legend(cellstr(num2str(round(omega1/(2*pi*1000),2)', '%0.2f')));
lgd.Location = 'northwest';
lgd.Box = 'off';
lgd.Title.String = 'RF amp attentuation \omega/2\pi, kHz';
savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
saveas(p1,[folder,'/Results/MT_Henck_offset.fig'])

% fit data
[parameters, ~, residuals, ~, ~, ~, jacobian] = ... 
lsqnonlin(@superLorentziancost,[20,0.8,1e-4,6,100],[0 0 0 0 0],[],[],...
           omega1,offsets,data);
savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
save([folder,'/Results/MT_Henck_fit'])
% p2=figure(2);
% semilogx(offsets/1000, data(:,:),'o',);
% title ('Magnetization transfer data. Henkelman approach. Fit')
% xlabel('\Delta\omega, kHz')
% ylabel('M_{free}/M_{free_{0}}')
% 
% lgd = legend(cellstr(num2str(round(omega1/1000,2)', '%0.2f')));
% lgd.Location = 'northwest';
% lgd.Box = 'off';
% lgd.Title.String = 'RF amp attentuation, kHz';
% savedir=[folder,'/Results'];
% if  ~exist(savedir)
%     mkdir(savedir)
% end
% saveas(p2,[folder,'/Results/CLL_henk_MT_offset_fit.fig'])


