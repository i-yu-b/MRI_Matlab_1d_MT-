function B1 = B1fit_new(folder,scan1,scan2,T1,T2)
% function fitting B1 values from MT offset resonance experiments 
% input: folder - folder name of current experiment
%        numscans - number of experiment series with different B1 power
%        firstscan - number of the first experiment
%        T1 - T1 value of solution, measured in separate experiment
%        T2 - T2 value of solution, measured in separate experiment
% output: B1 - fitted  B1 amplitude in angular frequency units, s-1

%load data
    scanname = ['mt_henk_1d_slice_select_0',num2str(scan1)];
    [RE,IM] = load_fid([folder, scanname]);
    fid = RE+1i*IM;
    procpar = parsepp([folder,scanname,'.fid/procpar']);
    coarse(:) = procpar.sat_rf_coarse_DG; % attenuation of MT pulse in db
    size(coarse)
    fine(:) = ones(length(coarse),1)*procpar.sat_rf_fine_DG; 
    % attenuation of MT pulse in db
    size(fine)
    offsets = procpar.cest_off_DG; %Hz
    fid=reshape(fid,32,length(coarse),length(offsets));
    fid=permute(fid,[1 3 2]);
    for i=1:length(coarse)
        data(:,i) = squeeze(abs(fid(1,:,i)))/max(squeeze(abs(fid(1,:,i))));
    end
    
    scanname = ['mt_henk_1d_slice_select_0',num2str(scan2)];
    [RE,IM] = load_fid([folder, scanname]);
    fid = RE+1i*IM;
    procpar = parsepp([folder,scanname,'.fid/procpar']);
    coarse(length(coarse)+1) = procpar.sat_rf_coarse_DG;
    size(coarse)
    fine(length(coarse)) = procpar.sat_rf_fine_DG; 
    size(fine)
    data(:,length(coarse)) = squeeze(abs(fid(1,1,:)))/...
                               max(squeeze(abs(fid(1,1,:))));
                           size(data)
 
attenuation_db = coarse(:)+20*log10(min(fine(:),4095)/4095);
attenuation_mag = db2mag(attenuation_db); %in magnitude values
% data = data/max(data(:));
savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
save([folder,'/Results/MT_Henckelman_raw.fig'])

% plot data
p1=figure(1);
semilogx(offsets/1000, data(:,:),'-o');
title ('Magnetization transfer data. Henkelman approach')
xlabel('\Delta\omega, kHz')
ylabel('M_{free}/M_{free_{0}}')

savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
saveas(p1,[folder,'/Results/MT_Henckelman.fig'])

numscans=length(coarse);
%fit data
B1=zeros(numscans,1);
for i = 1:numscans
    Mz = squeeze(data(:,i))/max(squeeze(data(:,i)));
    [B1_fit, ~, residuals, ~, ~, ~, jacobian] = lsqnonlin(@Mzfree_cost,...
                                            1000,0,[],[],offsets,Mz,T1,T2);
    display(['B1(',num2str(i),')= ', num2str(B1_fit), 's-1']);
    B1(i)=squeeze(B1_fit);
end

%plot fitting results


R1 = 1/T1;
fit=zeros(size(data));
for i = 1:numscans
    
    fit(:,i) = R1*(1+(2*pi*offsets*T2).^2)./ (R1*(1+(2*pi*offsets*T2).^2)+...
                                                              B1(i)^2*T2);
end

p2 = figure(2);
semilogx(offsets/1000,data(:,:),'o',offsets/1000,fit(:,:),'-k')
lgd = legend(cellstr(num2str(round(B1/1000/(2*pi),2), '%0.2f')));
lgd.Location = 'northwest';
lgd.Box = 'off';
lgd.Title.String = 'fitted B1 amplitude, kHz';
title ('Magnetization transfer data. Henkelman approach. B1 fitting')
xlabel('\Delta\omega, kHz')
ylabel('M_{free}/M_{free_{0}}')

savedir=[folder,'/Results'];
if  ~exist(savedir)
    mkdir(savedir)
end
saveas(p2,[folder,'/Results/B1_fitting.fig']);
save([folder,'/Results/B1_fit'],'B1');
