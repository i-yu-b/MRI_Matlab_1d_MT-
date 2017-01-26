function output = T2_mems_MET2(folder,scan,flag)
% function multiexponentional fitting of T2 data 
% input: folder - folder name of current experiment
%        scan - number of the experiment
%        flag - if flag = 1 monoexp fit, if flag =2 multiexponential fit
% output: if flag ==1 T2f - T2 of free pool, s
%         if flag ==2 [T2f T2r], where T2f - T2 of free pool, s 
%                                      T2r - T2 of restricted pool, s 


%% T2 longerecho
    scanname = ['mems_0',num2str(scan)];
    procpar = parsepp([folder,scanname,'.fid/procpar']);
    [RE,IM] = load_fid([folder,scanname]);
    s = RE+1i*IM;
    s = permute(s,[1 3 2]);
    NE = procpar.ne;
    TE = [1:NE]*procpar.te;
    for i = 1 : NE
      s_ap(:,:,i) = apodize2d(s(:,:,i),0.35);
      IMG(:,:,i) = abs(fftshift(ifft2(s(:,:,i),size(s,1),size(s,2)))); 
      IMG2(:,:,i) = abs(fftshift(ifft2(s_ap(:,:,i),size(s,1),size(s,2)))); 
    end 
    imagesc(IMG(:,:,1));
    roi=roipoly();
    for i=1:NE
        data(i) = mean2(IMG(:,:,i).*roi);
    end
    %plot data in semilogscale and save
    p1 = figure(1);
    plot(TE, data,'-o');
    title('T2 relaxation data')
    xlabel('TE, s')
    ylabel('M_{free}')
    savedir=[folder,'/Results'];
    if  ~exist(savedir)
        mkdir(savedir)
    end
    save([folder,'/Results/T2_spinecho_raw'])
    saveas(p1,[folder,'/Results/T2_spinecho.fig'])
%% fit data
if flag ==1 
    parameters_fit = lsqnonlin(@(parameters) ...
                           T2_mono_cost(TE,data,parameters), ...
                           [max(data),0.5],[0 0],[]);
    p2 = figure(2);
    M0 = parameters_fit(1); % M0 free
    T2f = parameters_fit(2); %T2 observed of free pool
    fit = M0*exp(-TE/T2f);
    
    p2 = figure(2);
    semilogy(TE,data,'ok',TE,fit,'-r',TE,abs(data-fit),'b--'); 
    title('T2 relaxation data, monoexp fit')
    xlabel('TE, s')
    ylabel('log(M_{free})')
    legend('Experimental data','Fitted data','Absolute Residuals', ...
           'Location','NorthEast')
    save([folder,'/Results/T2_spinecho_fit'])
    saveas(p2,[folder,'/Results/T2_spinecho_fit.fig'])
    output = T2f;
end
 
if flag ==2 
    cd('/Users/i_yu_b/Dropbox/Work/VUIIS/Matlab/Analysis/Irina/Tested')
    close all;
    p2 = figure(2);
    p3 = figure(3);
    Data.D = data; 
    Data.t = TE'; 
    fitting = []; 
    analysis.graph = 'y';
    out1D = MERA(Data,fitting,analysis)
    saveas(p2,[folder,'/Results/T2_spinecho_fit.fig'])
    saveas(p3,[folder,'/Results/T2_spinecho_fit_t-domain.fig'])
    cd('/Users/i_yu_b/Dropbox/Work/VUIIS/Matlab/Analysis/Irina/Tested/4.7T');

end
