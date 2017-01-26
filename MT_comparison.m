%%%     analysis code for comparison of MT offset and MT IR
%%%     12/28/16 
% function T1_IR_hard - fitting T1 data in hard pulse Inv-Recovery 
% input: folder - folder name of current experiment
%        scan - number of the experiment
%        flag - if flag=1 monoexp fit, if flag =2 biexp fit
% output: if flag ==1 T1f - T1 of free pool in s
%         if flag ==2 [T1f T1r], where T1f - T1 of free pool in s 
%                                      T1r - T1 of restricted pool in s

% function T2_MET2 = multiexponentional fitting of T2 data 
% input: folder - folder name of current experiment
%        scan - number of the experiment
%        flag - if flag = 1 monoexp fit, if flag =2 multiexponential fit
% output: if flag ==1 T2f - T2 of free pool, s
%         if flag ==2 [T2f T2r], where T2f - T2 of free pool, s 
%                                      T2r - T2 of restricted pool, s 

% function MT_Henckelman - fitting MT offset resonance experiments with 
%                          different B1 power
% input: folder - folder name of current experiment
%        numscans - number of experiment series with different B1 power
%        firstscan - number of the first experiment
% output: parameters - fitting parameters,
%         where parameters(1) - R, Henckelman exchange rate, s-1
%               parameters(2) - R1 rate of restricted pool, s-1
%               parameters(3) - T2 of restricted pool, s
%               parameters(4) - R*Mor/R1f, where Mor - Mo value 
%                                                      of restricted pool, 
%                                          R1f - R1 rate of free pool, s-1
%               parameters(5) - 1/(R1f*T2f), where T1f - T1 value of free 
%                                                        pool, s


% function MT_IR - fitting T1 data soft pulse Inv-Recovery echo pulse seq
% input: folder - folder name of current experiment
%        scan - number of the experiment
% output: [PSR kfr], where PSR - pool-size ratio, kfr - exchange rate from
%                    free to restricted pool
clear all
close all

%% determine T1, T2 MnCl2
% %no slice selection
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170102_02/';
% T1_MnCl2 = T1_IR_hard(folder,1,1);
% T2_MnCl2 = T2_MET2(folder,1,1);
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170112_01/';
T1_MnCl2 = T1_IR_hard_slice(folder,1,1);
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170118_01/';
T2_MnCl2 = T2_MET2(folder,1,1);

display(['T1_MnCl2 = ',num2str(T1_MnCl2*1000),' ms']);
display(['T2_MnCl2 = ',num2str(T2_MnCl2*1000),' ms']);
%% fit B1 values from MT Henckelman data for MnCl2
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170112_01/';
B1 = B1fit_new(folder,9,8,T1_MnCl2,T2_MnCl2); %B1 in angular freq units,s-1
% %no slice selection
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170102_02/';
% B1 = B1fit(folder,6,1,T1_MnCl2,T2_MnCl2); 
%% calculate T1_obs for free pool based on T1 Inversion_Recovery hard pulse 

% % slice selection Agar gel
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170120_01/';
% T1_parameters = T1_IR_hard_slice(folder,1,2);

% %no slice selection
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170105_04/';
% T1_parameters = T1_IR_hard(folder,2,1);

% slice select BSA
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170114_01/';
T1_parameters = T1_IR_hard_slice(folder,2,1);

R1_obs = 1 / T1_parameters; %observed relaxation rate for free pool, s-1
%% fit MT Henckelman offsets for BSA

% %slice select Agar gel
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170120_01/';
% parameters_MT_henck = MT_Henckelman_slice(folder,B1',2,5,1);

% slice select BSA
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170114_01/';
parameters_MT_henck = MT_Henckelman_slice(folder,B1',10,5,9);

% %no slice selection
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170105_04/';
% parameters_MT_henck = MT_Henckelman(folder,B1',6,1);
%% get parameters from MT Henckelman
% calculate R1r from R1obs
R = parameters_MT_henck(1); % Henckelman exchange rate, s-1
R1r = 1; % assume relaxation rate for restricted pool 1 s-1, it has weak 
         % influence on final results, and calculation error thus 1 s-1
         % also can be compared with parameters_MT_henck(2)
R1f = R1_obs / (1 + parameters_MT_henck(4)*(R1r - R1_obs)/(R1r-R1_obs+R));
kfr = parameters_MT_henck(4)*R1f; %exchage rate free-to-restricted pool,s-1
                                  %where  parameters_MT_henck(4) - RMor/R1a
PSR = kfr/R; 
T2r = parameters_MT_henck(3); % T2 of restricted pool, s-1
T2f = 1/(parameters_MT_henck(5)*R1f); % T2 of free pool, s-1


display(['Henkelman PSR using hard-inverted T1_obs: ', ...
        num2str(round(PSR*100)),'%']);
display(['Henkelman T2r = ',num2str(T2r*1e6),' us']);
display(['Henkelman exchange rate restricted-to-free pool = ', ...
          num2str(kfr),' s-1']);


%% MT Inversion Recovery
folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170114_01/';
parameters_MT_IR = MT_IR_slice(folder,3);
% %no slice selection
% folder = '/Users/i_yu_b/Dropbox/Work/VUIIS/Data/47data/s_20170105_04/';
% parameters_MT_IR = MT_IR(folder,1);

%% get parameters from Gochberg IR 
PSR = parameters_MT_IR(1);
kfr = parameters_MT_IR(2);

display('Soft inversion pulse:')
display(['MTIR PSR = ',num2str(round(PSR*100)),'%']);
display(['MTIR exchange rate restricted-to-free pool = ',num2str(kfr),...
        ' s-1'])
%% T2 parameters of BSA
