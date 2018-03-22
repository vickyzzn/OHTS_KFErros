%This code is used to obtain JackKnife estimates for AGIS CIGTS vs Japan

%%load data

%clear;
close all;
warning off;


fprintf('Loading data... \n');
% read training and testing data sets
trainACAll=readRawData('rawAGISCIGTS.xlsx');
trainOHTSRaw = readRawDataOHTS('OHTS_exclusion_censored.xlsx');

%lower and upper bounds for MD, IOP, PSD
thresholds = [-30, 5; 0, 60; 0, 16];

%get rid of physiologically impossible readings
[trainOHTSRaw_excl, excl_list_OHTS] = OHTS_exclusion(trainOHTSRaw, thresholds);

%interpolate data (capped)
trainACRaw=interpolateData(trainACAll);
[trainOHTS, excl_list_OHTSint] = interpolateDataOHTS_cap(trainOHTSRaw_excl);

%label data with 2 MD drop
AC_data=labelProgression(trainACRaw)
OHTS_data = labelProgression(trainOHTS)


[nAC, ~] = size(AC_data(2:end,:));
[nOHTS, ~] = size(OHTS_data(2:end,:));

AC_data(:,9) = AC_data(:,8);
AC_data{1,8} = 'Prog';

All_data = [AC_data; OHTS_data(2:end,:)];

[nAll,~] = size(All_data(2:end,:));


% %% Obtain all of the patient IDs for OHTS
% 
% eyeid=OHTS_data(2:end,1);
% pidAll={};
% for i=1:size(eyeid,1)
%     pid=eyeid{i,1};
%     pid=strrep(pid,'AGIS',''); 
%     pid=strrep(pid,'CIGTS','');
%     pid=strrep(pid,'L','');
%     pid=strrep(pid,'R','');
%     pidAll=[pidAll;pid];
% end
% 
% uniq_pid = unique(pidAll); %obtain the unique patient ids in the OHTS data
% [nOHTS_pt, ~] = size(uniq_pid);


%% Obtain all of the patient IDs for OHTS

eyeid=All_data(2:end,1);
pid_All={};
for i=1:size(eyeid,1)
    pid=eyeid{i,1};
    pid=strrep(pid,'L','');
    pid=strrep(pid,'R','');
    pid_All=[pid_All;pid];
end

uniq_all_pid = unique(pid_All); %obtain the unique patient ids in the OHTS data
[nAll_pt, ~] = size(uniq_all_pid);


%% Obtain Kalman Filter using ALL PATIENTS

fprintf('Obtaining apparent performance... \n');

kfreps = 100;

fprintf('Creating Kalman Filter for AGIS CIGTS...\n');
tic
% %AC
% [ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(AC_data,0);
% [A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC, LLAC] = learn_kalman_opt(AC_data(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,100);
% toc
% 
% fprintf('Finished AGIS CIGTS...\n Creating Kalman Filter for OHTS...\n');
% tic
% %OHTS MD IOP PSD
% [ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(OHTS_data,0);
% [A_OHTS, C_OHTS, Q_OHTS, R_OHTS, INITX_OHTS, INITV_OHTS, LLOHTS] = learn_kalman_opt(OHTS_data(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,kfreps);
% toc

%All
tic

[ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(All_data,0);
[A_All, C_All, Q_All, R_All, INITX_All, INITV_All, LLAll] = learn_kalman_opt(All_data(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,kfreps);
toc


%% begin making Leave-one-out Models
fprintf('Beginning LOOCV...\n');

% A_AC_JK = cell([nAC,1]);
% C_AC_JK = cell([nAC,1]);
% Q_AC_JK = cell([nAC,1]);
% R_AC_JK = cell([nAC,1]);
% INITX_AC_JK = cell([nAC,1]);
% INITV_AC_JK = cell([nAC,1]);
% LL = cell([nAC, 1]);
% %JK_data_AC = cell([nAC,1]); not necessary since we know it's just removing
% %the ith person from the dataset
% 
% A_OHTS_JK = cell([nOHTS_pt,1]);
% C_OHTS_JK = cell([nOHTS_pt,1]);
% Q_OHTS_JK = cell([nOHTS_pt,1]);
% R_OHTS_JK = cell([nOHTS_pt,1]);
% INITX_OHTS_JK = cell([nOHTS_pt,1]);
% INITV_OHTS_JK = cell([nOHTS_pt,1]);
%  
A_All_JK = cell([nAll_pt,1]);
C_All_JK = cell([nAll_pt,1]);
Q_All_JK = cell([nAll_pt,1]);
R_All_JK = cell([nAll_pt,1]);
INITX_All_JK = cell([nAll_pt,1]);
INITV_All_JK = cell([nAll_pt,1]);



%% LOOCV for AC
% 
% fprintf('Beginning Jackknife procedure for AGIS CIGTS...\n');
% 
% parfor i=1:nAC
%     
%     fprintf('Working on replication %.0f...\n', i);
%     tic
%     %take bootstrap sample of AC data
%     JK_AC = AC_data;
%     JK_AC(i+1,:) = []; %delete the ith person
%     %obtain bootstrap model of AC
% 	[ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(JK_AC, 0);
% 
% 	[A_AC_JK{i}, C_AC_JK{i}, Q_AC_JK{i}, R_AC_JK{i}, INITX_AC_JK{i}, INITV_AC_JK{i}, LL{i}] = learn_kalman_old(JK_AC(2:end,3), A0, C0, Q0, R0, INITX0, INITV0, 100);
%       toc
% end

%% LOOCV for OHTS

% fprintf('Finished with Jackknife for AGIS CIGTS...\n')
% fprintf('Beginning Jackknife procedure for OHTS...\n')
% 
% JK_loc_OHTS = cell([nOHTS_pt, 1]);
% parfor i=1:nOHTS_pt
%     
%     fprintf('Working on replication %.0f...\n', i);
%     tic
%     %take the first patient
%     pt_id = uniq_pid{i,1};
%     JK_loc_OHTS{i} = find(strcmp(pidAll(:,1), pt_id)) + 1;
%     
%     %take JK sample of OHTS data
%     JK_OHTS = OHTS_data;
%     JK_OHTS(JK_loc_OHTS{i},:) = []; %delete the ith person
%     
%     %obtain JK model of OHTS
% 	[ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(JK_OHTS, 0);
% 	[A_OHTS_JK{i}, C_OHTS_JK{i}, Q_OHTS_JK{i}, R_OHTS_JK{i}, INITX_OHTS_JK{i}, INITV_OHTS_JK{i}, ~] = learn_kalman_old(JK_OHTS(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,kfreps);
% 
%     toc
% end

% LOOCV for All

fprintf('Finished with Jackknife for OHTS...\n')
fprintf('Beginning Jackknife procedure for All...\n')

JK_loc_All = cell([nAll_pt, 1]);
parfor i=1:nAll_pt
    
    fprintf('Working on replication %.0f...\n', i);
    tic
    %take the first patient
    pt_id = uniq_all_pid{i,1};
    JK_loc_All{i} = find(strcmp(pid_All(:,1), pt_id)) + 1;
    
    %take JK sample of OHTS data
    JK_All = All_data;
    JK_All(JK_loc_All{i},:) = []; %delete the ith person
    
    %obtain JK model of OHTS
	[ A0, C0, Q0, R0, INITX0, INITV0 ]=initializeEM_JPOCT(JK_All, 0);
	[A_All_JK{i}, C_All_JK{i}, Q_All_JK{i}, R_All_JK{i}, INITX_All_JK{i}, INITV_All_JK{i}, ~] = learn_kalman_opt(JK_All(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,kfreps);

    toc
end


%% finished

fprintf('Finished with Jackknife. Saving results! \n');
%save output
fname = ['Jackknife_All_KFErrors_',date,'_', num2str(kfreps), 'reps_newJK_oldkf_physexcl_capint.mat'];
save(fname)
%send e-mail alert when code is done
recipient = 'vickyzzn@umich.edu';
subject = 'finished with Jackknife Kalman Filter - OHTS and AGIS/CIGTS';

finished_email(recipient, subject);



