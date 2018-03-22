%%This matlab file helps to look at errors from kalman filter

% Load data

fprintf('Loading data...\n');
load('Jackknife_AC_OHTS_KFErrors_02-Aug-2017_100reps_newJK_oldkf_excl.mat')

[numpatients_OHTS, ~] = size(OHTS_data);
numpatients_OHTS = numpatients_OHTS - 1; %remove one because of the first row

[numpatients_AC, ~] = size(AC_data);
numpatients_AC = numpatients_AC - 1;


%% Identify patients who have progressed via 2 MD

all_prog_list = zeros(numpatients_OHTS, 1);
all_prog_pid = {};
all_prog_pid_exact = {};
for b = 1:numpatients_OHTS
    pid = pidAll(b);
    all_prog_list(b) = max(OHTS_data{b+1,9}); %did this eye progress?
    if all_prog_list(b) == 1
        all_prog_pid = [all_prog_pid; pid]; %list those patient ids who have progressed
        all_prog_pid_exact = [all_prog_pid_exact; OHTS_data{b+1,1}];
    end
end

OHTS_data_prog = [OHTS_data(1,:), 'FirstProg']; %add a column specifying the first progression
OHTS_prog_loc = {};
row_num = 2;
for b = 1:length(all_prog_pid_exact)
    k = find(strcmp(OHTS_data(:,1), all_prog_pid_exact{b})); %find both eyes belonging to a patient
    for j = 1:length(k)
        prog_time = find(OHTS_data{k(j),9} == 1, 1); %find the first time the progression occurs
        new_row = [OHTS_data(k(j),:), prog_time];
        OHTS_data_prog = [OHTS_data_prog; new_row];
    end
    OHTS_prog_loc = [OHTS_prog_loc; [row_num:row_num+length(k)-1]];
    row_num = row_num + length(k);
end



num_prog = size(OHTS_data_prog,1) - 1;

%% Extract apparent error tables for OHTS patients who have progressed
fprintf('Extracting error tables for AGIS + CIGTS and OHTS...\n');
warmup = 3;
%errors on OHTS data
ERRORS_AC_OHTS = kalman_error_progression_v2(warmup, OHTS_data_prog, A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC, 0);

ERRORS_NULL_OHTS = null_errors_progression_v2(warmup, OHTS_data_prog, 0);

ERRORS_LR_OHTS = regression_errors_progression_v2(warmup, OHTS_data_prog, 0);

ERRORS_LR2_OHTS = regression2_errors_progression_v2(warmup, OHTS_data_prog, 0);


%% set up for error extraction
fprintf('Beginning extraction of errors on AC, Null, LR1, LR2...\n')

%set up variables
num_obs = [3,6];
vis_before = [1,4]; %visits before progression
V = length(vis_before);
N = length(num_obs);

%% extract errors for AC on OHTS and OHTS on AC
error_AC = NaN([N, V, num_prog, 3]);
for n = 1:N
    for v = 1:V
            [~, error_array_AC] = extract_errors_progression(ERRORS_AC_OHTS, num_obs(n), vis_before(v), 0);
        for m = 1:3
            error_AC(n,v,:,m) = squeeze(error_array_AC(:,m));
        end
    end
end

%% extract errors for Null on OHTS
error_NULL_OHTS = NaN([N, V, num_prog, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_NULL_OHTS] = extract_errors_progression(ERRORS_NULL_OHTS, num_obs(n), vis_before(v), 0);
        for m = 1:3
            error_NULL_OHTS(n,v,:,m) = squeeze(error_array_NULL_OHTS(:,m));
        end
    end
end

%% extract linear regression errors on OHTS
error_LR_OHTS = NaN([N, V, num_prog, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR_OHTS] = extract_errors_progression(ERRORS_LR_OHTS, num_obs(n), vis_before(v), 0);
        for m = 1:3
            error_LR_OHTS(n,v,:,m) = squeeze(error_array_LR_OHTS(:,m));
        end
    end
end

%% extract linear regression 2 errors on OHTS
error_LR2_OHTS = NaN([N, V, num_prog, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR2_OHTS] = extract_errors_progression(ERRORS_LR2_OHTS, num_obs(n), vis_before(v), 0);
        for m = 1:3
            error_LR2_OHTS(n,v,:,m) = squeeze(error_array_LR2_OHTS(:,m));
        end
    end
end
%% Find jackknife errors

JK_error_OHTS = NaN([N, V, num_prog, 3]); %visits ahead, number of patients, and MD/IOP/PSD

tic

fprintf('Working on OHTS Jackknife Errors...\n');
%obtain jackknife errors for OHTS
parfor b = 1:num_prog
    JK_data = [OHTS_data_prog(1,:); OHTS_data_prog(b+1,:)];
    %find the parameters corresponding to this patient
    pid = OHTS_data_prog{b+1,1};
    pid = strrep(pid, 'R','');
    pid = strrep(pid, 'L','');
    k = find(strcmp(pid, uniq_pid) == true);
    %performance of bootstrap model on bootstrap sample
    JK_ERRORS_OHTS = kalman_error_progression_v2(warmup, JK_data, A_OHTS_JK{k}, C_OHTS_JK{k}, Q_OHTS_JK{k}, R_OHTS_JK{k}, INITX_OHTS_JK{k}, INITV_OHTS_JK{k}, 0); %OHTS 
    for n = 1:N
        for v = 1:V
            %extract bootstrap model performance
            %OHTS
            [~,error_array_OHTS_JK] = extract_errors_progression(JK_ERRORS_OHTS, num_obs(n), vis_before(v), 0);
            for m = 1:3
                %obtain RMSE and optimism
                %OHTS
                JK_error_OHTS(n,v,b,m) = error_array_OHTS_JK(:,m);
            end
        end
    end
    
end
toc


%% Create RMSE Tables

%PRE RMSE
prog_pts = length(all_prog_pid) ;
RMSE_Pre_OHTS = NaN([N, V, prog_pts, 3]);
RMSE_prog_OHTS = NaN([N, V, prog_pts, 3]);


%find indexes related to progressing eye
prog_ids = OHTS_data_prog(2:end,1);
ids = NaN(length(all_prog_pid_exact),1);
for b = 1:length(all_prog_pid_exact)
    ids(b) = find(strcmp(prog_ids, all_prog_pid_exact{b}));
end


for b = 1:prog_pts
    
    for n = 1:N
        for v = 1:V
            for m = 1:3
                RMSE_Pre_OHTS(n,v,b,m) = (nanmean(squeeze(JK_error_OHTS(n,v,OHTS_prog_loc{b}-1,m)).^2)); %Obtain MSE for each CV fold
                RMSE_prog_OHTS(n,v,b,m) = (nanmean(squeeze(JK_error_OHTS(n,v,ids(b),m)).^2));
            end            
        end
    end
    
end



fprintf('Computing RMSE...\n');
mLabel = {'MD', 'IOP', 'PSD'};

RMSE_Table_OHTS = {'# Observations', 'Months Before Progression', 'Measure', '# Patients', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%OHTS Overall
RMSE_Table_OHTS_prog = {'# Observations', 'Months Before Progression', 'Measure', '# Patients', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%OHTS Progressing eye only





%OHTS testing set
RMSE_OHTS_OHTS = NaN([N, V, 3]);
RMSE_AC_OHTS = NaN([N, V, 3]);
RMSE_NULL_OHTS = NaN([N, V, 3]);
RMSE_LR_OHTS = NaN([N, V, 3]);
RMSE_LR2_OHTS = NaN([N, V, 3]);

RMSE_OHTS_OHTS_prog = NaN([N, V, 3]);
RMSE_AC_OHTS_prog = NaN([N, V, 3]);
RMSE_NULL_OHTS_prog = NaN([N, V, 3]);
RMSE_LR_OHTS_prog = NaN([N, V, 3]);
RMSE_LR2_OHTS_prog = NaN([N, V, 3]);
for n = 1:N
    for m = 1:3
        for v = 1:V

            RMSE_OHTS_OHTS(n,v,m) = sqrt(nanmean(squeeze(RMSE_Pre_OHTS(n,v,:,m))));
            RMSE_OHTS_OHTS_prog(n,v,m) = sqrt(nanmean(squeeze(RMSE_prog_OHTS(n,v,:,m))));
            
            RMSE_AC_OHTS(n,v,m) = sqrt(nanmean(error_AC(n,v,:,m).^2));
            RMSE_AC_OHTS_prog(n,v,m) = sqrt(nanmean(error_AC(n,v,ids,m).^2));
            
            RMSE_LR_OHTS(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,:,m).^2));
            RMSE_LR_OHTS_prog(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,ids,m).^2));
            
            RMSE_LR2_OHTS(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,:,m).^2));
            RMSE_LR2_OHTS_prog(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,ids,m).^2));
            
            RMSE_NULL_OHTS(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,:,m).^2));
            RMSE_NULL_OHTS_prog(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,ids,m).^2));
            
            numpats = sum(~isnan(error_AC(n,v,:,m)));
            
            numpats2 = sum(~isnan(error_AC(n,v,ids,m)));
            
            RMSE_Table_OHTS = [RMSE_Table_OHTS; {num_obs(n), vis_before(v)*6, mLabel{m}, numpats, RMSE_OHTS_OHTS(n,v,m), RMSE_AC_OHTS(n,v,m), ...
                RMSE_LR_OHTS(n,v,m), RMSE_LR2_OHTS(n,v,m), RMSE_NULL_OHTS(n,v,m)}];            
            
            RMSE_Table_OHTS_prog = [RMSE_Table_OHTS_prog; {num_obs(n), vis_before(v)*6, mLabel{m}, numpats2, RMSE_OHTS_OHTS_prog(n,v,m), RMSE_AC_OHTS_prog(n,v,m), ...
                RMSE_LR_OHTS_prog(n,v,m), RMSE_LR2_OHTS_prog(n,v,m), RMSE_NULL_OHTS_prog(n,v,m)}];
            
            
        end
    end
end

%% Save results
fprintf('Saving results to mat file and xlsx files\n')
fname = ['JK_results_AC_OHTS_LR_LR2_null_num_obs',num2str(num_obs),'_PROGRESSION_',date];
save(fname)
% 
%% Save RMSE Tables as XLSX 
xlsxfname = ['JK_Tables_AC_OHTS_LR_LR2_null_num_obs',num2str(num_obs),'_progressions_only_', date,'_2MD.xlsx'];
xlswrite(xlsxfname,RMSE_Table_OHTS,'OHTS Data')
xlswrite(xlsxfname,RMSE_Table_OHTS_prog,'OHTS Progressors')
%xlswrite(xlsxfname,RMSE_Table_OHTS_non,'OHTS Non Progressors')
%xlswrite(xlsxfname,RMSE_Table_AC,'AGIS CIGTS Data')



fprintf('Finished!\n')
%%
fprintf('Sending e-mail \n')
%%send e-mail alert when code is done
recipient = 'garciagg@umich.edu';
subject = 'finished with bootstrap RMSE code';

finished_email(recipient, subject);
