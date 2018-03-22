%%This matlab file helps to look at errors from kalman filter

% Load data

fprintf('Loading data...\n');
load('Jackknife_AC_OHTS_KFErrors_26-Jan-2018_100reps_newJK_oldkf_physexcl_capint_3MD.mat')


%% Extract apparent error tables for AGIS + CIGTS and OHTS
fprintf('Extracting error tables for AGIS + CIGTS and OHTS...\n');
warmup = 3;

%errors on Japan data
ERRORS_AC_OHTS = kalman_error_OCT(warmup, OHTS_data, A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC, 0);

ERRORS_OHTS_AC = kalman_error_OCT(warmup, AC_data, A_OHTS, C_OHTS, Q_OHTS, R_OHTS, INITX_OHTS, INITV_OHTS, 0);

ERRORS_NULL_AC = null_errors(warmup, AC_data, 0);
ERRORS_NULL_OHTS = null_errors(warmup, OHTS_data, 0);

ERRORS_NULL2_AC = null2_errors(warmup, AC_data, 0);
ERRORS_NULL2_OHTS = null2_errors(warmup, OHTS_data, 0);

ERRORS_LR_AC = regression_errors(warmup, AC_data, 0);
ERRORS_LR_OHTS = regression_errors(warmup, OHTS_data, 0);

ERRORS_LR2_AC = regression2_errors(warmup, AC_data, 0);
ERRORS_LR2_OHTS = regression2_errors(warmup, OHTS_data, 0);

[numpatients_OHTS, ~] = size(OHTS_data);
numpatients_OHTS = numpatients_OHTS - 1; %remove one because of the first row

[numpatients_AC, ~] = size(AC_data);
numpatients_AC = numpatients_AC - 1;




%% obtain optimism for mean error
fprintf('Beginning extraction of errors from Jackknife...\n')

%set up variables
num_obs = [3,6];
vis_ahead = [1,2,3,4];
first_visit = [4,1];
V = length(vis_ahead);
N = length(num_obs);

%% extract errors for AC on OHTS and OHTS on AC
error_JP = NaN([N, V, numpatients_AC, 3]);
error_AC = NaN([N, V, numpatients_OHTS, 3]);
for n = 1:N
    for v = 1:V
            [~, error_array_OHTS] = extract_errors_OCT(ERRORS_OHTS_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
            [~, error_array_AC] = extract_errors_OCT(ERRORS_AC_OHTS, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_JP(n,v,:,m) = squeeze(error_array_OHTS(:,m));
            error_AC(n,v,:,m) = squeeze(error_array_AC(:,m));
        end
    end
end

%% extract errors for Null on AC and Null on JP
error_NULL_AC = NaN([N, V, numpatients_AC, 3]);
error_NULL_OHTS = NaN([N, V, numpatients_OHTS, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_NULL_AC] = extract_errors_OCT(ERRORS_NULL_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_NULL_OHTS] = extract_errors_OCT(ERRORS_NULL_OHTS, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_NULL_AC(n,v,:,m) = squeeze(error_array_NULL_AC(:,m));
            error_NULL_OHTS(n,v,:,m) = squeeze(error_array_NULL_OHTS(:,m));
        end
    end
end

%% extract errors for Null2 on AC and Null on JP
error_NULL2_AC = NaN([N, V, numpatients_AC, 3]);
error_NULL2_OHTS = NaN([N, V, numpatients_OHTS, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_NULL2_AC] = extract_errors_OCT(ERRORS_NULL2_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_NULL2_OHTS] = extract_errors_OCT(ERRORS_NULL2_OHTS, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_NULL2_AC(n,v,:,m) = squeeze(error_array_NULL2_AC(:,m));
            error_NULL2_OHTS(n,v,:,m) = squeeze(error_array_NULL2_OHTS(:,m));
        end
    end
end

%% extract linear regression errors for on AC and on JP
error_LR_AC = NaN([N, V, numpatients_AC, 3]);
error_LR_OHTS = NaN([N, V, numpatients_OHTS, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR_AC] = extract_errors_OCT(ERRORS_LR_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_LR_OHTS] = extract_errors_OCT(ERRORS_LR_OHTS, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_LR_AC(n,v,:,m) = squeeze(error_array_LR_AC(:,m));
            error_LR_OHTS(n,v,:,m) = squeeze(error_array_LR_OHTS(:,m));
        end
    end
end

%% extract linear regression 2 errors for on AC and on JP
error_LR2_AC = NaN([N, V, numpatients_AC, 3]);
error_LR2_OHTS = NaN([N, V, numpatients_OHTS, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR2_AC] = extract_errors_OCT(ERRORS_LR2_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_LR2_OHTS] = extract_errors_OCT(ERRORS_LR2_OHTS, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_LR2_AC(n,v,:,m) = squeeze(error_array_LR2_AC(:,m));
            error_LR2_OHTS(n,v,:,m) = squeeze(error_array_LR2_OHTS(:,m));
        end
    end
end
%% Find jackknife errors

JK_error_OHTS = NaN([N, V, numpatients_OHTS, 3]); %visits ahead, number of patients, and MD/IOP/PSD

JK_error_AC = NaN([N, V, numpatients_AC, 3]); %visits ahead, bootstrap reps, and MD/IOP/PSD

tic

fprintf('Working on OHTS Jackknife Errors...\n');
%obtain jackknife errors for OHTS
parfor b = 1:numpatients_OHTS
    JK_data = [OHTS_data(1,:); OHTS_data(b+1,:)];
    %find the parameters corresponding to this patient
    pid = pidAll{b};
    k = find(strcmp(pid, uniq_pid) == true);
    %performance of bootstrap model on bootstrap sample
    JK_ERRORS_OHTS = kalman_error_OCT(warmup, JK_data, A_OHTS_JK{k}, C_OHTS_JK{k}, Q_OHTS_JK{k}, R_OHTS_JK{k}, INITX_OHTS_JK{k}, INITV_OHTS_JK{k}, 0); %OHTS 
    for n = 1:N
        for v = 1:V
            %extract bootstrap model performance
            %OHTS
            [~,error_array_JP_JK] = extract_errors_OCT(JK_ERRORS_OHTS, num_obs(n), vis_ahead(v), first_visit(n),0);
            for m = 1:3
                %obtain RMSE and optimism
                %OHTS
                JK_error_OHTS(n,v,b,m) = error_array_JP_JK(:,m);
            end
        end
    end
    
end
toc
%%
tic

fprintf('Working on AGIS CIGTS Jackknife Errors...\n');
%obtain jackknife errors for AGIS CIGTS
parfor b = 1:numpatients_AC
    JK_data = [AC_data(1,:); AC_data(b+1,:)];
    %performance of bootstrap model on bootstrap sample
    JK_ERRORS_AC = kalman_error_OCT(warmup, JK_data, A_AC_JK{b}, C_AC_JK{b}, Q_AC_JK{b}, R_AC_JK{b}, INITX_AC_JK{b}, INITV_AC_JK{b}, 0); %AGIS CIGTS  
    for n= 1:N
        for v = 1:V
            %extract bootstrap model performance
            %Japan
            [~,error_array_AC_JK] = extract_errors_OCT(JK_ERRORS_AC, num_obs(n), vis_ahead(v), first_visit(n),0);
            for m = 1:3
                %obtain RMSE and optimism
                %Japan
                JK_error_AC(n,v,b,m) = error_array_AC_JK(:,m);
            end
        end
    end
    
end
toc

%% Identify patients who have progressed via Endpoint Committee

all_prog_list_EC = zeros(numpatients_OHTS, 1);
all_prog_list_MD = zeros(numpatients_OHTS, 1);

for b = 1:numpatients_OHTS
    pid = pidAll(b);
    all_prog_list_EC(b) = max(OHTS_data{b+1,8}); %did this eye progress via EC?
    all_prog_list_MD(b) = max(OHTS_data{b+1,9}); %did this eye progress via MD?
end

prog_list_EC = zeros(nOHTS_pt, 1);
prog_list_MD = zeros(nOHTS_pt, 1);
for b = 1:nOHTS_pt
   k = find(strcmp(pidAll, uniq_pid{b}));
   prog_list_EC(b) = max(all_prog_list_EC(k)); %did at least one eye progress via EC?
   prog_list_MD(b) = max(all_prog_list_MD(k)); %did at least one eye progress via MD?
end



%% Create RMSE Tables

%PRE RMSE
RMSE_Pre_OHTS = NaN([N, V, nOHTS_pt, 3]);
for b = 1:nOHTS_pt
    
    for n = 1:N
        for v = 1:V
            for m = 1:3
                RMSE_Pre_OHTS(n,v,b,m) = (nanmean(squeeze(JK_error_OHTS(n,v,JK_loc_OHTS{b}-1,m)).^2)); %Obtain MSE for each CV fold
            end
            
            
        end
    end
    
end



fprintf('Computing RMSE...\n');
mLabel = {'MD', 'IOP', 'PSD'};

RMSE_Table_OHTS = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL', 'NULL2'};%OHTS Overall
RMSE_Table_OHTS_prog_EC = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL', 'NULL2'};%OHTS Progressors by endpoint committee
RMSE_Table_OHTS_non_EC = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL', 'NULL2'};%OHTS Non-progressors by endpoint committee
RMSE_Table_OHTS_prog_MD = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL','NULL2'};%OHTS Progressors by MD drop
RMSE_Table_OHTS_non_MD = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL','NULL2'};%OHTS Non-progressors by MD drop

RMSE_Table_AC = {'# Observations', 'Months Ahead', 'Measure', 'OHTS', 'AGIS/CIGTS*', 'Linear Regression', 'Linear Regression 2', 'NULL','NULL2'};%AGIS/CIGTS data

%OHTS testing set
RMSE_OHTS_OHTS = NaN([N, V, 3]);
RMSE_AC_OHTS = NaN([N, V, 3]);
RMSE_NULL_OHTS = NaN([N, V, 3]);
RMSE_NULL2_OHTS = NaN([N, V, 3]);
RMSE_LR_OHTS = NaN([N, V, 3]);
RMSE_LR2_OHTS = NaN([N, V, 3]);


RMSE_OHTS_OHTS_prog_EC = NaN([N, V, 3]);
RMSE_AC_OHTS_prog_EC = NaN([N, V, 3]);
RMSE_NULL_OHTS_prog_EC = NaN([N, V, 3]);
RMSE_NULL2_OHTS_prog_EC = NaN([N, V, 3]);
RMSE_LR_OHTS_prog_EC = NaN([N, V, 3]);
RMSE_LR2_OHTS_prog_EC = NaN([N, V, 3]);

RMSE_OHTS_OHTS_non_EC = NaN([N, V, 3]);
RMSE_AC_OHTS_non_EC = NaN([N, V, 3]);
RMSE_NULL_OHTS_non_EC = NaN([N, V, 3]);
RMSE_NULL2_OHTS_non_EC = NaN([N, V, 3]);
RMSE_LR_OHTS_non_EC = NaN([N, V, 3]);
RMSE_LR2_OHTS_non_EC = NaN([N, V, 3]);

RMSE_OHTS_OHTS_prog_MD = NaN([N, V, 3]);
RMSE_AC_OHTS_prog_MD = NaN([N, V, 3]);
RMSE_NULL_OHTS_prog_MD = NaN([N, V, 3]);
RMSE_NULL2_OHTS_prog_MD = NaN([N, V, 3]);
RMSE_LR_OHTS_prog_MD = NaN([N, V, 3]);
RMSE_LR2_OHTS_prog_MD = NaN([N, V, 3]);

RMSE_OHTS_OHTS_non_MD = NaN([N, V, 3]);
RMSE_AC_OHTS_non_MD = NaN([N, V, 3]);
RMSE_NULL_OHTS_non_MD = NaN([N, V, 3]);
RMSE_NULL2_OHTS_non_MD = NaN([N, V, 3]);
RMSE_LR_OHTS_non_MD = NaN([N, V, 3]);
RMSE_LR2_OHTS_non_MD = NaN([N, V, 3]);

 %AC Testing set
RMSE_OHTS_AC = NaN([N, V, 3]);
RMSE_AC_AC = NaN([N, V, 3]);
RMSE_NULL_AC = NaN([N, V, 3]);
RMSE_NULL2_AC = NaN([N, V, 3]);
RMSE_LR_AC = NaN([N, V, 3]);
RMSE_LR2_AC = NaN([N, V, 3]);

%find indices corresponding to progressing and non-progressing eyes
p_EC = find(all_prog_list_EC == 1);
non_EC = find(all_prog_list_EC == 0);
p_MD = find(all_prog_list_MD == 1);
non_MD = find(all_prog_list_MD == 0);

% Obtain RMSE associated with each group, measure, and number of months
% ahead
for n = 1:N
    for m = 1:3
        for v = 1:V
            %RMSE_JP_JP(n,v,m) = sqrt(nanmean(JK_error_JP(n,v,:,m).^2));

            RMSE_OHTS_OHTS(n,v,m) = sqrt(nanmean(squeeze(RMSE_Pre_OHTS(n,v,:,m))));
%             RMSE_OHTS_OHTS_prog_EC(n,v,m) = sqrt(nanmean(squeeze(RMSE_Pre_OHTS(n,v,p_EC,m))));
%             RMSE_OHTS_OHTS_non_EC(n,v,m) = sqrt(nanmean(squeeze(RMSE_Pre_OHTS(n,v,non_EC,m))));
%            RMSE_OHTS_OHTS(n,v,m) = sqrt(nanmean(squeeze(JK_error_OHTS(n,v,:,m)).^2));
            RMSE_OHTS_OHTS_prog_EC(n,v,m) = sqrt(nanmean(squeeze(JK_error_OHTS(n,v,p_EC,m)).^2));
            RMSE_OHTS_OHTS_non_EC(n,v,m) = sqrt(nanmean(squeeze(JK_error_OHTS(n,v,non_EC,m)).^2));
            RMSE_OHTS_OHTS_prog_MD(n,v,m) = sqrt(nanmean(squeeze(JK_error_OHTS(n,v,p_MD,m)).^2));
            RMSE_OHTS_OHTS_non_MD(n,v,m) = sqrt(nanmean(squeeze(JK_error_OHTS(n,v,non_MD,m)).^2));
            RMSE_AC_AC(n,v,m) = sqrt(nanmean(JK_error_AC(n,v,:,m).^2));
            
            RMSE_AC_OHTS(n,v,m) = sqrt(nanmean(error_AC(n,v,:,m).^2));
            RMSE_AC_OHTS_prog_EC(n,v,m) = sqrt(nanmean(error_AC(n,v,p_EC,m).^2));
            RMSE_AC_OHTS_non_EC(n,v,m) = sqrt(nanmean(error_AC(n,v,non_EC,m).^2));
            RMSE_AC_OHTS_prog_MD(n,v,m) = sqrt(nanmean(error_AC(n,v,p_MD,m).^2));
            RMSE_AC_OHTS_non_MD(n,v,m) = sqrt(nanmean(error_AC(n,v,non_MD,m).^2));
            RMSE_OHTS_AC(n,v,m) = sqrt(nanmean(error_JP(n,v,:,m).^2));
            
            RMSE_LR_OHTS(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,:,m).^2));
            RMSE_LR_OHTS_prog_EC(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,p_EC,m).^2));
            RMSE_LR_OHTS_non_EC(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,non_EC,m).^2));
            RMSE_LR_OHTS_prog_MD(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,p_MD,m).^2));
            RMSE_LR_OHTS_non_MD(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,non_MD,m).^2));
            RMSE_LR_AC(n,v,m) = sqrt(nanmean(error_LR_AC(n,v,:,m).^2));
            
            RMSE_LR2_OHTS(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,:,m).^2));
            RMSE_LR2_OHTS_prog_EC(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,p_EC,m).^2));
            RMSE_LR2_OHTS_non_EC(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,non_EC,m).^2));
            RMSE_LR2_OHTS_prog_MD(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,p_MD,m).^2));
            RMSE_LR2_OHTS_non_MD(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,non_MD,m).^2));
            RMSE_LR2_AC(n,v,m) = sqrt(nanmean(error_LR2_AC(n,v,:,m).^2));
            
            RMSE_NULL_OHTS(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,:,m).^2));
            RMSE_NULL_OHTS_prog_EC(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,p_EC,m).^2));
            RMSE_NULL_OHTS_non_EC(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,non_EC,m).^2));
            RMSE_NULL_OHTS_prog_MD(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,p_MD,m).^2));
            RMSE_NULL_OHTS_non_MD(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,non_MD,m).^2));
            RMSE_NULL_AC(n,v,m) = sqrt(nanmean(error_NULL_AC(n,v,:,m).^2));
            
            RMSE_NULL2_OHTS(n,v,m) = sqrt(nanmean(error_NULL2_OHTS(n,v,:,m).^2));
            RMSE_NULL2_OHTS_prog_EC(n,v,m) = sqrt(nanmean(error_NULL2_OHTS(n,v,p_EC,m).^2));
            RMSE_NULL2_OHTS_non_EC(n,v,m) = sqrt(nanmean(error_NULL2_OHTS(n,v,non_EC,m).^2));
            RMSE_NULL2_OHTS_prog_MD(n,v,m) = sqrt(nanmean(error_NULL2_OHTS(n,v,p_MD,m).^2));
            RMSE_NULL2_OHTS_non_MD(n,v,m) = sqrt(nanmean(error_NULL2_OHTS(n,v,non_MD,m).^2));
            RMSE_NULL2_AC(n,v,m) = sqrt(nanmean(error_NULL2_AC(n,v,:,m).^2));

            RMSE_Table_OHTS = [RMSE_Table_OHTS; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS(n,v,m), RMSE_AC_OHTS(n,v,m), ...
                RMSE_LR_OHTS(n,v,m), RMSE_LR2_OHTS(n,v,m), RMSE_NULL_OHTS(n,v,m), RMSE_NULL2_OHTS(n,v,m)}];            
            
            RMSE_Table_OHTS_prog_EC = [RMSE_Table_OHTS_prog_EC; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS_prog_EC(n,v,m), RMSE_AC_OHTS_prog_EC(n,v,m), ...
                RMSE_LR_OHTS_prog_EC(n,v,m), RMSE_LR2_OHTS_prog_EC(n,v,m), RMSE_NULL_OHTS_prog_EC(n,v,m), RMSE_NULL2_OHTS_prog_EC(n,v,m)}];

            RMSE_Table_OHTS_non_EC = [RMSE_Table_OHTS_non_EC; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS_non_EC(n,v,m), RMSE_AC_OHTS_non_EC(n,v,m), ...
                RMSE_LR_OHTS_non_EC(n,v,m), RMSE_LR2_OHTS_non_EC(n,v,m), RMSE_NULL_OHTS_non_EC(n,v,m), RMSE_NULL2_OHTS_non_EC(n,v,m)}];
            
            RMSE_Table_OHTS_prog_MD = [RMSE_Table_OHTS_prog_MD; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS_prog_MD(n,v,m), RMSE_AC_OHTS_prog_MD(n,v,m), ...
                RMSE_LR_OHTS_prog_MD(n,v,m), RMSE_LR2_OHTS_prog_MD(n,v,m), RMSE_NULL_OHTS_prog_MD(n,v,m), RMSE_NULL2_OHTS_prog_MD(n,v,m)}];

            RMSE_Table_OHTS_non_MD = [RMSE_Table_OHTS_non_MD; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS_non_MD(n,v,m), RMSE_AC_OHTS_non_MD(n,v,m), ...
                RMSE_LR_OHTS_non_MD(n,v,m), RMSE_LR2_OHTS_non_MD(n,v,m), RMSE_NULL_OHTS_non_MD(n,v,m), RMSE_NULL2_OHTS_non_MD(n,v,m)}];
            
            RMSE_Table_AC = [RMSE_Table_AC; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_AC(n,v,m), RMSE_AC_AC(n,v,m), ...
                RMSE_LR_AC(n,v,m), RMSE_LR2_AC(n,v,m), RMSE_NULL_AC(n,v,m), RMSE_NULL2_AC(n,v,m)}];
        end
    end
end

%% Save results
fprintf('Saving results to mat file and xlsx files\n')
fname = ['JK_results_AC_OHTS_LR_LR2_null_num_obs',num2str(num_obs),'_first_vis_',num2str(first_visit),'_3MD_',date];
save(fname)
% 
%% Save RMSE Tables as XLSX 
xlsxfname = ['JK_Tables_AC_OHTS_LR_LR2_null2_num_obs',num2str(num_obs),'_first_vis_', num2str(first_visit), '_oldkf_', date,'_v2.xlsx'];
xlswrite(xlsxfname,RMSE_Table_OHTS,'OHTS Data')
xlswrite(xlsxfname,RMSE_Table_OHTS_prog_EC,'OHTS Progressors EC')
xlswrite(xlsxfname,RMSE_Table_OHTS_non_EC,'OHTS Non Progressors EC')
xlswrite(xlsxfname,RMSE_Table_OHTS_prog_MD,'OHTS Progressors MD')
xlswrite(xlsxfname,RMSE_Table_OHTS_non_MD,'OHTS Non Progressors MD')
xlswrite(xlsxfname,RMSE_Table_AC,'AGIS CIGTS Data')



