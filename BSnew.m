% load bootsrap data and the data for outlier patients

load('bootstrap_AC_JP_reduced.mat')
load('lowestfive.mat')

warmup = 3;

%Run though 100 BS 
ERRORS_AC_BS = cell(1,length(A_AC_BS));
ERRORS_JP_BS = cell(1,length(A_JP_BS));

parfor i=1:bootstrap_reps
ERRORS_AC_BS{i}= kalman_error(warmup, lowestfive, A_AC_BS{i}, C_AC_BS{i}, Q_AC_BS{i}, R_AC_BS{i}, INITX_AC_BS{i}, INITV_AC_BS{i});
ERRORS_JP_BS{i} = kalman_error(warmup, lowestfive, A_JP_BS{i}, C_JP_BS{i}, Q_JP_BS{i}, R_JP_BS{i}, INITX_JP_BS{i}, INITV_JP_BS{i});
end

save('ErrorTables_bootstrapSample.mat')

%% Obtain the data you need for graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%These are the variables you want to change%%%%%
num_obs = 3;
vis_ahead = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MD_errors_AC_BS = NaN([numpatients, length(vis_ahead),bootstrap_reps]);
MD_errors_JP_BS = NaN([numpatients, length(vis_ahead),bootstrap_reps]);

for v = 1:length(vis_ahead)
    parfor i = 1:bootstrap_reps
        %find the entries in the OHTS table which try to predict v units ahead
        %with num_obs observations as a warm up

        [~, error_array_AC_BS] = extract_errors(ERRORS_AC_BS{i}, num_obs, vis_ahead(v),first_visit);
        [~, error_array_JP_BS] = extract_errors(ERRORS_JP_BS{i}, num_obs, vis_ahead(v),first_visit); 

        %take just the first column of the error_array.
        %column 1 = MD, column 2 = IOP, column 3 = PSD
        MD_array_AC_BS = error_array_AC_BS(:,1);
        MD_array_JP_BS = error_array_JP_BS(:,1);

        %need to know number of entries in this array
        n = length(MD_array_AC_BS);
        m = length(MD_array_JP_BS);

        %record MD errors from JP in this array
        MD_errors_AC_BS(1:n, v,i) = MD_array_AC_BS;
        MD_errors_JP_BS(1:m, v,i) = MD_array_JP_BS;
    end
end


%%
%to obtain a boxplot of the bootstrap estimates for error of patient p from lowest five,
%for v visits ahead, you would make boxplots of
%(1) MD_errors_AC_BS(p,v,:) and
%(2) MD_errors_JP_BS(p,v,:)

