%%This matlab file helps to look at errors from kalman filter

% Load data

fprintf('Loading data...\n');
load('Jackknife_AC_JP_KFErrors_11-Jun-2017.mat')


%% Extract apparent error tables for AGIS + CIGTS and Japan
fprintf('Extracting error tables for AGIS + CIGTS and Japan...\n');
warmup = 3;

%errors on Japan data
ERRORS_AC_JP = kalman_error_OCT(warmup, JP_data, A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC, 0);

ERRORS_JP_AC = kalman_error_OCT(warmup, AC_data, A_JP, C_JP, Q_JP, R_JP, INITX_JP, INITV_JP, 0);

ERRORS_NULL_AC = null_errors(warmup, AC_data, 0);
ERRORS_NULL_JP = null_errors(warmup, JP_data, 0);

ERRORS_LR_AC = regression_errors(warmup, AC_data, 0);
ERRORS_LR_JP = regression_errors(warmup, JP_data, 0);

ERRORS_LR2_AC = regression2_errors(warmup, AC_data, 0);
ERRORS_LR2_JP = regression2_errors(warmup, JP_data, 0);

[numpatients_JP, ~] = size(JP_data);
numpatients_JP = numpatients_JP - 1; %remove one because of the first row

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

%% extract errors for AC on JP and JP on AC
error_JP = NaN([N, V, numpatients_AC, 3]);
error_AC = NaN([N, V, numpatients_JP, 3]);
for n = 1:N
    for v = 1:V
            [~, error_array_JP] = extract_errors_OCT(ERRORS_JP_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
            [~, error_array_AC] = extract_errors_OCT(ERRORS_AC_JP, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_JP(n,v,:,m) = squeeze(error_array_JP(:,m));
            error_AC(n,v,:,m) = squeeze(error_array_AC(:,m));
        end
    end
end

%% extract errors for Null on AC and Null on JP
error_NULL_AC = NaN([N, V, numpatients_AC, 3]);
error_NULL_JP = NaN([N, V, numpatients_JP, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_NULL_AC] = extract_errors_OCT(ERRORS_NULL_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_NULL_JP] = extract_errors_OCT(ERRORS_NULL_JP, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_NULL_AC(n,v,:,m) = squeeze(error_array_NULL_AC(:,m));
            error_NULL_JP(n,v,:,m) = squeeze(error_array_NULL_JP(:,m));
        end
    end
end

%% extract linear regression errors for on AC and on JP
error_LR_AC = NaN([N, V, numpatients_AC, 3]);
error_LR_JP = NaN([N, V, numpatients_JP, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR_AC] = extract_errors_OCT(ERRORS_LR_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_LR_JP] = extract_errors_OCT(ERRORS_LR_JP, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_LR_AC(n,v,:,m) = squeeze(error_array_LR_AC(:,m));
            error_LR_JP(n,v,:,m) = squeeze(error_array_LR_JP(:,m));
        end
    end
end

%% extract linear regression 2 errors for on AC and on JP
error_LR2_AC = NaN([N, V, numpatients_AC, 3]);
error_LR2_JP = NaN([N, V, numpatients_JP, 3]);

for n = 1:N
    for v = 1:V
        [~, error_array_LR2_AC] = extract_errors_OCT(ERRORS_LR2_AC, num_obs(n), vis_ahead(v), first_visit(n), 0);
        [~, error_array_LR2_JP] = extract_errors_OCT(ERRORS_LR2_JP, num_obs(n), vis_ahead(v), first_visit(n), 0);
        for m = 1:3
            error_LR2_AC(n,v,:,m) = squeeze(error_array_LR2_AC(:,m));
            error_LR2_JP(n,v,:,m) = squeeze(error_array_LR2_JP(:,m));
        end
    end
end
%% Find jackknife errors

JK_error_JP = NaN([N, V, numpatients_JP, 3]); %visits ahead, number of patients, and MD/IOP/PSD

JK_error_AC = NaN([N, V, numpatients_AC, 3]); %visits ahead, bootstrap reps, and MD/IOP/PSD

tic

fprintf('Working on Japan Jackknife Errors...\n');
%obtain jackknife errors for Japan
parfor b = 1:numpatients_JP
    JK_data = [JP_data(1,:); JP_data(b+1,:)];
    %performance of bootstrap model on bootstrap sample
    JK_ERRORS_JP = kalman_error_OCT(warmup, JK_data, A_JP_JK{b}, C_JP_JK{b}, Q_JP_JK{b}, R_JP_JK{b}, INITX_JP_JK{b}, INITV_JP_JK{b}, 0); %Japan   
    for n = 1:N
        for v = 1:V
            %extract bootstrap model performance
            %Japan
            [~,error_array_JP_JK] = extract_errors_OCT(JK_ERRORS_JP, num_obs(n), vis_ahead(v), first_visit(n),0);
            for m = 1:3
                %obtain RMSE and optimism
                %Japan
                JK_error_JP(n,v,b,m) = error_array_JP_JK(:,m);
            end
        end
    end
    
end
toc

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

%% BOXPLOTS
% Plot results

%plot errors for JP testing set
fprintf('Creating boxplots for Japan testing set...\n');
measure_label = {'MD', 'IOP', 'PSD'};
mlim(1,:) = [-10,10]; %y axis limits for MD plots
mlim(2,:) = [-10,10]; %y axis limits for IOP plots
mlim(3,:) = [-5,5]; %y axis limits for PSD plots

%comparing Japan with AGIS + CIGTS
for n = 1:N
    for m = 1:3
        g = figure;
        hold on
        for v = 1:V
            p(v) = subplot(2,2,v); %assumes length of visits ahead is 4
            x = [squeeze(JK_error_JP(n,v,:,m)),squeeze(error_AC(n,v,:,m)), squeeze(error_NULL_JP(n,v,:,m)), squeeze(error_LR_JP(n,v,:,m)),squeeze(error_LR2_JP(n,v,:,m))];
            %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
            boxplot(x,'labels', {'JP','AC','Null','LR1','LR2'})
            title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
            ylim(mlim(m,:))
        end
        suptitle({'Japan Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        hold off
        print(g,['Figures/Boxplot_JPTest_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf')

    end
end

%%
% %% Plot differences between Japan and AGIS+CIGTS on Japan
% for n = 1:N
%     for m = 1:3
%         g = figure;
%         hold on
%         for v = 1:V
%             p(v) = subplot(1,V,v);
%             x = [squeeze(JK_error_JP(n,v,:,m))-squeeze(error_AC(n,v,:,m))];
%             %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
%             %distributionPlot(x)
%             boxplot(x,'labels', {'JP-AC'})
%             %xlabel('JP - AC')
%             title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%             ylim(mlim(m,:))
%         end
%         suptitle({'Japan Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         hold off
%         print(g,['Figures/Boxplot_JPTest_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf', '-bestfit')
% 
%     end
% end

% %% violin plot
% for n = 1:N
%     for m = 1:3
%         g = figure;
%         hold on
%         for v = 1:V
%             p(v) = subplot(1,V,v);
%             x = [squeeze(JK_error_JP(n,v,:,m))-squeeze(error_AC(n,v,:,m))];
%             %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
%             distributionPlot(x, 'xNames', {'JP-AC'})
%             %boxplot(x,'labels', {'JP-AC'})
%             title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%             ylim(mlim(m,:))
%         end
%         suptitle({'Japan Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         hold off
%         print(g,['Figures/Violinplot_JPTest_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf', '-bestfit')
% 
%     end
% end
% %% comparing Japan with 3 obs vs. 6 obs on Japan testing set
% for m = 1:3
%     g = figure;
%     hold on
%     for v = 1:V
%         p(v) = subplot(1,V,v);
%         x=[];
%         labels = cell([1 N]);
%         for n = 1:N
%             x = [x, squeeze(JK_error_JP(n,v,:,m))];
%             labels{n} = num2str(num_obs(n));
%         end
%         boxplot(x, 'labels', labels);
%         title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%         ylim(mlim(1,:))
%         xlabel('Number of observations')
%     end
%     suptitle({'Japan Training Set, Japan Testing Set', [measure_label{m}]})
%     set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%     hold off
%     print(g,['Figures/Boxplot_JPTrain_JPTest_comparing_', measure_label{m}, '_numobs_', num2str(num_obs),'_',date,'.pdf'], '-dpdf', '-bestfit')
% end
% 
% %% comparing AC with 3 obs vs. 6 obs on Japan testing set
% for m = 1:3
%     g = figure;
%     hold on
%     for v = 1:V
%         p(v) = subplot(1,V,v);
%         x=[];
%         labels = cell([1 N]);
%         for n = 1:N
%             x = [x, squeeze(error_AC(n,v,:,m))];
%             labels{n} = num2str(num_obs(n));
%         end
%         boxplot(x, 'labels', labels);
%         title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%         ylim(mlim(1,:))
%         xlabel('Number of observations')
%     end
%     suptitle({'AGIS CIGTS Training Set, Japan Testing Set', [measure_label{m}]})
%     set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%     hold off
%     print(g,['Figures/Boxplot_ACTrain_JPTest_comparing_', measure_label{m}, '_numobs_', num2str(num_obs),'_',date,'.pdf'], '-dpdf', '-bestfit')
% end
% %% plot errors for AC testing set
% fprintf('Creating boxplots for AGIS CIGTS testing set...\n');
% 
% for n = 1:N
%     for m = 1:3
%         g = figure;
%         hold on
%         for v = 1:V
%             p(v) = subplot(1,V,v);
%             x = [squeeze(error_JP(n,v,:,m)),squeeze(JK_error_AC(n,v,:,m))];
%             %y = [ones(size(error_JP(n,v,:,m))), 2*ones(size(JK_error_AC(n,v,:,m)))];
%             boxplot(x, 'labels', {'JP', 'AC'})
%             title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%             ylim(mlim(m,:))
%         end
%         suptitle({'AGIS + CIGTS Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         hold off
%         print(g,['Figures/Boxplot_ACTest_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf', '-bestfit')
% 
%     end
% end
% %% Plot differences between Japan and AGIS+CIGTS on AGIS+CIGTS
% for n = 1:N
%     for m = 1:3
%         g = figure;
%         hold on
%         for v = 1:V
%             p(v) = subplot(1,V,v);
%             x = [squeeze(error_JP(n,v,:,m))-squeeze(JK_error_AC(n,v,:,m))];
%             %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
%             %distributionPlot(x)
%             boxplot(x,'labels', {'JP-AC'})
%             %xlabel('JP - AC')
%             title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%             ylim(mlim(m,:))
%         end
%         suptitle({'AGIS CIGTS Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         hold off
%         print(g,['Figures/Boxplot_ACTest_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf', '-bestfit')
% 
%     end
% end
% 
% %% violin plot for AC differences
% for n = 1:N
%     for m = 1:3
%         g = figure;
%         hold on
%         for v = 1:V
%             p(v) = subplot(1,V,v);
%             x = [squeeze(error_JP(n,v,:,m))-squeeze(JK_error_AC(n,v,:,m))];
%             %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
%             distributionPlot(x, 'xNames', {'JP-AC'})
%             %boxplot(x,'labels', {'JP-AC'})
%             title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%             ylim(mlim(m,:))
%         end
%         suptitle({'AGIS CIGTS Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         hold off
%         print(g,['Figures/Violinplot_ACTest_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.pdf'], '-dpdf', '-bestfit')
% 
%     end
% end
% 
% %% comparing AC with 3 obs vs. 6 obs on AC Testing set
% for m = 1:3
%     g = figure;
%     hold on
%     for v = 1:V
%         p(v) = subplot(1,V,v);
%         x=[];
%         labels = cell([1 N]);
%         for n = 1:N
%             x = [x, squeeze(JK_error_AC(n,v,:,m))];
%             labels{n} = num2str(num_obs(n));
%         end
%         boxplot(x, 'labels', labels);
%         title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%         ylim(mlim(1,:))
%         xlabel('Number of observations')
%     end
%     suptitle({'AGIS CIGTS Training Set, AGIS CIGTS Testing Set', [measure_label{m}]})
%     set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%     hold off
%     print(g,['Figures/Boxplot_ACTrain_ACTest_comparing_', measure_label{m}, '_numobs_', num2str(num_obs),'_',date,'.pdf'], '-dpdf', '-bestfit')
% end
% 
% %% comparing JP with 3 obs vs. 6 obs on AC testing set
% for m = 1:3
%     g = figure;
%     hold on
%     for v = 1:V
%         p(v) = subplot(1,V,v);
%         x=[];
%         labels = cell([1 N]);
%         for n = 1:N
%             x = [x, squeeze(error_JP(n,v,:,m))];
%             labels{n} = num2str(num_obs(n));
%         end
%         boxplot(x, 'labels', labels);
%         title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
%         ylim(mlim(1,:))
%         xlabel('Number of observations')
%     end
%     suptitle({'Japan Training Set, AGIS CIGTS Testing Set', [measure_label{m}]})
%     set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%     hold off
%     print(g,['Figures/Boxplot_JPTrain_ACTest_comparing_', measure_label{m}, '_numobs_', num2str(num_obs),'_',date,'.pdf'], '-dpdf', '-bestfit')
% end

%% Create RMSE Tables

fprintf('Computing RMSE...\n');
mLabel = {'MD', 'IOP', 'PSD'};

RMSE_Table_JP = {'# Observations', 'Months Ahead', 'Measure', 'Japan*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%Japan Data

RMSE_Table_AC = {'# Observations', 'Months Ahead', 'Measure', 'Japan', 'AGIS/CIGTS*', 'Linear Regression', 'Linear Regression 2', 'NULL'};%AGIS/CIGTS data

%Japan testing set
RMSE_JP_JP = NaN([N, V, 3]);
RMSE_AC_JP = NaN([N, V, 3]);
RMSE_NULL_JP = NaN([N, V, 3]);
RMSE_LR_JP = NaN([N, V, 3]);
RMSE_LR2_JP = NaN([N, V, 3]);
 %AC Testing set
RMSE_JP_AC = NaN([N, V, 3]);
RMSE_AC_AC = NaN([N, V, 3]);
RMSE_NULL_AC = NaN([N, V, 3]);
RMSE_LR_AC = NaN([N, V, 3]);
RMSE_LR2_AC = NaN([N, V, 3]);

for n = 1:N
    for m = 1:3
        for v = 1:V
            %RMSE_JP_JP(n,v,m) = sqrt(nanmean(JK_error_JP(n,v,:,m).^2));
            RMSE_JP_JP(n,v,m) = sqrt(nanmean(JK_error_JP(n,v,:,m).^2));
            RMSE_AC_AC(n,v,m) = sqrt(nanmean(JK_error_AC(n,v,:,m).^2));

            RMSE_AC_JP(n,v,m) = sqrt(nanmean(error_AC(n,v,:,m).^2));
            RMSE_JP_AC(n,v,m) = sqrt(nanmean(error_JP(n,v,:,m).^2));
            
            RMSE_LR_JP(n,v,m) = sqrt(nanmean(error_LR_JP(n,v,:,m).^2));
            RMSE_LR_AC(n,v,m) = sqrt(nanmean(error_LR_AC(n,v,:,m).^2));
 
            RMSE_LR2_JP(n,v,m) = sqrt(nanmean(error_LR2_JP(n,v,:,m).^2));
            RMSE_LR2_AC(n,v,m) = sqrt(nanmean(error_LR2_AC(n,v,:,m).^2));
            
            RMSE_NULL_JP(n,v,m) = sqrt(nanmean(error_NULL_JP(n,v,:,m).^2));
            RMSE_NULL_AC(n,v,m) = sqrt(nanmean(error_NULL_AC(n,v,:,m).^2));

            RMSE_Table_JP = [RMSE_Table_JP; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_JP_JP(n,v,m), RMSE_AC_JP(n,v,m), ...
                RMSE_LR_JP(n,v,m), RMSE_LR2_JP(n,v,m), RMSE_NULL_JP(n,v,m)}];
            RMSE_Table_AC = [RMSE_Table_AC; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_JP_AC(n,v,m), RMSE_AC_AC(n,v,m), ...
                RMSE_LR_AC(n,v,m), RMSE_LR2_AC(n,v,m), RMSE_NULL_AC(n,v,m)}];
        end
    end
end
%% 
% %% P-values for differences
% 
% fprintf('Computing p-values for different observations...\n')
% PVAL_Table_JP = {'Training set', 'Months Ahead', 'Measure', 'Mean (3 obs)', 'SE (3 obs)', 'Mean (6 obs)', 'SE (6 obs)', 'p-value', '95% Confidence',' Interval'};
% PVAL_Table_AC = {'Training set', 'Months Ahead', 'Measure', 'Mean (3 obs)', 'SE (3 obs)', 'Mean (6 obs)', 'SE (6 obs)', 'p-value', '95% Confidence',' Interval'};
% 
% %Jackknife estimates
% if N >= 2
%     for m = 1:3
%         for v = 1:V
%             %make pvalue table for Japan
%             [~, pval, ci] = ttest(JK_error_JP(1,v,:,m), JK_error_JP(2,v,:,m)); 
%             mean1 = nanmean(JK_error_JP(1,v,:,m));
%             mean2 = nanmean(JK_error_JP(2,v,:,m));
%             ste1 = nanstd(JK_error_JP(1,v,:,m))/sqrt(length(squeeze(JK_error_JP(1,v,:,m))));
%             ste2 = nanstd(JK_error_JP(2,v,:,m))/sqrt(length(squeeze(JK_error_JP(1,v,:,m))));
%             PVAL_Table_JP = [PVAL_Table_JP; {'Japan*', vis_ahead(v)*6, mLabel{m}, mean1, ste2, mean2, ste2, pval, ci(1), ci(2)}];
% 
%             %make pvalue table for AGIS CIGTS
%             [~, pval, ci] = ttest(JK_error_AC(1,v,:,m), JK_error_AC(2,v,:,m)); 
%             mean1 = nanmean(JK_error_AC(1,v,:,m));
%             mean2 = nanmean(JK_error_AC(2,v,:,m));
%             ste1 = nanstd(JK_error_AC(1,v,:,m))/sqrt(length(squeeze(JK_error_AC(1,v,:,m))));
%             ste2 = nanstd(JK_error_AC(2,v,:,m))/sqrt(length(squeeze(JK_error_AC(1,v,:,m))));        
%             PVAL_Table_AC = [PVAL_Table_AC; {'AGIS CIGTS*', vis_ahead(v)*6, mLabel{m}, mean1, ste2, mean2, ste2, pval, ci(1), ci(2)}];
% 
%         end
%     end
% else
%     fprintf('Only have one observation - nothing to compare...\n')
% end
% 
% %Full set estimates
% if N >= 2
%     for m = 1:3
%         for v = 1:V
%             %make pvalue table for Japan
%             [~, pval, ci] = ttest(error_AC(1,v,:,m), error_AC(2,v,:,m)); 
%             mean1 = nanmean(error_AC(1,v,:,m));
%             mean2 = nanmean(error_AC(2,v,:,m));
%             ste1 = nanstd(error_AC(1,v,:,m))/sqrt(length(squeeze(error_AC(1,v,:,m))));
%             ste2 = nanstd(error_AC(2,v,:,m))/sqrt(length(squeeze(error_AC(1,v,:,m))));
%             PVAL_Table_JP = [PVAL_Table_JP; {'AGIS CIGTS', vis_ahead(v)*6, mLabel{m}, mean1, ste2, mean2, ste2, pval, ci(1), ci(2)}];
% 
%             %make pvalue table for AGIS CIGTS
%             [~, pval, ci] = ttest(error_JP(1,v,:,m), error_JP(2,v,:,m)); 
%             mean1 = nanmean(error_JP(1,v,:,m));
%             mean2 = nanmean(error_JP(2,v,:,m));
%             ste1 = nanstd(error_JP(1,v,:,m))/sqrt(length(squeeze(error_JP(1,v,:,m))));
%             ste2 = nanstd(error_JP(2,v,:,m))/sqrt(length(squeeze(error_JP(1,v,:,m))));        
%             PVAL_Table_AC = [PVAL_Table_AC; {'Japan', vis_ahead(v)*6, mLabel{m}, mean1, ste2, mean2, ste2, pval, ci(1), ci(2)}];
% 
%         end
%     end
% else
%     fprintf('Only have one observation - nothing to compare...\n')
% end

%% Save results
fprintf('Saving results to mat file and xlsx files\n')
fname = ['JK_results_AC_JP_LR_LR2_null_num_obs',num2str(num_obs),'_first_vis_',num2str(first_visit),'_',date];
save(fname)
% 
% 
xlsxfname = ['JK_Tables_AC_JP_LR_LR2_null_num_obs',num2str(num_obs),'_first_vis_', num2str(first_visit), '_oldkf_', date,'_v2.xlsx'];
xlswrite(xlsxfname,RMSE_Table_JP,'Japan Data')
xlswrite(xlsxfname,RMSE_Table_AC,'AGIS CIGTS Data')
%xlswrite(xlsxfname,PVAL_Table_JP,'Japan p-values')
%xlswrite(xlsxfname,PVAL_Table_AC,'AGIS CIGTS p-values')


fprintf('Finished!\n')
%%
% fprintf('Sending e-mail \n')
% %%send e-mail alert when code is done
% recipient = 'garciagg@umich.edu';
% subject = 'finished with bootstrap RMSE code';
% 
% finished_email(recipient, subject);
