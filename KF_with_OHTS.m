%This matlab file helps to look at errors from kalman filter

% Load data

%% Extract error tables for AGIS + CIGTS and OHTS
load('1repSameSize_OHTS_3MD_040117.mat') %have to load the data right before this function because we have repeat name variables

warmup = 3;
ERRORS_AC = kalman_error(warmup, OHTSTest, A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC);
ERRORS_OHTS = kalman_error(warmup, OHTSTest, A_OHTS, C_OHTS, Q_OHTS, R_OHTS, INITX_OHTS, INITV_OHTS);

save('ErrorTables_OHTS.mat')

%% For a fixed number of observations and fixed visits ahead (6 month units), find the errors

load('ErrorTables_OHTS.mat')
[numpatients, ~] = size(ERRORS_OHTS);

%% Obtain the data you need for graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%These are the variables you want to change%%%%%
num_obs = 3;
vis_ahead = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MD_errors_AC = NaN([numpatients, length(vis_ahead)]);
MD_errors_OHTS = NaN([numpatients, length(vis_ahead)]);

for v = 1:length(vis_ahead)

    %find the entries in the OHTS table which try to predict v units ahead
    %with num_obs observations as a warm up
    
    [~, error_array_AC1] = extract_errors(ERRORS_AC, num_obs, vis_ahead(v));
    [~, error_array_OHTS1] = extract_errors(ERRORS_OHTS, num_obs, vis_ahead(v)); 
    
    %take just the first column of the error_array.
    %column 1 = MD, column 2 = IOP, column 3 = PSD
    MD_array_AC1 = error_array_AC1(:,1);
    MD_array_OHTS1 = error_array_OHTS1(:,1);
    
    %need to know number of entries in this array
    n = length(MD_array_AC1);
    m = length(MD_array_OHTS1);
    
    %record MD errors from JP in this array
    MD_errors_AC1(1:n, v) = MD_array_AC1;
    MD_errors_OHTS1(1:m, v) = MD_array_OHTS1;
    
    MD_array_AC1 = rmmissing(MD_array_AC1)
    MD_array_OHTS1 = rmmissing(MD_array_OHTS1)
    
    Y = prctile(MD_array_AC1,25);
    Z = prctile(MD_array_AC1,75);
    I = find(MD_array_AC1 < -.7241 );
    J = find(MD_array_AC1 > 0.4379);
   
    
 
    % Paired t-test 
    %if h=0 then there is no difference between the 2 data sets
    %p is the p-value, significance level (alpha = 0.05)
%     [h,p]=ttest(MD_array_AC1,MD_array_OHTS1);
%     
%     P_value=[h,p];
%     disp(P_value);
%     
%     L= mean(MD_array_AC1);
%     disp(L);
%     
%     M= mean(MD_array_OHTS1);
%     disp(M);
%     
%     O = median (MD_array_AC1);
%     disp(O);
    
%     N = median (MD_array_OHTS1);
%     disp(N);

end
% %% Obtain the data you need for graphs
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%These are the variables you want to change%%%%%
% num_obs = 6;
% vis_ahead = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MD_errors_AC = NaN([numpatients, length(vis_ahead)]);
% MD_errors_OHTS = NaN([numpatients, length(vis_ahead)]);
% 
% for v = 1:length(vis_ahead)
%     %find the entries in the OHTS table which try to predict v units ahead
%     %with num_obs observations as a warm up
%     
%     [~, error_array_AC2] = extract_errors(ERRORS_AC, num_obs, vis_ahead(v));
%     [~, error_array_OHTS2] = extract_errors(ERRORS_OHTS, num_obs, vis_ahead(v)); 
%     
%     %take just the first column of the error_array.
%     %column 1 = MD, column 2 = IOP, column 3 = PSD
%     MD_array_AC2 = error_array_AC2(:,3);
%     MD_array_OHTS2 = error_array_OHTS2(:,3);
%     
%     %need to know number of entries in this array
%     n = length(MD_array_AC2);
%     m = length(MD_array_OHTS2);
%     
%     %record MD errors from JP in this array
%     MD_errors_AC2(1:n, v) = MD_array_AC2;
%     MD_errors_OHTS2(1:m, v) = MD_array_OHTS2;
%    
%     MD_array_AC2 = rmmissing(MD_array_AC2)
%     MD_array_OHTS2 = rmmissing(MD_array_OHTS2)
%     
%     % Paired t-test 
%     %if h=0 then there is no difference between the 2 data sets
%     %p is the p-value, significance level (alpha = 0.05)
%     [h,p]=ttest(MD_array_AC2,MD_array_OHTS2);
%     
%     P_value=[h,p];
%     disp(P_value);
%     
%     L= mean(MD_array_AC2);
%     disp(L);
%     
%     M= mean(MD_array_OHTS2);
%     disp(M);
%     
%     O = median (MD_array_AC2);
%     disp(O);
%     
%     N = median (MD_array_OHTS2);
%     disp(N);
%     
% end

% %% Obtain the data you need for graphs
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%These are the variables you want to change%%%%%
% num_obs = 6;
% vis_ahead = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MD_errors_AC = NaN([numpatients, length(vis_ahead)]);
% MD_errors_OHTS = NaN([numpatients, length(vis_ahead)]);
% 
% for v = 1:length(vis_ahead)
%     %find the entries in the OHTS table which try to predict v units ahead
%     %with num_obs observations as a warm up
%     
%     [~, error_array_AC4] = extract_errors(ERRORS_AC, num_obs, vis_ahead(v));
%     [~, error_array_OHTS4] = extract_errors(ERRORS_OHTS, num_obs, vis_ahead(v)); 
%     
%     %take just the first column of the error_array.
%     %column 1 = MD, column 2 = IOP, column 3 = PSD
%     MD_array_AC4 = error_array_AC4(:,1);
%     MD_array_OHTS4 = error_array_OHTS4(:,1);
%     
%     %need to know number of entries in this array
%     n = length(MD_array_AC4);
%     m = length(MD_array_OHTS4);
%     
%     %record MD errors from JP in this array
%     MD_errors_AC4(1:n, v) = MD_array_AC4;
%     MD_errors_OHTS4(1:m, v) = MD_array_OHTS4;
%    
%     MD_array_AC4 = rmmissing(MD_array_AC4)
%     MD_array_OHTS4 = rmmissing(MD_array_OHTS4)
%     
%     % Paired t-test 
%     %if h=0 then there is no difference between the 2 data sets
%     %p is the p-value, significance level (alpha = 0.05)
%     [h,p]=ttest(MD_array_AC4,MD_array_OHTS4);
%     
%     P_value=[h,p];
%     disp(P_value);
%     
%     L= mean(MD_array_AC4);
%     disp(L);
%     
%     M= mean(MD_array_OHTS4);
%     disp(M);
%     
%     O = median (MD_array_AC4);
%     disp(O);
%     
%     N = median (MD_array_OHTS4);
%     disp(N);
%     
% end
%% 

% Plot boxplots
% box plots of error
figure
hold on

visit_label = {'24 months ahead'}; %strings that we can use to label our plots
    g = [repmat({'MD_errors_AC1'}, 824, 1) ; repmat({'MD_errors_OHTS1'}, 824, 1);repmat({'MD_errors_AC4'}, 824, 1); repmat({'MD_errors_OHTS4'}, 824, 1)];
    boxplot([MD_errors_AC1; MD_errors_OHTS1; MD_errors_AC4; MD_errors_OHTS4],g,'labels',{'AGIS/CIGTS','OHTS','AGIS/CIGTS','OHTS'});
    ylim([-5,5]);
    yl=ylim;
    title(visit_label{v})
  

suptitle(strcat('KF Error for PSD'))

hold off