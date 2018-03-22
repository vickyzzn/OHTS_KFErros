%This matlab file helps to look at errors from kalman filter

% Load data

%% Extract error tables for AGIS + CIGTS and OHTS
 %have to load the data right before this function because we have repeat name variables
load('1repSameSize_AC_Test.mat')

warmup = 3;
ERRORS_AC = kalman_error(warmup, testAC, A_AC, C_AC, Q_AC, R_AC, INITX_AC, INITV_AC);
ERRORS_JP = kalman_error(warmup, testAC, A_JP, C_JP, Q_JP, R_JP, INITX_JP, INITV_JP);

save('ErrorTables_bootstrapSample.mat')

%% For a fixed number of observations and fixed visits ahead (6 month units), find the errors

load('ErrorTables_bootstrapSample.mat')
[numpatients, ~] = size(ERRORS_JP);

%% Obtain the data you need for graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%These are the variables you want to change%%%%%
num_obs = 3;
vis_ahead = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MD_errors_AC = NaN([numpatients, length(vis_ahead)]);
MD_errors_JP = NaN([numpatients, length(vis_ahead)]);

for v = 1:length(vis_ahead)

    %find the entries in the OHTS table which try to predict v units ahead
    %with num_obs observations as a warm up
    
    [~, error_array_AC] = extract_errors(ERRORS_AC, num_obs, vis_ahead(v));
    [~, error_array_JP] = extract_errors(ERRORS_JP, num_obs, vis_ahead(v)); 
    
    %take just the first column of the error_array.
    %column 1 = MD, column 2 = IOP, column 3 = PSD
    MD_array_AC = error_array_AC(:,1);
    MD_array_JP = error_array_JP(:,1);
    
    %need to know number of entries in this array
    n = length(MD_array_AC);
    m = length(MD_array_JP);
    
    %record MD errors from JP in this array
    MD_errors_AC(1:n, v) = MD_array_AC;
    MD_errors_JP(1:m, v) = MD_array_JP;
    
    %Remove cells with NaN data
%     MD_array_AC = rmmissing(MD_array_AC)
%     MD_array_JP = rmmissing(MD_array_JP)

    %Find upper and lower values
    Y = prctile(MD_array_AC,25);
    Z = prctile(MD_array_AC,75);
    
    %Table I shows everything below the 25th percentile
    % and table J shows everything above the 75th
    I = find(MD_array_AC < Y );
    J = find(MD_array_AC > Z);
    
    %Create matrix with five patients from table I (manually put in)
    lower1=testAC(13,:);
    lower2=testAC(25,:);
    lower3=testAC(32,:);
    lower4=testAC(33,:);
    lower5=testAC(37,:);
    header={'ID','VisNum','InterpolatedReadings','baseline','age','race','sex'};

    lowestfive = [header; lower1; lower2; lower3; lower4; lower5];
    save('lowestfive.mat')
    
    %Create matrix with five patients from table J (manually put in)
    highest1=testAC(2,:);
    highest2=testAC(5,:);
    highest3=testAC(7,:);
    highest4=testAC(8,:);
    highest5=testAC(12,:);
    header={'ID','VisNum','InterpolatedReadings','baseline','age','race','sex'};

    highestfive = [header; highest1; highest2; highest3; highest4; highest5];
    save('highestfive.mat')
   
    
 
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
    boxplot([MD_errors_AC1; MD_errors_JP; MD_errors_AC4; MD_errors_OHTS4],g,'labels',{'AGIS/CIGTS','OHTS','AGIS/CIGTS','OHTS'});
    ylim([-5,5]);
    yl=ylim;
    title(visit_label{v})
  

suptitle(strcat('KF Error for PSD'))

hold off