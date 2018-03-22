%% Plot results



%plot errors for OHTS testing set
fprintf('Creating boxplots for OHTS testing set...\n');
measure_label = {'MD', 'IOP', 'PSD'};
mlim(1,:) = [-12,12]; %y axis limits for MD plots
mlim(2,:) = [-15,15]; %y axis limits for IOP plots
mlim(3,:) = [-10,10]; %y axis limits for PSD plots

%comparing OHTS with AGIS + CIGTS
for n = 1:N
    for m = 1:3
        g = figure;
        hold on
        for v = 1 %run 6 and 24 months ahead
            %p(v) = subplot(2,2,v); %assumes length of visits ahead is 4
            x = [squeeze(JK_error_JP(n,v,:,m)), squeeze(error_AC(n,v,:,m)), squeeze(error_NULL_JP(n,v,:,m)), squeeze(error_LR_JP(n,v,:,m)),squeeze(error_LR2_JP(n,v,:,m))];
            %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
            boxplot(x,'labels', {'JP(w/ JK)','AC','Null','LR1','LR2'})
            title(strcat(num2str(6*vis_ahead(v)), ' months ahead'))
            ylim(mlim(m,:))
            yticks([-10 -5 -2.5 -1 -0.5 0 0.5 1 2.5 5 10])
            set(gca,'fontsize',8)
            line1=refline(0,0);
                line1.Color='k';
            line2=refline(0,.5)
                line2.Color='r'
                line2.LineStyle='-'
            line3=refline(0,-.5)
                line3.Color='r'
                line3.LineStyle='-'
            line4=refline(0,1)
                line4.Color='k'
                line4.LineStyle=':'
            line5=refline(0,-1)
                line5.Color='k'
                line5.LineStyle=':'
            line6=refline(0,2.5)
                line6.Color='b'
                line6.LineStyle='-.'
            line7=refline(0,-2.5)
                line7.Color='b'
                line7.LineStyle='-.'
        end
        suptitle({'JP Testing Set',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
        %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        hold off
        print(g,['Figures/boxplot_jp_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_firstvis_', num2str(first_visit(n)),'_',date,'.png'], '-dpng')
    end
end
%% %% violin plot
fprintf('Creating violin for OHTS testing set...\n');
measure_label = {'MD', 'IOP', 'PSD'};
axis_label= {'OHTS(w/ JK)','AC','Null','LR1','LR2'};
mlim(1,:) = [-12,12]; %y axis limits for MD plots
mlim(2,:) = [-15,15]; %y axis limits for IOP plots
mlim(3,:) = [-10,10]; %y axis limits for PSD plots
vis_before = [1,4];
V = length(vis_before);
for n = 1:N
    for m = 1:3
        g = figure;
        hold on
        for v = 1%run 6 and 24 months ahead
            %p(v) = subplot(1,V,v);
            x = [squeeze(JK_error_OHTS(n,v,p_MD,m)), squeeze(error_AC(n,v,p_MD,m)), squeeze(error_NULL_OHTS(n,v,p_MD,m)), squeeze(error_NULL2_OHTS(n,v,p_MD,m)),squeeze(error_LR_OHTS(n,v,p_MD,m)),squeeze(error_LR2_OHTS(n,v,p_MD,m))];
            %y = [ones(size(squeeze(JK_error_JP(n,v,:,m)))), 2*ones(size(squeeze(error_AC(n,v,:,m))))];
            distributionPlot(x, 'globalNorm', 2)
            %xlabels=(axis_label(v))
            %boxplot(x,'labels', {'JP-AC'})
            title(strcat(num2str(6*vis_before(v)), ' months ahead'))
            ylim(mlim(m,:))
            yticks([-10 -5 -2.5 -1 -0.5 0 0.5 1 2.5 5 10])
            set(gca,'fontsize',8)
            line1=refline(0,0);
                line1.Color='k';
            line2=refline(0,.5)
                line2.Color='r'
                line2.LineStyle='-'
            line3=refline(0,-.5)
                line3.Color='r'
                line3.LineStyle='-'
            line4=refline(0,1)
                line4.Color='k'
                line4.LineStyle=':'
            line5=refline(0,-1)
                line5.Color='k'
                line5.LineStyle=':'
            line6=refline(0,2.5)
                line6.Color='b'
                line6.LineStyle='-.'
            line7=refline(0,-2.5)
                line7.Color='b'
                line7.LineStyle='-.'
                end
        suptitle({'OHTS Prog (MD)',[measure_label{m},' with ',num2str(num_obs(n)), ' observations']})
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        hold off
        print(g,['Figures/Violinplot_OHTS_differences_JK_',measure_label{m},'_numobs_',num2str(num_obs(n)), '_',date,'.png'], '-dpng')

    end
end
%% 


%Ranges for violin plots
labeling= {'Type','# Observations', 'Months Ahead','OHTS','AC','NULL','NULL2','LR1','LR2'}
mLabel = {'MD', 'IOP', 'PSD'};
vis_before = [1,4];
V = length(vis_before);
    for n=1:2
        for m=1:3
            for v=1:2%Change ranges to: +-.5, +-1, +-2.5
                OHTS=JK_error_OHTS(n,v,non_EC,m)
                numberOHTS=OHTS>=-1& OHTS<=1%becomes logical 1 if true or 0 if false
                countOHTS(n,v,m)=sum(numberOHTS) %count number of trues (sum won't count any zeros)
                percentageOHTS(n,v,m)=countOHTS(n,v,m)/numel(numberOHTS); %divide number of trues by the total number trues and falses

                %             All=JK_error_All_OHTS(n,v,:,m)
    %             numberAll=All>=-.5 & All<=.5 %becomes logical 1 if true or 0 if false
    %             countAll(n,v,m)=sum(numberAll) %count number of trues (sum won't count any zeros)
    %             percentageAll(n,v,m)=countAll(n,v,m)/numel(numberAll); %di

                AC=error_AC(n,v,non_EC,m)
                numberAC=AC>=-1& AC<=1
                countAC(n,v,m)=sum(numberAC)
                percentageAC(n,v,m)=countAC(n,v,m)/numel(numberAC);

                NULL=error_NULL_OHTS(n,v,non_EC,m)
                numberNULL=NULL>=-1& NULL<=1
                countNULL(n,v,m)=sum(numberNULL)
                percentageNULL(n,v,m)=countNULL(n,v,m)/numel(numberNULL);


                NULL2=error_NULL2_OHTS(n,v,non_EC,m)
                numberNULL2=NULL2>=-1& NULL2<=1
                countNULL2(n,v,m)=sum(numberNULL2)
                percentageNULL2(n,v,m)=countNULL2(n,v,m)/numel(numberNULL2);


                LR1=error_LR_OHTS(n,v,non_EC,m)
                numberLR1=LR1>=-1& LR1<1
                countLR1(n,v,m)=sum(numberLR1)
                percentageLR1(n,v,m)=countLR1(n,v,m)/numel(numberLR1);

                LR2=error_LR2_OHTS(n,v,non_EC,m)
                numberLR2=LR2>=-1& LR2<=1
                countLR2(n,v,m)=sum(numberLR2)
                percentageLR2(n,v,m)=countLR2(n,v,m)/numel(numberLR2);


                labeling=[labeling;{mLabel{m},num_obs(n), vis_before(v)*6,percentageOHTS(n,v,m),percentageAC(n,v,m),percentageNULL(n,v,m),percentageNULL2(n,v,m),percentageLR1(n,v,m),percentageLR2(n,v,m)}];
                disp(labeling) 
            end
        end
    end
%% repeat outliers
table=NaN([263,5])

labels={'JP', 'AC', 'Null', 'LR1', 'LR2'}
for a= 1:262
    jp_element(a)=numberJP(:,:,a)
    ac_element(a)=numberAC(:,:,a)
    null_element(a)=numberNULL(:,:,a)
    lr1_element(a)=numberLR1(:,:,a)
    lr2_element(a)=numberLR2(:,:,a)
    labels=[labels; {jp_element(a), ac_element(a), null_element(a), lr1_element(a), lr2_element(a)}]
table=[labels;{numberJP(:,:,a), numberAC(:,:,a),numberNULL(:,:,a), numberLR1(:,:,a), numberLR2(:,:,a)}]
table=horzcat(numberJP(:,:,a), numberAC(:,:,a),numberNULL(:,:,a), numberLR1(:,:,a), numberLR2(:,:,a))
end



%% Create RMSE Tables

% %PRE RMSE
% RMSE_Pre_OHTS = NaN([N, V, nOHTS_pt, 3]);
% for b = 1:nOHTS_pt
%     
%     for n = 1:N
%         for v = 1:V
%             for m = 1:3
%                 RMSE_Pre_OHTS(n,v,b,m) = (nanmean([JK_error_OHTS(n,v,JK_loc_OHTS{b}-1,m)].^2)); %Obtain MSE for each CV fold
%             end
%         end
%     end
%     
% end

% 
% 
% fprintf('Computing RMSE...\n');
% mLabel = {'MD', 'IOP', 'PSD'};
% 
% RMSE_Table_OHTS = {'# Observations', 'Months Ahead', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%Japan Data
% 
% RMSE_Table_AC = {'# Observations', 'Months Ahead', 'Measure', 'OHTS', 'AGIS/CIGTS*', 'Linear Regression', 'Linear Regression 2', 'NULL'};%AGIS/CIGTS data
% 
% %OHTS testing set
% RMSE_OHTS_OHTS = NaN([N, V, 3]);
% RMSE_AC_OHTS = NaN([N, V, 3]);
% RMSE_NULL_OHTS = NaN([N, V, 3]);
% RMSE_LR_OHTS = NaN([N, V, 3]);
% RMSE_LR2_OHTS = NaN([N, V, 3]);
%  %AC Testing set
% RMSE_OHTS_AC = NaN([N, V, 3]);
% RMSE_AC_AC = NaN([N, V, 3]);
% RMSE_NULL_AC = NaN([N, V, 3]);
% RMSE_LR_AC = NaN([N, V, 3]);
% RMSE_LR2_AC = NaN([N, V, 3]);
% 
% for n = 1:N
%     for m = 1:3
%         for v = 1:V
%             %RMSE_JP_JP(n,v,m) = sqrt(nanmean(JK_error_JP(n,v,:,m).^2));
%             RMSE_OHTS_OHTS(n,v,m) = sqrt(nanmean(squeeze(RMSE_Pre_OHTS(n,v,:,m))));
%             RMSE_AC_AC(n,v,m) = sqrt(nanmean(JK_error_AC(n,v,:,m).^2));
% 
%             RMSE_AC_OHTS(n,v,m) = sqrt(nanmean(error_AC(n,v,:,m).^2));
%             RMSE_OHTS_AC(n,v,m) = sqrt(nanmean(error_OHTS(n,v,:,m).^2));
%             
%             RMSE_LR_OHTS(n,v,m) = sqrt(nanmean(error_LR_OHTS(n,v,:,m).^2));
%             RMSE_LR_AC(n,v,m) = sqrt(nanmean(error_LR_AC(n,v,:,m).^2));
%  
%             RMSE_LR2_OHTS(n,v,m) = sqrt(nanmean(error_LR2_OHTS(n,v,:,m).^2));
%             RMSE_LR2_AC(n,v,m) = sqrt(nanmean(error_LR2_AC(n,v,:,m).^2));
%             
%             RMSE_NULL_OHTS(n,v,m) = sqrt(nanmean(error_NULL_OHTS(n,v,:,m).^2));
%             RMSE_NULL_AC(n,v,m) = sqrt(nanmean(error_NULL_AC(n,v,:,m).^2));
% 
%             RMSE_Table_OHTS = [RMSE_Table_OHTS; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_OHTS(n,v,m), RMSE_AC_OHTS(n,v,m), ...
%                 RMSE_LR_OHTS(n,v,m), RMSE_LR2_OHTS(n,v,m), RMSE_NULL_OHTS(n,v,m)}];
%             RMSE_Table_AC = [RMSE_Table_AC; {num_obs(n), vis_ahead(v)*6, mLabel{m}, RMSE_OHTS_AC(n,v,m), RMSE_AC_AC(n,v,m), ...
%                 RMSE_LR_AC(n,v,m), RMSE_LR2_AC(n,v,m), RMSE_NULL_AC(n,v,m)}];
%         end
%     end
% end
% 
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

RMSE_Table_OHTS = {'# Observations', 'Months Before Progression', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%OHTS Overall
RMSE_Table_OHTS_prog = {'# Observations', 'Months Before Progression', 'Measure', 'OHTS*', 'AGIS/CIGTS', 'Linear Regression', 'Linear Regression 2', 'NULL'};%OHTS Progressing eye only





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

            RMSE_Table_OHTS = [RMSE_Table_OHTS; {num_obs(n), vis_before(v)*6, mLabel{m}, RMSE_OHTS_OHTS(n,v,m), RMSE_AC_OHTS(n,v,m), ...
                RMSE_LR_OHTS(n,v,m), RMSE_LR2_OHTS(n,v,m), RMSE_NULL_OHTS(n,v,m)}];            
            
            RMSE_Table_OHTS_prog = [RMSE_Table_OHTS_prog; {num_obs(n), vis_before(v)*6, mLabel{m}, RMSE_OHTS_OHTS_prog(n,v,m), RMSE_AC_OHTS_prog(n,v,m), ...
                RMSE_LR_OHTS_prog(n,v,m), RMSE_LR2_OHTS_prog(n,v,m), RMSE_NULL_OHTS_prog(n,v,m)}];

        end
    end
end

%% Save results
fprintf('Saving results to mat file and xlsx files\n')
fname = ['JK_results_AC_OHTS_LR_LR2_null_num_obs',num2str(num_obs),'_PROGRESSION_',date];
save(fname)
% 
% Save RMSE Tables as XLSX 
xlsxfname = ['JK_Tables_AC_OHTS_LR_LR2_null_num_obs',num2str(num_obs),'_first_vis_', num2str(first_visit), '_oldkf_', date,'_v2.xlsx'];
xlswrite(xlsxfname,RMSE_Table_OHTS,'OHTS Data')
xlswrite(xlsxfname,RMSE_Table_OHTS_prog,'OHTS Progressors')
xlswrite(xlsxfname,RMSE_Table_OHTS,'OHTS Non Progressors')
xlswrite(xlsxfname,RMSE_Table_AC,'AGIS CIGTS Data')
