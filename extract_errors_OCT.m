function [ID_col, error_array] = extract_errors_OCT(ERRORS, num_obs, vis_ahead, first_visit, type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       For a fixed number of observations and visits ahead, we extract the
%       error of our kalman filter predictions
%INPUTS: %   type: 0 corresponds to MD, IOP, PSD
%         1 corresponds to MD, IOP, PSD, OCT
%         2 corresponds to MD, IOP, PSD, OCT, DH
%ERRORS is a table which has patient IDs in the first column and the entire
%   Error Table saved in the second column. Obtained using kalman_error
%   function
%num_obs is the number of observations before making kalman filter
%   predictions
%vis_ahead is the number of visits ahead we want to predict
%first_visit is the visit number of the first visit used to generate errors
%   a value of 0 implies take all
%OUTPUTS:
%ID_col is the column of patient IDs and error_array is a (num_patients)x9
%array which specifies the error in kalman filter prediction for each
%patient and measure



errors = ERRORS(2:end,:); %take everything except the current headers

N = size(errors,1); %number of patients

if type == 0
    numrows = 9;
elseif type == 1
    numrows = 10;
else
    numrows = 11;
end



counter = 1; %row counter
for n = 1:N
    data = errors{n,2};   
 
    if ~ischar(data)
        
        data = data(2:end,:);

        obs_col = [data{:,1}]; %number of observations
        vis_col = [data{:,2}]; %number of visits ahead
        first_col = [data{:,4}]; %first visit number

        if length(obs_col) > 1    
            
            if first_visit ~= 0
                row_num = find(obs_col == num_obs & vis_col == vis_ahead & first_col == first_visit);
            else
                row_num = find(obs_col == num_obs & vis_col == vis_ahead);
            end
        else %the double && is used for scalar conditions (only one row)
            if first_visit ~= 0
                row_num = find(obs_col == num_obs && vis_col == vis_ahead && first_col == first_visit);
            else
                row_num = find(obs_col == num_obs && vis_col == vis_ahead);
            end
        end

        if ~isempty(row_num) %if we were able to find such a row
            for r =1:length(row_num)
                error_row = data{row_num(r),3}';
                error_array(counter,:) = error_row;
                ID_col{counter,1} = errors{n,1}; %record patient id
                counter = counter + 1;
            end
        else
            error_array(counter,:) = NaN([1, numrows]);
            ID_col{counter,1} = errors{n,1}; %record patient id
            fprintf(['Not enough observations for patient ', errors{n,1},'\n'])
            counter = counter + 1;
        end
    else
        error_array(counter,:) =  NaN([1, numrows]);
        ID_col{counter,1} = errors{n,1}; %record patient id
        fprintf(['Not enough observations for patient ', errors{n,1}, '\n'])
        counter = counter + 1;
        
    end
end



end