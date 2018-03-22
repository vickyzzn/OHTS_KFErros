function [ERRORS] = kalman_error_new(warmup, test, A, C, Q, R, init_x, init_V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the error in kalman filter predictions of MD, PSD,
% IOP.

%Inputs:
%Warmup is the number of time periods for warmup (6 month increments)
%test is the entire testing set of data. It is assumed that the first row
%   has text which describes what belongs to that column.
%A, C, Q, R, init_X, and init_V are all Kalman filter parameters

%outputs:
%provides the array of errors

id_col = find(strcmp(test(1,:), 'ID') == 1, 1); %find the column w/ IDs
readings_col = find(strcmp(test(1,:), 'InterpolatedReadings') == 1, 1); %find the columns w/ readings

data = [test(2:end,id_col), test(2:end,readings_col)]; %take just the IDs and readings



[N, ~] = size(data); %number of patients

%create errors table
ERRORS{1,1} = 'ID';
ERRORS{1,2} = 'Error table';


for n = 1:N
    y = data{n,2}; %get data
    T = size(y,2); %number of visits
    num_measurements = size(y,1); %number of measurements
    error = cell([1,4]);
    error{1,1} = 'Num obs'; %create a table which we can use to store errors for patient n
    error{1,2} = 'Visits ahead (6 mos)';
    error{1,3} = 'Errors';
    error{1,4} = 'First visit';
    step = 2; %start from second row
    if T > warmup
        for t = warmup:T %number of observations
            for j = 0:T-t-1
                x = [y(:,j+1:(j+t)), NaN([num_measurements, T - (t+j)])];%predict all visits from t until the end of the horizon
                [z,~,~,~]=kalman_filter(x,A,C,Q,R,init_x,init_V);
                for s = j+t+1:T %months ahead
                    err_vec = y(:,s) - z(:,s-j); %observed - predicted
                    error{step, 1} = t; %record number of observations
                    error{step, 2} = s-(t+j); %record number of months ahead
                    error{step, 3} = err_vec; %record the errors
                    error{step, 4} = j + 1; %record the first visit
                    step = step + 1;
                end
            end
        end
        ERRORS{n, 1} = data{n,1}; %record patient ID
        ERRORS{n, 2} = error; %record error table
    else
        ERRORS{n, 1} = data{n, 1}; %record patient ID
        ERRORS{n, 2} = 'Not enough observations';
        fprintf(['Not enough observations for patient ', data{n,1},'\n'])
    end
end

ERRORS = [{'ID', 'Error Table'}; ERRORS];
end