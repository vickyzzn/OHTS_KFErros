function [DATA, exclude_list] = OHTS_exclusion(data, thresholds)

%This function throws out readings based on MD, IOP, and PSD thresholds
%Thresholds lists the thresholds for MD, IOP, PSD as:
% [MDLB, MDUB; IOPLB, IOPUB; PSDLB, PSDUB]

numPat=size(data,1)-1;
outPatCount=2;
DATA=data;

excludeCount=0;

exclude_list = cell(1,4);
exclude_list(1,:) = {'Patient ID', 'Reading', 'Value','Visit'};

measures = {'MD', 'IOP', 'PSD'};

for i = 2:numPat+1
    readings = data{i,3}; %old readings
    new_readings = data{i,3}; %new readings
    [~, num_readings] = size(readings);
    visits = data{i,2};
    for j = 1:3
        for n = 1:num_readings
            if readings(j,n) <= thresholds(j,1) | readings(j,n) >= thresholds(j,2)
                new_readings(j,n) = NaN;
                exclude_list = [exclude_list; {data{i,1}, measures{j}, readings(j,n), data{i,2}(n)}];
            end
            
        end
    end
    DATA{i,3} = new_readings;
end



end
