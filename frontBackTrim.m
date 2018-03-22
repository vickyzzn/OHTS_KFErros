function [ data ] = frontBackTrim( data )
%FRONTTRIM Summary of this function goes here
%   Detailed explanation goes here
for i=2:size(data,1)
    while data{i,2}(2)- data{i,2}(1)>24
        for j=2:size(data,2)
            if j==3
                data{i,3}=data{i,3}(:,2:end);
                continue;
            end
            data{i,j}=data{i,j}(2:end);
        end
        if  length(data{i,2})==1
            break;
        end
    end
    if  length(data{i,2})==1
        break;
    end
    
    while data{i,2}(end)- data{i,2}(end-1)>24
        for j=2:size(data,2)
            if j==3
                data{i,3}=data{i,3}(:,2:end);
                continue;
            end
            data{i,j}=data{i,j}(2:end);
        end
        if  length(data{i,2})==1
            break;
        end
    end
    
    
    
end

end

