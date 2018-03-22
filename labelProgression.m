function [ data ] = labelProgression( data )
%LABELPROGRESSION Summary of this function goes here
%   Detailed explanation goes here
numPat=size(data,1)-1;
data{1,end+1}='MDProg';
for i=2:numPat+1
    baseLine=data{i,4};
    readings=data{i,3};
    numVis=size(baseLine,2);
    MDprog=zeros(1,numVis);
    k=1;
    while k<numVis
        diff = (baseLine(1,:) - readings(1,:))>2; % drop in MD
        if diff(k) == 1
            if sum(diff(k:end))> 0 % confirmation
                MDprog(k) = 1;
                if k<numVis-1
                    baseLine = (readings(1:3,k+1)+readings(1:3,k+2))*0.5*ones(1,numVis);
                else
                    baseLine = (readings(1:3,k)+readings(1:3,k+1))*0.5*ones(1,numVis);
                end
                data{i,4}(:,k+1:end)=baseLine(:,k+1:end);                
            end
        end
        k = k + 1;
        
    end
    data{i,end}=MDprog;
end

fprintf('Total num of progressions: %d\n',sum([data{2:end,end}]));