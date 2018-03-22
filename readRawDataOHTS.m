function  [ DATA ] = readRawDataOHTS(filename)

[numbers,strings]=xlsread(filename);
ids=unique(strings(:,1));
ids=ids(1:end-1);
numPat=length(ids);
DATA=cell(numPat+1,7); 
DATA(1,:)={'ID','VisNum','RawReadings','age','race','sex','MDprog'};
%race sex coding information
%sex – sex (1 = male, 2 = female)
%race – race (1 = white, 2 = black, 3 others)
for i=1:numPat
    readingsIndex=find(ismember(strings(:,1),ids{i}));
    DATA{i+1,1}=ids{i}; %id
    DATA{i+1,2}=numbers(readingsIndex-1,1)'; % timing
    DATA{i+1,3}=numbers(readingsIndex-1,2:4)'; %3d readings
    DATA{i+1,4}=numbers(readingsIndex-1,5)';% age
    DATA{i+1,5}=numbers(readingsIndex-1,6)';% race
    DATA{i+1,6}=numbers(readingsIndex-1,7)';% sex
    DATA{i+1,7}=numbers(readingsIndex-1,8)';%progression
end
DATA=frontBackTrim(DATA);
end