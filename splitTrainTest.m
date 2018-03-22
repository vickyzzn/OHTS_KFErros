function [ train,test ] = splitTrainTest( data, ratio )
%SPLITTRAINTEST Summary of this function goes here
%   Detailed explanation goes here
train=data(1,:);
test=data(1,:);

rng(1);
for i=2:size(data,1)
if rand()<=ratio
    train=[train;data(i,:)];
else
     test=[test;data(i,:)];
end
end

