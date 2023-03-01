function [data_1,features, yyfit] = MLmodel(data, model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(model)


b=200;
points=floor(length(data)/b)*b;
data_1=[0 0 0];

for i = 1:b:points
    data_1=[data_1; data(i:i+b-1,1:3)];
end

data_1=data_1(2:length(data_1),:);

features=featurecalc1(data_1,b);

%% prediction
yyfit = trainedModel.predictFcn(feature);




end