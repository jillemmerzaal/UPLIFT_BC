%this subject is a right arm, flip left arm data to right arm.
%50 Hz data
%coding is 1=>nonfunctional, 2=>functional, 0=>unknown

load('model.mat');

data2=readmatrix('RSE_026_NP.csv');

b=200;
points=floor(length(data)/b)*b;
data_1=[0 0 0];
data_2=[0 0 0];
k=1;
j=1;
for i=1:b:points
    temp0=length(find(data(i:i+b-1,4)==0));
    temp1=length(find(data(i:i+b-1,4)==1));
    temp2=length(find(data(i:i+b-1,4)==2));

    if temp0>b/2
    elseif temp1>b*3/4
        data_1=[data_1; data(i:i+b-1,1:3)];
    elseif temp2>b*3/4
        data_2=[data_2; data(i:i+b-1,1:3)];
    end
end

data_1=data_1(2:length(data_1),:);
data_2=data_2(2:length(data_2),:);

feature1=featurecalc1(data_1,b);  %nonfunctional
feature2=featurecalc1(data_2,b);  %functional
%recoding to 0,1 for nonfunctional, functional
output=[zeros(length(feature1(:,1)),1);ones(length(feature2(:,1)),1)];
feature=[feature1;feature2];

yyfit = trainedModel.predictFcn(feature);
CC(:,:)=confusionmat(output,yyfit);
Accuracy=(CC(1,1)+CC(2,2))/(CC(1,1)+CC(1,2)+CC(2,1)+CC(2,2));
