%test sampen
clear all; close all; clc

df = [9 4 8 1 2 4 5 8 7 7 0 2 1 5 9 8 8 0 1 0 4 2 10 0 0 1];

hydra_sampen = sampen(df, 2, 1);

figure; plot(df)

N = length(df);
m = 2;
r = 1;

matches = NaN(m,N);
for i = 1:m
    matches(i,1:N+1-i) = df(i:end);
end
matches = matches';

for i = 1:size(matches,1)-1
    lower_bounds = matches(i,:) - r;
    upper_bounds = matches(i,:) + r;
%     lower_row2 = matches(i,2) - r;
%     upper_row2 = matches(i,2) + r;

    
%     match_is_in_range(i) = sum(matches(:,1) >= lower_row1 & ...
%         matches(:,1) <= upper_row1 & ...
%         matches(:,2) >= lower_row2 & ...
%         matches(:,2) <= upper_row2)-1;

    match_is_in_range = sum(matches >= lower_bounds & matches <= upper_bounds,2);
    matches_total(i) = sum(match_is_in_range == m) - 1; 
end

B_i = (sum(matches_total) / (N-m)) / (N-m);

clear matches matches_total

matches = NaN(m+1,N);
for i = 1:1:m+1
    matches(i,1:N+1-i) = df(i:end);
end

matches = matches';

for i = 1:size(matches,1)-2
    lower_bounds = matches(i,:) - r;
    upper_bounds = matches(i,:) + r;
  
    matches_in_range = sum(matches >= lower_bounds & matches <= upper_bounds,2);
    matches_total(i) = sum(matches_in_range == (m+1))-1;
    
end


A_i = (sum(matches_total) / (N-(m+1)))/(N-m);

A = (((N-m-1) * (N-m)) / 2) * A_i;
B = (((N-m-1) * (N-m)) / 2) * B_i;

SampEn = -log(A/B);
