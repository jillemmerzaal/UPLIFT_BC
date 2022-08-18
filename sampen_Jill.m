function [SampEn] = sampen_Jill(df, m, r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% figure; plot(df)
N = length(df);

if isempty(m)
    m = 2;
end
if isempty(r)
    r = 1;
end

%% calculate B_i
matches = NaN(m,N);
for i = 1:m
    matches(i,1:N+1-i) = df(i:end);
end
matches = matches';

for i = 1:size(matches,1)-1
    lower_bounds = matches(i,:) - r;
    upper_bounds = matches(i,:) + r;

    match_is_in_range = sum(matches >= lower_bounds & matches <= upper_bounds,2);
    matches_total(i) = sum(match_is_in_range == m) - 1; 
end
B_i = (sum(matches_total) / (N-m)) / (N-m);

clear matches matches_total

%% calculate A_i

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

%% calculate sample entropy
A = (((N-m-1) * (N-m)) / 2) * A_i;
B = (((N-m-1) * (N-m)) / 2) * B_i;

SampEn = -log(A/B);
end