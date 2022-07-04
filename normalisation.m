function norm_data = normalisation(input_data, x)

if isempty(x)
    x = 1:length(input_data);
else
    input_data = input_data(x,:);
end


if size(input_data, 1) < 100
    nf = 102;
else
    nf = 101;
end


 % Time normalised phases; phase 1
 %--------------------------------
norm_data = interp1([1:size(input_data,1)],...
     input_data', [1:(size(input_data,1))/nf:size(input_data,1)], 'spline');


