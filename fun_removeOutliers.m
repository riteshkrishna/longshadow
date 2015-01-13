% Function to remove outliers from data
function [Data_nooutlier,I] = fun_removeOutliers(Data)

z = 3;
n = size(Data,1);

mean_data = mean(Data);
std_data = std(Data);
I = abs(Data - repmat(mean_data,n,1)) <= z .* repmat(std_data,n,1);

idx = find(all(I,2));

Data_nooutlier = Data(idx,:);
end