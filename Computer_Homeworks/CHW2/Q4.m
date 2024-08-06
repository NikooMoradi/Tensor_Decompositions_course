% Nikoo Moradi
% 400101934
% CHW 2 
% Q4
%% 

clc; clear; close all;

house_dataset=load('house_dataset.mat');
input = house_dataset.houseInputs; % 13 x 506 matrix (13 features and 506 datapoint)
target = house_dataset.houseTargets;

% Calculate the average price
average_price = mean(target);

% Classify each datapoint
% If price is higher than the average, class is 1, otherwise class is 0
class_labels = target > average_price;

% Display the first few classifications to verify
disp('First few classifications:');
disp(class_labels(1:10));
disp(target(1:10));
disp('Mean: ');
disp(average_price);

num_class_1 = sum(class_labels);
num_class_0 = sum(~class_labels);

fprintf('Number of houses in class 1 (price > average): %d\n', num_class_1);
fprintf('Number of houses in class 0 (price <= average): %d\n', num_class_0);

% Optional: Visualize the classification on the first two principal components
% Centering the data
features_centered = (input - mean(input,2)) ./ std(input,0,2);

% Performing SVD
[U, S, V] = svd(features_centered, 'econ');

% Projecting the features onto the first two principal components
projected_features = U' * features_centered ;



% Display the principal component vectors (loadings)
disp('Principal Component Vectors (Loadings):');
disp(U);

% Plot the loadings for the first two principal components
figure;
num = 1;
for i = 1:13
    subplot(3,5,num);
    num = num + 1;
    bar(U(:,i));
    xlabel(sprintf('Feature %d', i));
    ylabel('Loading');
    title(sprintf('%d th Principal Component',i));
end

singular_values = diag(S);

% Total variance
total_variance = sum(singular_values .^ 2);

% Variance captured by each principal component
explained_variance = (singular_values .^ 2) / total_variance * 100;

% Display the explained variance
disp('Explained Variance by each Principal Component (%):');
disp(explained_variance);

% Plot the explained variance
figure;
bar(explained_variance);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Variance Captured by Each Principal Component');
grid on;


% Plot the classification
figure('Name', 'Classification based on House Prices');

num = 1;

for i= 1:5
    for j = i+1 : 5
        
        
        subplot(5,2,num);
        hold on
        scatter(projected_features(i, class_labels), projected_features(j, class_labels),'r');
        scatter(projected_features(i, ~class_labels), projected_features(j, ~class_labels),'b');
        
        xlabel(sprintf('PC %d', i));
        ylabel(sprintf('PC %d', j));
        legend('Class 1 (price > average)', 'Class 0 (price <= average)');
        grid on;
        num = num + 1;
        
    end
end



