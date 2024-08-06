% Nikoo Moradi
% 400101934
% CHW 2 
% Q3
%% A

clc; clear; close all;

PCAdata = load('PCAdata.mat');
PCAdata = PCAdata.PCAdata;
PCAdata = (PCAdata - mean(PCAdata,2)) ./ std(PCAdata,0,2);


X = PCAdata(1,:);
Y = PCAdata(2,:);
Z = PCAdata(3,:);


[U,S,V] = svd(PCAdata, 'econ');

figure;
scatter3(X,Y,Z);
hold on
p = plot3(5*[-U(1,1),U(1,1)],5*[-U(1,2),U(1,2)],5*[-U(1,3),U(1,3)],'r');
p.LineWidth = 2;
p = plot3(3*[-U(2,1),U(2,1)],3*[-U(2,2),U(2,2)],3*[-U(2,3),U(2,3)],'g');
p.LineWidth = 2;
p = plot3(7*[-U(3,1),U(3,1)],7*[-U(3,2),U(3,2)],7*[-U(3,3),U(3,3)],'c');
p.LineWidth = 2;

xlabel('X');
ylabel('Y');
zlabel('Z');
legend('datapoints','U1','U2','U3');
title('Normalized data');
%% B


% Scatter of the data
total_variance = sum(diag(S) .^ 2);

% Data variance explained by each principal component
data_variance_explained_svd = (diag(S) .^ 2) / total_variance * 100

% Plot the explained variance
figure;
bar(data_variance_explained_svd);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Variance Captured by Each Principal Component _ using SVD');
grid on;


% New directions (principal components)
new_directions_svd = U

figure;
num = 1;
for i = 1:3
    subplot(3,1,num);
    num = num + 1;
    bar(U(:,i));
    xlabel(sprintf('PC %d', i));
    ylabel('Loading');
    title(sprintf('%d th Principal Component _ using SVD',i));
end


% Whitened data (data projected onto principal components)
whitened_data = U' * PCAdata;

figure;
hold on
grid on
scatter3(whitened_data(1,:),whitened_data(2,:),whitened_data(3,:));
title('whitened data using SVD');
xlabel('U1');
ylabel('U2');
zlabel('U3');

%% C

[coeff, score, ~, ~, explained] = pca(PCAdata');

% Data variance explained by each principal component
data_variance_explained_pca = explained

data_var_pca = var(PCAdata,0,2);

% Plot the explained variance
figure;
bar(data_variance_explained_pca);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Variance Captured by Each Principal Component _ usnig PCA');
grid on;

% New directions (principal components)
new_directions_pca = coeff

figure;
num = 1;
for i = 1:3
    subplot(3,1,num);
    num = num + 1;
    bar(coeff(:,i));
    xlabel(sprintf('PC %d', i));
    ylabel('Loading');
    title(sprintf('%d th Principal Component _ using PCA',i));
end

% Whitened data (data projected onto principal components)
whitened_data_pca = score';


figure;
hold on
grid on
scatter3(whitened_data_pca(1,:),whitened_data_pca(2,:),whitened_data_pca(3,:));
title('whitened data using PCA');
xlabel('coeff1');
ylabel('coeff2');
zlabel('coeff3');
