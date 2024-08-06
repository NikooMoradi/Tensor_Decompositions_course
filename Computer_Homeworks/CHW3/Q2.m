% Nikoo Moradi
% 400101934
% CHW 2 
% Q2
%%
clear all; clc; close all;

load('swimmer.mat');
L = length(A);

Y = zeros(L, numel(A{1}));
for i = 1:L
    Y(i, :) = reshape(A{i}, 1, []);
end

% Finding the minimum number of needed sources for both algorithms
J_range = 1:50;
errors_mult = zeros(1, length(J_range));
errors_als = zeros(1, length(J_range));

for J = J_range
    % Mult
    [B_mult, C_mult] = nnmf(Y, J, 'algorithm', 'mult');
    E_mult = Y - B_mult * C_mult;
    errors_mult(J) = norm(E_mult, 'fro');
    
    % ALS 
    [B_als, C_als] = nnmf(Y, J, 'algorithm', 'als');
    E_als = Y - B_als * C_als;
    errors_als(J) = norm(E_als, 'fro');
end

% Plotting error versus number of sources for both algorithms
figure;
plot(J_range, errors_mult, 'b', J_range, errors_als, 'r');
% p(2).LineWidth = 2;
title('Error vs Number of Sources');
grid on;
xlabel('Number of sources (J)');
ylabel('Frobenius norm of E');
legend('Multiplicative algorithm', 'ALS algorithm');

% Plotting Y 
figure;
imagesc(Y);
title('Y');

% Best number of sources for both algorithms
best_J_mult = 20; 
best_J_als = 1;   

% Mult with best J
[B_mult, C_mult] = nnmf(Y, best_J_mult, 'algorithm', 'mult');
E_mult = Y - B_mult * C_mult;

% ALS with best J
[B_als, C_als] = nnmf(Y, best_J_als, 'algorithm', 'als');
E_als = Y - B_als * C_als;

% Plotting B 
figure;
subplot(1, 2, 1);
imagesc(B_mult);
title('B (Multiplicative Algorithm)');

subplot(1, 2, 2);
imagesc(B_als);
title('B (ALS Algorithm)');

% Plotting C for Mult
figure;
for j = 1:best_J_mult
    subplot(4,5,j);
    imagesc(reshape(C_mult(j, :), 9, 14));
    title(['C (Multiplicative) ', num2str(j)]);
end

% Plotting C for ALS 
for j = 1:best_J_als
    figure;
    imagesc(reshape(C_als(j, :), 9, 14));
    title(['C (ALS) ', num2str(j)]);
end

% Plotting E 
figure;
subplot(1, 2, 1);
imagesc(E_mult);
title('Error (Multiplicative Algorithm)');

subplot(1, 2, 2);
imagesc(E_als);
title('Error (ALS Algorithm)');

