% Nikoo Moradi
% 400101934
% CHW 4
%% Q1 _ ALS algorithm for CP decomposition 
clc;clear;close all;

addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TensorLab'); 
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TMSE.m'); 

if ~exist('cpd', 'file')
    error('Tensorlab is not installed or not added to the MATLAB path');
end


TMSE_func = @TMSE; 

I = 5; J = 5; K = 5;
R = 3;

U1_0 = randn(I, R);
U2_0 = randn(J, R);
U3_0 = randn(K, R);
U0 = {U1_0, U2_0, U3_0};
 
T = randn(I, J, K);
% T = zeros(I, J, K);
% for r = 1:R
%     T = T + outer_product(U1_0(:,r), U2_0(:,r), U3_0(:,r));
% end

U_my_cp_als = cp_als(T, U1_0, U2_0, U3_0, 100,1e-5);

U_my_cp_als1 = U_my_cp_als{1};
U_my_cp_als2 = U_my_cp_als{2};
U_my_cp_als3 = U_my_cp_als{3};

T_hat = zeros(I, J, K);
for i = [1: I]
    for t = [1: J]
        for q = [1: K]
            s = 0;
            for j = [1: R]
                s = s + U_my_cp_als1(i, j) * U_my_cp_als2(t, j) * U_my_cp_als3(q, j);
            end
            T_hat(i, t, q) = s;
        end
    end
end

% TMSE_func(U_my_cp_als,U0)
error = TMSE_func(U_my_cp_als,U0);
disp(['TMSE = ', num2str(error)]);

%% Q2 _ A and B

clc; clear; close all;

addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TensorLab'); 
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TMSE.m'); 


TMSE_func = @TMSE; 

I = 5; J = 5; K = 5;
R = 3; 
num_tensors = 30;
snr_values = [0, 20, 40, 60];
TMSE_results = zeros(num_tensors, 4, 4);
TMSE_results_with_HOSVD = zeros(num_tensors, 4, 4);
count = 0;

for snr = snr_values
    count = count + 1
    
    for n = 1:num_tensors
        U1_org = normrnd(0, 1 ,I, R);
        U2_org = normrnd(0, 1 ,J, R);
        U3_org = normrnd(0, 1 ,K, R);
        U_org = {U1_org, U2_org, U3_org};
        
        
        T = zeros(I, J, K);
        for r = 1:R
            T = T + outer_product(U1_org(:,r), U2_org(:,r), U3_org(:,r));
        end
        
        N = normrnd(0,1,size(T));
        T_noisy = T + (10^(-snr/20)) * norm_fro(T) / norm_fro(N) * N;
        
        U1_0 = randn(I, R);
        U2_0 = randn(J, R);
        U3_0 = randn(K, R);
        U0 = {U1_0, U2_0, U3_0};
        
        U_cpd_als = cpd_als(T_noisy, U0);
        
        U_cpd_minf = cpd_minf(T_noisy, U0);
        
        U_cpd3_sd = cpd3_sd(T_noisy, U0);
        
        U_my_cp_als = cp_als(T_noisy, U1_0, U2_0, U3_0, 50, 1e-5);
        
        
        
        U0 = HOSVD(T_noisy);
        
        U_cpd_als_HOSVD = cpd_als(T_noisy, U0);
        
        U_cpd_minf_HOSVD = cpd_minf(T_noisy, U0);
        
        U_cpd3_sd_HOSVD = cpd3_sd(T_noisy, U0);
        
        U_my_cp_als_HOSVD = cp_als(T_noisy, U0{1}, U0{2}, U0{3}, 50, 1e-5);
        
        TMSE_results_with_HOSVD(n, count, 1) = TMSE_func(U_org, U_cpd_als_HOSVD);
        TMSE_results_with_HOSVD(n, count, 2) = TMSE_func(U_org, U_cpd_minf_HOSVD);
        TMSE_results_with_HOSVD(n, count, 3) = TMSE_func(U_org, U_cpd3_sd_HOSVD);
        TMSE_results_with_HOSVD(n, count, 4) = TMSE_func(U_org, U_my_cp_als_HOSVD);
        
        
        TMSE_results(n, count, 1) = TMSE_func(U_org, U_cpd_als);
        TMSE_results(n, count, 2) = TMSE_func(U_org, U_cpd_minf);
        TMSE_results(n, count, 3) = TMSE_func(U_org, U_cpd3_sd);
        TMSE_results(n, count, 4) = TMSE_func(U_org, U_my_cp_als);
    end

    fprintf('\nInitialization randomly\n');
    fprintf('SNR: %d dB\n', snr);
    fprintf('CPD: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, count, 1)), std(TMSE_results(:, count,1)));
    fprintf('CPD_MINF: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, count,2)), std(TMSE_results(:, count,2)));
    fprintf('CPD3_SD: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, count,3)), std(TMSE_results(:, count,3)));
    fprintf('My CP ALS: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, count,4)), std(TMSE_results(:, count,4)));
    
    
    fprintf('\n\nInitialization with HOSVD\n');
    fprintf('SNR: %d dB\n', snr);
    fprintf('CPD: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, count, 1)), std(TMSE_results_with_HOSVD(:, count,1)));
    fprintf('CPD_MINF: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, count,2)), std(TMSE_results_with_HOSVD(:, count,2)));
    fprintf('CPD3_SD: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, count,3)), std(TMSE_results_with_HOSVD(:, count,3)));
    fprintf('My CP ALS: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, count,4)), std(TMSE_results_with_HOSVD(:, count,4)));
end


%% Plotting results of parts A and B

CP_Mean = squeeze(mean(TMSE_results, 1));

figure;
hold on;
p(1)=plot(snr_values, CP_Mean(:, 1), 'b');
p(1).LineWidth = 2;
p(2)=plot(snr_values, CP_Mean(:, 2), 'r');
p(2).LineWidth = 2;
p(3)=plot(snr_values,CP_Mean(:, 3), 'g');
p(3).LineWidth = 2;
p(4)=plot(snr_values,CP_Mean(:, 4));
p(4).LineWidth = 2;
xlabel('SNR');
ylabel('TMSE Error');
title('Error vs SNR  (Initialization : Randomly)');
legend('ALS with TensorLab', 'Minf with Tensor Lab', 'SD with Tensor Lab', 'My CP ALS');
grid on;

HOSVD_Mean = squeeze(mean(TMSE_results_with_HOSVD, 1));


figure;
hold on;
p(1)=plot(snr_values, HOSVD_Mean(:, 1), 'b');
p(1).LineWidth = 2;
p(2)=plot(snr_values, HOSVD_Mean(:, 2), 'r');
p(2).LineWidth = 2;
p(3)=plot(snr_values,HOSVD_Mean(:, 3), 'g');
p(3).LineWidth = 2;
p(4)=plot(snr_values,HOSVD_Mean(:, 4));
p(4).LineWidth = 2;
xlabel('SNR');
ylabel('TMSE Error');
title('Error vs SNR (Initialization : with HOSVD)');
legend('ALS with TensorLab', 'Minf with Tensor Lab', 'SD with Tensor Lab', 'My CP ALS');
grid on;

%% Q2 _ C and D
clc; clear; close all;

addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TensorLab'); 
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TMSE.m'); 


TMSE_func = @TMSE;

I = 5; J = 5; K = 5;
R = 3; 
num_tensors = 30;
snr_values = [0, 20, 40, 60];

TMSE_results_with_HOSVD = zeros(num_tensors, length(snr_values), 4);
TMSE_results = zeros(num_tensors, length(snr_values), 4);

for snr_idx = 1:length(snr_values)
    snr = snr_values(snr_idx);
    
    for n = 1:num_tensors
        u1_1 = randn(I, 1);
        v1_2 = randn(I, 1);
        u1_2 = u1_1 + 0.5 * v1_2;
        u1_3 = randn(I, 1);
        U1_org = [u1_1, u1_2, u1_3];
        
        U3_org = randn(K, R);
        
        u2_1 = randn(J, 1);
        v2_2 = randn(J, 1);
        u2_2 = u2_1 + 0.5 * v2_2;
        u2_3 = randn(J, 1);
        U2_org = [u2_1, u2_2, u2_3];
        
        U_org = {U1_org, U2_org, U3_org};
        T = zeros(I, J, K);
        for r = 1:R
            T = T + outer_product(U1_org(:,r), U2_org(:,r), U3_org(:,r));
        end
        
        N = randn(size(T));
        T_noisy = T + (10^(-snr/20)) * norm_fro(T) / norm_fro(N) * N;
        
        U1_0 = randn(I, R);
        U2_0 = randn(J, R);
        U3_0 = randn(K, R);
        U0 = {U1_0, U2_0, U3_0};
        
        U_cpd_als = cpd_als(T_noisy, U0);
        
        U_cpd_minf = cpd_minf(T_noisy, U0);
        
        U_cpd3_sd = cpd3_sd(T_noisy, U0);
        
        U_my_cp_als = cp_als(T_noisy, U1_0, U2_0, U3_0, 50, 1e-5);
        
        TMSE_results(n, snr_idx, 1) = TMSE_func(U_org, U_cpd_als);
        TMSE_results(n, snr_idx, 2) = TMSE_func(U_org, U_cpd_minf);
        TMSE_results(n, snr_idx, 3) = TMSE_func(U_org, U_cpd3_sd);
        TMSE_results(n, snr_idx, 4) = TMSE_func(U_org, U_my_cp_als);
        
        U0 = HOSVD(T_noisy);
        
        U_cpd_als_HOSVD = cpd_als(T_noisy, U0);
        
        U_cpd_minf_HOSVD = cpd_minf(T_noisy, U0);
        
        U_cpd3_sd_HOSVD = cpd3_sd(T_noisy, U0);
        
        U_my_cp_als_HOSVD = cp_als(T_noisy, U0{1}, U0{2}, U0{3}, 50, 1e-5);
        
        TMSE_results_with_HOSVD(n, snr_idx, 1) = TMSE_func(U_org, U_cpd_als_HOSVD);
        TMSE_results_with_HOSVD(n, snr_idx, 2) = TMSE_func(U_org, U_cpd_minf_HOSVD);
        TMSE_results_with_HOSVD(n, snr_idx, 3) = TMSE_func(U_org, U_cpd3_sd_HOSVD);
        TMSE_results_with_HOSVD(n, snr_idx, 4) = TMSE_func(U_org, U_my_cp_als_HOSVD);

    end
end

for snr_idx = 1:length(snr_values)
    snr = snr_values(snr_idx);
    
    fprintf('\nInitialization randomly dependently\n');
    fprintf('SNR: %d dB\n', snr);
    fprintf('CPD: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, snr_idx, 1)), std(TMSE_results(:, snr_idx,1)));
    fprintf('CPD_MINF: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, snr_idx,2)), std(TMSE_results(:, snr_idx,2)));
    fprintf('CPD3_SD: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, snr_idx,3)), std(TMSE_results(:, snr_idx,3)));
    fprintf('My CP ALS: TMSE = %.2e ± %.2e\n', mean(TMSE_results(:, snr_idx,4)), std(TMSE_results(:, snr_idx,4)));
    
    
    fprintf('\n\nInitialization with HOSVD\n');
    fprintf('SNR: %d dB\n', snr);
    fprintf('CPD: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, snr_idx, 1)), std(TMSE_results_with_HOSVD(:, snr_idx,1)));
    fprintf('CPD_MINF: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, snr_idx,2)), std(TMSE_results_with_HOSVD(:, snr_idx,2)));
    fprintf('CPD3_SD: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, snr_idx,3)), std(TMSE_results_with_HOSVD(:, snr_idx,3)));
    fprintf('My CP ALS: TMSE = %.2e ± %.2e\n', mean(TMSE_results_with_HOSVD(:, snr_idx,4)), std(TMSE_results_with_HOSVD(:, snr_idx,4)));
end


%% Plotting results of parts C and D

CP_Mean = squeeze(mean(TMSE_results, 1));

figure;
hold on;
p(1)=plot(snr_values, CP_Mean(:, 1), 'b');
p(1).LineWidth = 2;
p(2)=plot(snr_values, CP_Mean(:, 2), 'r');
p(2).LineWidth = 2;
p(3)=plot(snr_values,CP_Mean(:, 3), 'g');
p(3).LineWidth = 2;
p(4)=plot(snr_values,CP_Mean(:, 4));
p(4).LineWidth = 2;
xlabel('SNR');
ylabel('TMSE Error');
title('Error vs SNR  (Initialization : Randomly Depentently)');
legend('ALS with TensorLab', 'Minf with Tensor Lab', 'SD with Tensor Lab', 'My CP ALS');
grid on;

HOSVD_Mean = squeeze(mean(TMSE_results_with_HOSVD, 1));


figure;
hold on;
p(1)=plot(snr_values, HOSVD_Mean(:, 1), 'b');
p(1).LineWidth = 2;
p(2)=plot(snr_values, HOSVD_Mean(:, 2), 'r');
p(2).LineWidth = 2;
p(3)=plot(snr_values,HOSVD_Mean(:, 3), 'g');
p(3).LineWidth = 2;
p(4)=plot(snr_values,HOSVD_Mean(:, 4));
p(4).LineWidth = 2;
xlabel('SNR');
ylabel('TMSE Error');
title('Error vs SNR (Initialization : with HOSVD)');
legend('ALS with TensorLab', 'Minf with TensorLab', 'SD with TensorLab', 'My CP ALS');
grid on;

%% Q3 _ A

clc; clear; close all;

addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TensorLab'); 
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/4/SimHW4/TMSE.m'); 

TMSE_func = @TMSE;

% Load the data
data = load('amino.mat');
T_temp = data.X; 
I = data.DimX(1);
J = data.DimX(2);
K = data.DimX(3);

T = [];
for i = [1: 1: K]
    T(:, :, i) = T_temp(:, [(i - 1) * J + 1: i * J]);
end

ranks = [2, 3, 4, 5];

results = struct();
algorithms = {'ALS__with__TensorLab', 'Minf__with__TensorLab', 'SD__with__TensorLab', 'My__CP__ALS'};

for r = ranks
    
    % Initialization for decomposition
%         U1_0 = randn(I, r);
%         U2_0 = randn(J, r);
%         U3_0 = randn(K, r);
%         U0 = {U1_0, U2_0, U3_0};
        U0 = cpd_rnd(size(T), r);
        
        
        U_cpd_als = cpd_als(T, U0);
        results.("ALS__with__TensorLab")(r).U = U_cpd_als;
        
        U_cpd_minf = cpd_minf(T, U0);
        results.("Minf__with__TensorLab")(r).U = U_cpd_minf;
        
        U_cpd3_sd = cpd3_sd(T, U0);
        results.("SD__with__TensorLab")(r).U = U_cpd3_sd;
        
        U_my_cp_als = cp_als(T, U0{1}, U0{2}, U0{3}, 50, 1e-5);
        results.("My__CP__ALS")(r).U = U_my_cp_als;
        
    
end

% Plotting results
for r = ranks
    figure('Name', ['Rank ' num2str(r)], 'NumberTitle', 'off');
    
    for alg_idx = 1:length(algorithms)
        algorithm = algorithms{alg_idx}
        U = results.(algorithm)(r).U;
       
        subplot(2, 2, alg_idx);
        hold on;
        grid on;
        for mode = 1:r
          
            plot([250: 450],U{2}(:,mode));
            hold on;
            title((algorithm));
        end
        hold off;
    end
end

%% Q3 _ B

algorithms = {'ALS__with__TensorLab', 'Minf__with__TensorLab', 'SD__with__TensorLab', 'My__CP__ALS'};

for r = 2:5
    for alg_idx = 1:length(algorithms)
        algorithm = algorithms{alg_idx};
        U = results.(algorithm)(r).U;
        
        figure;
        [Consistency, G, stdG, Target] = corcond(T, U, [], 1);
        title([algorithm, ' R = ', num2str(r), ', Core Consistency: ', num2str(Consistency)]);
        grid on;
    end
end


%% functions

function U = cp_als(T, U1_0, U2_0, U3_0, maxIter, tol)
   
    U1 = U1_0;
    U2 = U2_0;
    U3 = U3_0;
    
    I = size(T, 1);
    J = size(T, 2);
    K = size(T, 3);
    R = size(U1, 2);
    
    for iter = 1:maxIter
        % Update U1
        kr_U3_U2 = kr(U3, U2);
        T1 = unfold(T, 1);
        U1_new = (T1 * kr_U3_U2) / (kr_U3_U2' * kr_U3_U2);
        
        % Update U2
        kr_U3_U1 = kr(U3, U1_new);
        T2 = unfold(T, 2);
        U2_new = (T2 * kr_U3_U1) / (kr_U3_U1' * kr_U3_U1);
        
        % Update U3
        kr_U2_U1 = kr(U2_new, U1_new);
        T3 = unfold(T, 3);
        U3_new = (T3 * kr_U2_U1) / (kr_U2_U1' * kr_U2_U1);
        
        if norm(U1_new - U1, 'fro') < tol && ...
           norm(U2_new - U2, 'fro') < tol && ...
           norm(U3_new - U3, 'fro') < tol
            break;
        end
        
        U1 = U1_new;
        U2 = U2_new;
        U3 = U3_new;
    end
    U = {U1, U2, U3};
end

function T_unfold = unfold(T, mode)
    
    sz = size(T);
    switch mode
        case 1
            T_unfold = reshape(permute(T, [1, 2, 3]), sz(1), []);
        case 2
            T_unfold = reshape(permute(T, [2, 1, 3]), sz(2), []);
        case 3
            T_unfold = reshape(permute(T, [3, 1, 2]), sz(3), []);
    end
end

function K = kr(A, B)
    
    [I, R] = size(A);
    [J, ~] = size(B);
    K = zeros(I*J, R);
    for r = 1:R
        K(:, r) = kron(A(:, r), B(:, r));
    end
end


function tmse = TMSE(U,esU)

N = length(U) ; 
P = size(U{1},2) ; 

for n=1:N
    for p=1:P
        newU{n}(:,p) = U{n}(:,p)/norm(U{n}(:,p)) ;
        new_esU{n}(:,p) = esU{n}(:,p)/norm(esU{n}(:,p)) ;
    end
end
    
permVec = perms(1:P) ;

tempdist = [] ;
for i = 1: size(permVec,1)
    for n=1:N
        A = newU{n} ;
        esA = new_esU{n} ;
     
        newA = esA(:,permVec(i,:)) ;
        diffA = [] ;
        for p = 1:P
            diffA(:,p) = A(:,p) - newA(:,p)'*A(:,p)/(newA(:,p)'*newA(:,p))*newA(:,p) ;
        end           
        tempdist(i,n) = norm(diffA,'fro')^2/norm(A,'fro')^2 ;
    end
    alldist = sum(tempdist,2) ;
end
        
tmse = min(alldist) ;

end

function T = outer_product(U1, U2, U3)
    T = zeros(length(U1), length(U2), length(U3));
    for i = 1:length(U1)
        for j = 1:length(U2)
            for k = 1:length(U3)
                T(i,j,k) = U1(i) * U2(j) * U3(k);
            end
        end
    end
end



function [s] = norm_fro(T)
    [I, J, K] = size(T);
    s = 0;
    for i = [1: 1: I]
        for j = [1: 1: J]
            for k = [1: 1: K]
                s = s + T(i, j, k) * T(i, j, k);
            end
        end
    end
end

function U = HOSVD(T)
   

    T1 = unfold(T, 1);
    T2 = unfold(T, 2);
    T3 = unfold(T, 3);

    [U1, ~, ~] = svd(T1, 'econ');
    [U2, ~, ~] = svd(T2, 'econ');
    [U3, ~, ~] = svd(T3, 'econ');
    
    U = {U1 , U2, U3};
end


