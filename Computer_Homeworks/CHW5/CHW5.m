% Nikoo Moradi
% 400101934
% CHW 5
%% 1
I = 4;
J = 3;
K = 2;
X = tensor(rand(I, J, K), [I, J, K])
R = [3 2 2];

[G, U1, U2, U3] = HOOI(X, R);

X_pred = ttm(G, {U1, U2, U3}, [1 2 3])


%% 2

clc; close all; clear;
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/tensor_toolbox-v3.2.1/@ttensor');
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/tensor_toolbox-v3.2.1/@tenmat');
addpath('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/tensor_toolbox-v3.2.1/@tensor');

% Load the ORL dataset images and create a tensor
image_folder = '/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/ORL';
num_people = 5;
num_images_per_person = 10;
image_height = 112;
image_width = 92;

% Initialize tensor
T = zeros(image_height, image_width, num_people * num_images_per_person);

% Read and store images in the tensor
for person = 1:num_people
    for img = 1:num_images_per_person
        image_index = (person - 1) * num_images_per_person + img;
        image_path = fullfile(image_folder, sprintf('s%d', person),sprintf('%d.pgm', img));
        T(:, :, image_index) = imread(image_path);
    end
end


% Perform Tucker decomposition
% R = 5;

size(T)

T = tensor(T);

colors = {'blue', 'cyan', 'green', 'red', 'magenta'};

% tucker_als
R = [50 50 5];
[G_tucker, factors] = tucker_als(T, R);


figure('Name','Using tucker-als function');
hold on;
A3 = G_tucker.U{3};
for person = 1:5
    subplot(6,1,person);
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 2;
    legend(sprintf('Person %d',person));
    grid on;
    
    subplot(6,1,6);
    hold on;
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 1.5;
    grid on;
end

% non-negative tucker
R = [50 50 5];
opts.maxit = 500;
opts.tol = 1e-5;
[G_ntd, C1, Out1] = NNTD(T, R, opts);

figure('Name','Using non-negative tucker(ntd) function, R = [50 50 5]');

A3 = G_ntd{3};
for person = 1:5
    subplot(6,1,person);
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 2;
    legend(sprintf('Class %d',person));
    grid on;
    
    subplot(6,1,6);
    hold on;
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 1.5;
    grid on;
end

% non-negative tucker
R = [5 5 5];
opts.maxit = 500;
opts.tol = 1e-5;
[G_ntd, C1, Out1] = NNTD(T, R, opts);

figure('Name','Using non-negative tucker(ntd) function, R = [5 5 5]');

A3 = G_ntd{3};
for person = 1:5
    subplot(6,1,person);
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 2;
    legend(sprintf('Class %d',person));
    grid on;
    
    subplot(6,1,6);
    hold on;
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 1.5;
    grid on;
end


% HOOI
R = [50 50 5];
[G_HOOI, U1, U2, U3] = HOOI(T,R);

figure('Name','Using my HOOI function');

A3 = U3;
for person = 1:5
    subplot(6,1,person);
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 2;
    legend(sprintf('Class %d',person));
    grid on;
    
    subplot(6,1,6);
    hold on;
    p = plot(A3(:, person), 'Color', colors{person});
    p.LineWidth = 1.5;
    grid on;
end


%% Q3 _ A


clear all; clc; close all;

% Parameters
Down = 4;
Num = [10, 25, 50, 70];
persons = 10;
images_per_person = 9;
data = zeros(persons, images_per_person, 192*168/(Down^2));


% Reading Data
for i = 1:persons
    k = 1;
    if i < 10
        base_path = ['/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/Illumination_Yale/yaleB0', num2str(i), '/yaleB0', num2str(i)];
    else
        base_path = ['/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/Illumination_Yale/yaleB', num2str(i), '/yaleB', num2str(i)];
    end
    data(i, k, :) = read_and_downsample_image([base_path, '_P00A+000E+00.pgm'], Down);
    k = k + 1;
    for j = 1:4
        data(i, k, :) = read_and_downsample_image([base_path, '_P00A-0', num2str(Num(j)), 'E+00.pgm'], Down);
        k = k + 1;
        data(i, k, :) = read_and_downsample_image([base_path, '_P00A+0', num2str(Num(j)), 'E+00.pgm'], Down);
        k = k + 1;
    end
end

% Tucker Decomposition
T1 = tucker_als(tensor(double(data)), [10 1 90]);
T3 = tucker_als(tensor(double(data)), [10 3 90]);
T5 = tucker_als(tensor(double(data)), [10 5 90]);

% Reconstruct
P_pred1 = ttm(T1.core, {T1.U{1}, T1.U{2}, T1.U{3}}, [1 2 3]);
P_pred3 = ttm(T3.core, {T3.U{1}, T3.U{2}, T3.U{3}}, [1 2 3]);
P_pred5 = ttm(T5.core, {T5.U{1}, T5.U{2}, T5.U{3}}, [1 2 3]);


% Plotting
for person = 1:4:persons
    plot_images(person, data, P_pred1, P_pred3, P_pred5, Down);
end

%% Q3 _ B 

clear all; clc; close all; 

% Parameters
Down = 4;
Num = [10, 25, 50, 70];
persons = 10;
images_per_person = 9;
data = zeros(persons, images_per_person, 192*168/(Down^2));

% Reading Data
for i = 1:persons
    k = 1;
    if i < 10
        base_path = ['/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/Illumination_Yale/yaleB0', num2str(i), '/yaleB0', num2str(i)];
    else
        base_path = ['/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/5/SimHW5/Illumination_Yale/yaleB', num2str(i), '/yaleB', num2str(i)];
    end
    data(i, k, :) = read_and_downsample_image([base_path, '_P00A+000E+00.pgm'], Down);
    k = k + 1;
    for j = 1:4
        data(i, k, :) = read_and_downsample_image([base_path, '_P00A-0', num2str(Num(j)), 'E+00.pgm'], Down);
        k = k + 1;
        data(i, k, :) = read_and_downsample_image([base_path, '_P00A+0', num2str(Num(j)), 'E+00.pgm'], Down);
        k = k + 1;
    end
end

data_mat = transpose(tenmat(data, 3));
data_mat = data_mat.data;

% Calculating SVD of denormalized
[m, n] = size(data_mat);
p = min(m, n);
[U, S, V] = svd(data_mat);

% Approximation
data_mat_approximation = zeros(m, n, p);
error = zeros(1, p);
for k = 1:p
    data_mat_approximation(:, :, k) = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
    error(k) = norm(data_mat - data_mat_approximation(:, :, k));
end

data_mat_approximation10 = data_mat_approximation(:, :, 10);
data_mat_approximation30 = data_mat_approximation(:, :, 30);
data_mat_approximation50 = data_mat_approximation(:, :, 50);

data_tensor_approximation10 = reshape(transpose(data_mat_approximation10), [10 9 2016]);
data_tensor_approximation30 = reshape(transpose(data_mat_approximation30), [10 9 2016]);
data_tensor_approximation50 = reshape(transpose(data_mat_approximation50), [10 9 2016]);


% Plotting
for person = 1:4:persons
    plot_images2(person, data, data_tensor_approximation10, data_tensor_approximation30, data_tensor_approximation50, Down);
end






%% functions 

function [G, U1, U2, U3] = HOOI(X, R)
    
    number_of_itterations = 2000;   % Number of Itteraion of outer for loop
    
    % HOSVD
    for i = [1: 1: 3]
        temp = tenmat(X, i);
        [U, ~, ~] = svd(temp.data);
        A{i} = U(:, 1:R(i));
    end
    
    % HOOI:
    for steps = [1: 1: number_of_itterations]
        A_old = A;
        
        Z = ttm(X,{A{2}', A{3}'},[2, 3]);
        temp = tenmat(Z, 1);
        [U, ~, ~] = svd(temp.data);
        A{1} = U(:, 1:R(1));

        Z = ttm(X, {A{1}', A{3}'}, [1, 3]);
        temp = tenmat(Z, 2);
        [U, ~, ~] = svd(temp.data);
        A{2} = U(:, 1:R(2));

        Z = ttm(X, {A{1}', A{2}'}, [1 2]);
        temp = tenmat(Z, 3);
        [U, ~, ~] = svd(temp.data);
        A{3} = U(:, 1:R(3));
    end
    
    U1 = A{1};
    U2 = A{2};
    U3 = A{3};
    G = ttm(X,{A{1}',A{2}',A{3}'},[1 2 3]);
    
end



function [G, U1, U2, U3] = my_HOOI(T, R1, R2, R3)
    
    % Perform HOSVD to initialize U1, U2, U3
    [U1, U2, U3] = HOSVD(T, R1, R2, R3);
    
    % Initialize G using the initialized U1, U2, U3
    G = ttm(T, {U1', U2', U3'}, [1 2 3]);
    
    % Iterate to optimize U1, U2, U3 and G
    maxIter = 100; 
    tol = 1e-6;    
    
    for iter = 1:maxIter
        
         % Update U1
        Z = ttm(G, {U2, U3}, [2 3]);
        Z1 = tenmat(Z, 1);
        T1 = tenmat(T, 1);
        U1 = nvecs(double(T1 * pinv(double(Z1))), R1, 'LR');
        
        % Update U2
        Z = ttm(G, {U1, U3}, [1 3]);
        Z2 = tenmat(Z, 2);
        T2 = tenmat(T, 2);
        U2 = nvecs(double(T2 * pinv(double(Z2))), R2, 'LR');
        
        % Update U3
        Z = ttm(G, {U1, U2}, [1 2]);
        Z3 = tenmat(Z, 3);
        T3 = tenmat(T, 3);
        U3 = nvecs(double(T3 * pinv(double(Z3))), R3, 'LR');
        
        % Update G
        G_new = ttm(T, {U1', U2', U3'}, [1 2 3]);
        
        % Check for convergence
        if norm(G_new(:) - G(:)) < tol
            break;
        end
        
        G = G_new;
    end
end

function [U1, U2, U3] = HOSVD(T, R1, R2, R3)
    % Higher-Order SVD for initialization
    [U1, ~, ~] = svds(double(tenmat(T, 1)), R1);
    [U2, ~, ~] = svds(double(tenmat(T, 2)), R2);
    [U3, ~, ~] = svds(double(tenmat(T, 3)), R3);
end

% Function to read and downsample images
function img_vec = read_and_downsample_image(filepath, Down)
    img = imread(filepath);
    img = img(1:Down:end, 1:Down:end);
    img_vec = img(:);
end

% Function to plot images
function plot_images(person_index, data, P_pred1, P_pred3, P_pred5, Down)
    figure('Name',['Person ', num2str(person_index)]);
    
    for j = 1:9
        P0 = reshape(double(data(person_index, j, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P0));
        if j == 5
            xlabel('Original');
        end
    end
    for j = 28:36
        P5 = reshape(double(P_pred5(person_index, j-27, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P5));
        if j == 32
            xlabel('Illumination = 5');
        end
    end
    for j = 19:27
        P3 = reshape(double(P_pred3(person_index, j-18, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P3));
        if j == 23
            xlabel('Illumination = 3');
        end
    end
    for j = 10:18
        P1 = reshape(double(P_pred1(person_index, j-9, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P1));
        if j == 14
            xlabel('Illumination = 1');
        end
    end
    
    
end


% Function to plot images
function plot_images2(person_index, data, data_tensor_approximation10, data_tensor_approximation30, data_tensor_approximation50, Down)
    figure;
    for j = 1:9
        P0 = reshape(double(data(person_index, j, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P0));
        if j == 5
            xlabel('Original');
        end
    end
    for j = 10:18
        P1 = reshape(double(data_tensor_approximation10(person_index, j-9, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P1));
        if j == 14
            xlabel('Rank = 10');
        end
    end
    for j = 19:27
        P3 = reshape(double(data_tensor_approximation30(person_index, j-18, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P3));
        if j == 23
            xlabel('Rank = 30');
        end
    end
    for j = 28:36
        P5 = reshape(double(data_tensor_approximation50(person_index, j-27, :)), 192/Down, 168/Down);
        subplot(4, 9, j);
        imshow(uint8(P5));
        if j == 32
            xlabel('Rank = 50');
        end
    end
    annotation('textbox', [0.5, 0.98, 0, 0], 'string', ['Person ', num2str(person_index)], 'FitBoxToText', 'on', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'LineStyle', 'none');
end

