% Nikoo Moradi
% 400101934
% CHW 2 
% Q2
%% 
clc; clear; close all;

EEGdata = load('EEGdata.mat');
Xorg = EEGdata.Xorg;
Xnoise = EEGdata.Xnoise;

[U,S,V] = svd(Xnoise);
disp(size(S));
ranks = 1:1:32;
errors = zeros(size(ranks));

for i = 1:length(ranks)
    k = ranks(i);
    
    % Reconstruct the image with k rank
    reconstructed_X = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
    
    % Calculate error (root mean squared error)
    errors(i) = norm(Xorg - reconstructed_X, 'fro') ;
       
end

figure;
% Display the error plot vs rank
p(1)=plot(ranks, errors,'*-');
p(1).LineWidth = 2;
grid on;
xlabel('Rank');
ylabel('Frobinious Error');
title('Reconstruction Error ')
legend('error');

[min_error, min_source] = min(errors);

% Display the minimum value and its index
disp(['Minimum error: ', num2str(min_error)]);
disp(['Optimal number of non-noise sources: ', num2str(min_source)]);

