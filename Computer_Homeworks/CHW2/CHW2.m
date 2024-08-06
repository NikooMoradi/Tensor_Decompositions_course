% Nikoo Moradi
% 400101934
% CHW 2 
% Q1
%%


%% using normalization

clc; clear; close all;
image = imread('/Users/nikoomoradi/Desktop/SlUT/6/Tensor/CHW/2/SimHW2/Q1/cameraman.jpg');
image = rgb2gray(image);

image = double(image);


% min = min(image,[],"all");
temp = min(image,[],2);
min_val = min(temp, [], 1);
max_val = max(max(image,[],2),[],1);

normalized_image = (image - min_val) / (max_val - min_val);

[U,S,V] = svd(normalized_image);

ranks = [1,5,15,35,50,70,100, 200];
errors_with = zeros(size(ranks));
compress = zeros(size(ranks));
temp = (256*2 + 1) / (256 * 256);

figure('Name','With normalization');
subplot(3, 3, 1);
imshow(normalized_image);
title('Original Image');

% Iterate over each rank
for i = 1:length(ranks)
    k = ranks(i);
    
    % Reconstruct the image with k rank
    reconstructed_image = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
    
    % Calculate error (root mean squared error)
    errors_with(i) = sqrt(mean((normalized_image(:) - reconstructed_image(:)).^2));
    compress(i) = temp * k
    
    % Display the reconstructed image
    subplot(3, 3, i+1);
    imshow(uint8(255*reconstructed_image));
    title(['Rank = ' num2str(k)]);
end
figure;
% Display the error plot vs rank
p(1)=plot(ranks, errors_with,'*-');
p(1).LineWidth = 2;
xlabel('Rank');
ylabel('Root Mean Squared Error');
title('Reconstruction Error - with normalization')
legend('error: MSE');

figure;
% Display the error plot vs compress
p(1)=plot(compress,errors_with, 'o-');
p(1).LineWidth = 2;
xlabel('Amount of compression');
ylabel('Root Mean Squared Error');
legend('error2: MSE');
title('Reconstruction Error - with normalization');

%% Not using normalization
% close all;

image = imread('cameraman.jpg');
image = rgb2gray(image);
image = double(image);


[U,S,V] = svd(image);

ranks = [1,5,15,35,50,70,100, 200];

% errors = zeros(size(ranks));
errors_without = zeros(size(ranks));
compress = zeros(size(ranks));
temp = (256*2 + 1) / (256 * 256);

figure('Name','Without normalization');
subplot(3,3 , 1);
imshow(uint8(image));
title('Original Image');

% Iterate over each rank
for i = 1:length(ranks)
    k = ranks(i);
    
    % Reconstruct the image with k rank
    reconstructed_image = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
    
    % Calculate error (root mean squared error)
%     errors(i) = norm(image - reconstructed_image, 'fro') / norm(image, 'fro');
    errors_without(i) = sqrt(mean((image(:) - reconstructed_image(:)).^2));
    compress(i) = temp * k;
    
    % Display the reconstructed image
    subplot(3, 3, i+1);
    imshow(uint8(reconstructed_image));
    title(['Rank = ' num2str(k)]);
end

figure;
% Display the error plot vs rank
p(1)=plot(ranks, errors_without, 'o-');
p(1).LineWidth = 2;
xlabel('Rank');
ylabel('Root Mean Squared Error');
title('Reconstruction Error - without normalization')
legend('error2: MSE');

figure;
% Display the error plot vs compress
p(1)=plot(compress,errors_without,'o-');
p(1).LineWidth = 2;
xlabel('Amount of compression');
ylabel('Root Mean Squared Error');
legend('error2: MSE');
title('Reconstruction Error - without normalization');

%%

figure;
hold on;
p(1)=plot(ranks, errors_without, 'ro-');
p(2)=plot(ranks, errors_with, 'b*-');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
xlabel('Rank');
ylabel('Root Mean Squared Error');
legend('error without normalization','error with normalization');
title('Reconstruction Error (MSE)');


figure;
hold on;
p(3)=plot(compress, errors_without, 'g.-');
p(4)=plot(compress, errors_with, 'c--');
p(3).LineWidth = 2;
p(4).LineWidth = 2;
legend('error without normalization','error with normalization');
xlabel('Amount of compression');
ylabel('Root Mean Squared Error');
title('Reconstruction Error (MSE)');

