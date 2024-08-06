% Nikoo Moradi
% 400101934
% CHW 2 
% Q1
%% A
clc; clear; close all;
A = randseed(400,6, 4) % Example non-negative matrix A
j = 3; % Rank of the factorization
B0 = randseed(400,6, j); % Initial B matrix
C0 = randseed(400,j, 4); % Initial C matrix


maxIter = 100; % Maximum number of iterations
tol = 1e-3; % Tolerance for convergence

[B, C, E] = als_algorithm(A, B0, C0, j, maxIter, tol);
A_pred_als = B*C 
error_my_als = norm(E, 'fro')

[B,C] = nnmf(A,j,'algorithm','als','w0',B0,'h0',C0);
A_pred_matlab_als = B*C
error_nnmf_als = norm(A - B * C, 'fro')

%% B
clc; clear; close all;

A = randseed(400,6, 4) % Example non-negative matrix A
j = 3; % Rank of the factorization
B0 = randseed(400,6, j); % Initial B matrix
C0 = randseed(400,j, 4); % Initial C matrix
E0 = randseed(400,j, 4);

maxIter = 100; % Maximum number of iterations
tol = 1e-4; % Tolerance for convergence

[B, C, E] = multiplicative_algorithm(A, B0, C0, j, maxIter, tol);

A_pred_mult = B * C 
error_my_mult = norm(E, 'fro')
[B,C] = nnmf(A,j,'algorithm','mult','w0',B0,'h0',C0);
A_pred_matlab_mult =B*C
error_nnmf_mult = norm(A - B * C, 'fro')

%% C 
%% 1
clc; clear; close all;

% A = randseed(400,6, 4) % Example non-negative matrix A

j = 3; % Rank of the factorization
B0 = rand(6, j); % Initial B matrix
C0 = rand(j, 4); % Initial C matrix
E0 = rand(6, 4);

%% 2


SNR = [-10,0,10,30,50];
a = zeros(1,5);

A = zeros(6,4,5);

for i = 1 : size(SNR,2)
    W = norm(B0 * C0, 'fro');
    H = norm(E0, 'fro');
    help = 10 ^ (SNR(i) / 20);
    a(i) = W / (help * H);
    A(:,:,i) = B0 * C0 + a(i) * E0;
end



%% 3

E_my_als =  zeros(1,5);
E_my_mult =  zeros(1,5);
E_nnmf_als =  zeros(1,5);
E_nnmf_mult =  zeros(1,5);

j = 3;
B0 = rand(6, j); % Initial B matrix
C0 = rand(j, 4); % Initial C matrix

maxIter = 10000; % Maximum number of iterations
tol = 1e-4; % Tolerance for convergence

for i = 1 : 5
    
            
    [B, C, E_temp1] = als_algorithm(A(:,:,i), B0, C0, j, maxIter, tol);
    E_my_als(i) = norm(E_temp1,'fro');
    
    [B, C, E_temp2] = multiplicative_algorithm(A(:,:,i), B0, C0, j, maxIter, tol);
    E_my_mult(i) = norm(E_temp2,'fro');
    
    [B,C] = nnmf(A(:,:,i),j,'algorithm','als','w0',B0,'h0',C0);
    E_nnmf_als(i) = norm(A(:,:,i) - B * C, 'fro');
    
    [B,C] = nnmf(A(:,:,i),j,'algorithm','mult','w0',B0,'h0',C0);
    E_nnmf_mult(i) = norm(A(:,:,i) - B * C,'fro');
end

E_my_als
E_my_mult
E_nnmf_als
E_nnmf_mult

%% 4
clc; clear; close all;

SNR = [-10,0,10,30,50];
maxIter = 100; % Maximum number of iterations
tol = 1e-3; % Tolerance for convergence

fig = figure;
% set(fig, 'Position', [100, 100, 800, 800]);

E_my_als_per_j =  zeros(3,5);
E_my_mult_per_j =  zeros(3,5);

for j = 2 : 4
    
    j
    A = create_A(SNR,j);
    
    E_my_als_avg =  zeros(1,5);
    E_my_mult_avg =  zeros(1,5);
    E_nnmf_als_avg =  zeros(1,5);
    E_nnmf_mult_avg =  zeros(1,5);

    for i = 1 : 5
        for n = 1 : 10
       
            B0 = rand(6, j); % Initial B matrix
            C0 = rand(j, 4); % Initial C matrix
        
            [B, C, E_temp1] = als_algorithm(A(:,:,i), B0, C0, j, maxIter, tol);
            E_my_als_avg(i) =  E_my_als_avg(i) + norm(E_temp1,'fro');
         
            [B, C, E_temp2] = multiplicative_algorithm(A(:,:,i), B0, C0, j, maxIter, tol);
            E_my_mult_avg(i) =  E_my_als_avg(i) + norm(E_temp2 ,'fro');
        
            [B,C] = nnmf(A(:,:,i),j,'algorithm','als','w0',B0,'h0',C0);
            E_nnmf_als_avg(i) = E_nnmf_als_avg(i) + norm(A(:,:,i) - B * C,'fro');
        
            [B,C] = nnmf(A(:,:,i),j,'algorithm','mult','w0',B0,'h0',C0);
            E_nnmf_mult_avg(i) = E_nnmf_mult_avg(i) + norm(A(:,:,i) - B * C,'fro');
        end
        E_my_als_avg(i) =  E_my_als_avg(i) / 10;
        E_my_mult_avg(i) =  E_my_als_avg(i) / 10;
        E_nnmf_als_avg(i) = E_nnmf_als_avg(i) / 10;
        E_nnmf_mult_avg(i) = E_nnmf_mult_avg(i) / 10;
    end
    subplot(1,3,j-1);
    hold on

    p(1)=plot(SNR,E_my_als_avg,'*-');
    p(1).LineWidth = 2;
    
    p(2)=plot(SNR,E_my_mult_avg,'s-');
    p(2).LineWidth = 2;
    
    p(3)=plot(SNR,E_nnmf_als_avg,'+-');
    p(3).LineWidth = 2;
    
    p(4)=plot(SNR,E_nnmf_mult_avg,'o-');
    p(4).LineWidth = 2;
    
    grid on;
    xlabel('SNR');
    ylabel('Errors');
    title(sprintf('j =  %d', j));
    legend('my als','my Mult','nnmf als','nnmf Mult');
    
    E_my_als_per_j(j-1,:) =  E_my_als_avg;
    E_my_mult_per_j(j-1,:) =  E_my_mult_avg;
    
end

figure;
subplot(1,2,1);
hold on
p(1)=plot(SNR,E_my_als_per_j(1,:),'*-');
p(1).LineWidth = 2;
    
p(2)=plot(SNR,E_my_als_per_j(2,:),'s-');
p(2).LineWidth = 2;
    
p(3)=plot(SNR,E_my_als_per_j(3,:),'+-');
p(3).LineWidth = 2;
grid on;
xlabel('SNR');
ylabel('Errors');
title('My ALS for differnet j');
legend('j = 2','j = 3','j = 4');


subplot(1,2,2);
hold on
p(1)=plot(SNR,E_my_mult_per_j(1,:),'*-');
p(1).LineWidth = 2;
    
p(2)=plot(SNR,E_my_mult_per_j(2,:),'s-');
p(2).LineWidth = 2;
    
p(3)=plot(SNR,E_my_mult_per_j(3,:),'+-');
p(3).LineWidth = 2;
grid on;
xlabel('SNR');
ylabel('Errors');
title('My Mult for differnet j');
legend('j = 2','j = 3','j = 4');



%% Functions

function [B, C, E] = als_algorithm(A, B0, C0, j, maxIter, tol)
   
    % Initialize B and C
    B = B0;
    C = C0;
    eps = 1e-16;
    
        

        for iter = 1: maxIter
            C_prime = pinv(B' * B) * (B' * A);
            C = max(eps,C_prime);
        
            B_prime = (A * C') * pinv(C * C');
            B = max(eps,B_prime);
        
            % Calculate error
            E = A - B * C;
            errorNorm = norm(E, 'fro');
        
            % Check for convergence
            if errorNorm < tol
                break;
            end
        
        end


end



function [B, C, E] = multiplicative_algorithm(A, B0, C0, j, maxIter, tol)

    % Initialize B and C
    B = B0;
    C = C0;
    n = size(A, 1);
    m = size(A, 2);
    
    for iter = 1:maxIter
        % Update C
        C = C .* (B' * A) ./ (B' * B * C + eps);
        
        % Update B
        B = B .* (A * C') ./ (B * (C * C') + eps);
        
        % Calculate error
        E = A - B * C;
        errorNorm = norm(E, 'fro');
        
        % Check for convergence
        if errorNorm < tol
            break;
        end
    end
end


function A = create_A(SNR,j)
    B0 = rand(6, j); % Initial B matrix
    C0 = rand(j, 4); % Initial C matrix
    E0 = rand(6, 4);
    a = zeros(1,5);

    A = zeros(6,4,5);

    for i = 1 : size(SNR,2)
        W = norm(B0 * C0, 'fro');
        H = norm(E0, 'fro');
        help = 10 ^ (SNR(i) / 20);
        a(i) = W / (help * H);
        A(:,:,i) = B0 * C0 + a(i) * E0;
    end
    
end
