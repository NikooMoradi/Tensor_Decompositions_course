% Nikoo Moradi _ 400101934 
% CHW 1
%% Question 1
clc; clear;

A1 = [1 sqrt(2) 2; sqrt(2) 3 sqrt(2); 2 sqrt(2) 1];
A2 = [1 1 1; 1 1 1; 1 1 1];
A3 = [7 2 1; 2 5 1; 1 1 9];

A1
[V1, D1] = Jacobi_eig(A1)
[V11, D11] = eig(A1)

A2
[V2, D2] = Jacobi_eig(A2)
[V22, D22] = eig(A2)

A3
[V3, D3] = Jacobi_eig(A3)
[V33, D33] = eig(A3)

%% Question 2
clc; clear;

A1 = [1 2 3; 1 5 3; 1 0 0; 0 2 2]
[U1,S1,V1] = Jacobi_svd_2sided(A1)
[U11,S11,V11] = svd(A1)

A2 = [1 0 3; 1 5 3; 1 0 7; 1 sqrt(5) 2]
[U2,S2,V2] = Jacobi_svd_2sided(A2)
[U22,S22,V22] = svd(A2)

A3 = [0 9 2; 7 1 4; 5 8 sqrt(7); 3 sqrt(6) 1]
[U3,S3,V3] = Jacobi_svd_2sided(A3)
[U33,S33,V33] = svd(A3)

%% Question 3
clc; clear;

A1 = [1 2 ; 1 5 ; 1 0 ; 0 2 ]
[U1,S1,V1] = Jacobi_svd_1sided(A1)
[U11,S11,V11] = svd(A1)

A2=[1,0 ; 0,1; 1,1]
[U2,S2,V2] = Jacobi_svd_1sided(A2)
[U22,S22,V22] = svd(A2)

% A3 = [1 5 2; 3 2 4; 5 6 sqrt(2); 2 sqrt(3) 7]
A3 = [1  8; 2  1; 5  sqrt(2); 2 6]
[U3,S3,V3] = Jacobi_svd_1sided(A3)
[U33,S33,V33] = svd(A3)
%% Functions 3

function [c,s] = making_orthogonal(x,y)
    a = norm(x)^2;
    b = norm(y)^2;
    c = dot(x,y);
    
    if c ~= 0
        tau = (b - a)/(2*c);
        if tau >= 0
            t = 1/(tau + sqrt(1 + tau^2));
        end
        if tau < 0
            t = 1/(tau + sqrt(1 + tau^2));
        end
        
    c = 1/sqrt(1 + t^2);
    s = c * t;
    
    else 
        c = 1;
        s = 0;
    end
end

function [U,S,V]=Jacobi_svd_1sided(A)
    [m,n] = size(A);
    V = eye(n,n);
    S = zeros(m,n);
    U = zeros(m,m);
    e = 10^(-30);
    err = 0;
    not_orthogonal = true;
    
    while not_orthogonal
        for p = 1:n-1
            for q = p+1 : n
                
                x = A(:,p);
                y = A(:,q);
                [cos,sin] = making_orthogonal(x,y);
                
                J = eye(n,n);   
                J([p,q],[p,q]) = [cos, sin; -sin, cos];
                
                A = A * J;
                V = J' * V;
                
                
            end
        end
        err = 0;
        for p = 1:n
            for q = p+1 : n
                err = err + (dot(A(:,p),A(:,q)));
            end
        end
        
        if abs(err) < e
            not_orthogonal = false;
        else
            not_orthogonal = true;
        end
        
    end 
    
    for j = 1:min(m,n) 
        S(j,j) = norm(A(:,j)); 
        U(:,j) = A(:,j) / norm(A(:,j));
    end
    
end
%% Functions 1

function [V,D] = Jacobi_eig(A)
    
    % Initiating D and V
    D = A;
    V = eye(size(A));
    
    e = 10^(-10);
    
    
    while off(A) > e
        
        [row, col] = size(A);
        
        % Finding the maximum off-diagonal element of a (p,q)
        p = 0;
        q = 0;
        max_val = 0;
        for i = 1:row
            for j = i+1:col  
                current_val = abs(A(i,j)); 
                if current_val > max_val
                    max_val = current_val;
                    p = i;
                    q = j;
                end
            end
        end
        
        % Calculating theta and sin and cos 
        apq = D(p,q);
        app = D(p,p);
        aqq = D(q,q);
        
        theta = atan((2*apq)/(app - aqq))/2;
        
        c = cos(theta);
        s = sin(theta);
        
        % Forming J (Givens matrix)
        J = eye(size(A));   
        J([p,q],[p,q]) = [c, -s; s, c];

        
        % Calculating new D and V
        D = (J')*D*J;
        V = V*J;
        A = D;
        

    end
 

end

%% Functions 2


function [c,s] = symSchur2(A,p,q)
    if A(p,q) ~= 0
        tau = (A(q,q)-A(p,p))/(2*A(p,q));
        if tau >= 0
            t = 1/(tau + sqrt(1 + tau^2));
        else
            t = 1/(tau - sqrt(1 + tau^2));
        end
        c = 1/sqrt(1 + t^2);
        s = c * t;
    
    else 
        c = 1;
        s = 0;
    end

end

function [c,s] = making_symmetric(A,p,q)
        app = A(p,p);
        apq = A(p,q);
        aqp = A(q,p);
        aqq = A(q,q);
                
        t = (aqp - apq)/(app + aqq);
        c = 1/sqrt(1 + t^2);
        s = t*c;    
end

function [U,S,V] = Jacobi_svd_2sided(A)
    [m,n] = size(A);
%     S = A;
    V = eye(n,n);
    U = eye(m,m);
    e = 10^(-20);
    
    while off(A) > e
        for p = 1:n
            for q = p+1 : n
        
                [cos,sin] = making_symmetric(A,p,q);
                
                symmetric = eye(m,m);    
                symmetric([p,q],[p,q]) = [cos, sin; -sin, cos];
                
                A = symmetric * A;
                
                [c1,s1] = symSchur2(A,p,q);
                J1 = eye(m,m);   
                J1([p,q],[p,q]) = [c1, s1; -s1, c1];
                J2 = eye(n,n);   
                J2([p,q],[p,q]) = [c1, s1; -s1, c1];
                
                A = (J1') * A * J2;
                U = U * symmetric' * J1;
                V = V * J2;
                
            end
            for q = n+1 : m
                
                [c3,s3] = setting_to_zero(A,p,q);
                
                J3 = eye(m,m);   
                J3([p,q],[p,q]) = [c3, -s3; s3, c3];
                A = J3' * A;
                U = U * J3;
            end
             
        end
        
        

    end
    S = A;
end

function [c,s] = setting_to_zero(A,p,q)
          app = A(p,p);
          aqp = A(q,p);
          t3 = (aqp / app);
          c3 = 1/sqrt(1 + (t3)^2);
          s3 = t3*c3; 
end

function Off = off(A)

    [row, col] = size(A);
     Off = 0;
    for i = 1:min(row,col)
        for j = 1:min(row, col)
            if i ~= j
                Off = Off + A(i,j)^2;
            end
        end
    end
    Off = sqrt(Off);
    
end