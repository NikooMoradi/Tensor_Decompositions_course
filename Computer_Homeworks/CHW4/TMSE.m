function tmse = TMSE(U,esU)

N = length(U) ; %N-way array
P = size(U{1},2) ; %Rank

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