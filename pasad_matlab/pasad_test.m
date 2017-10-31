function D = pasad_test(series,N,L,U)

X = hankel(series(1:L),series(L:length(series)));
T = size(X,2);
K = N-L+1;
D = zeros(T-K-1,1);
tic;
    for n = N-L+2:T
        Z = X(:,n);
        D(n-K) = (Z'*Z - (Z'*U)*(U'*Z))/L;
    end

disp('Testing complete.');
toc;
end