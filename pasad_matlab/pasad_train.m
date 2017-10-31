function U = pasad_train(series,N,L)

close all;
K = N-L+1;
X = hankel(series(1:L),series(L:length(series)));
Xb = X(:,1:K);

disp('SVD Decomposition started...');tic
[t,e,~] = svd(Xb);
d = diag(e);
disp('SVD Decomposition complete');toc

G = (d(2:end)./sum(d(2:end)))*100;
figure
plot(G,'color',[.4 .4 .4],'linewidth',2),hold on,plot(G,'rx','color',[1 .4 .2]);
xlabel('Number of eigenvalues');
ylabel('Eigenvalue share')
title('Screen plot');
set(gca,'fontsize',16);
r = input('Choose the statistical dimension ');
close all
I = (1:r);
U = t(:,I);

% Write U in a file
fileID = fopen('U.txt','w');
for ii= 1:size(U,1)
    fprintf(fileID,'%g\t',U(ii,:));
    fprintf(fileID,"\n");
end
fclose(fileID);

% Write s in a file
fileID = fopen('s.txt','w');
for ii= 1:size(series,1)
    fprintf(fileID,'%g\t',series(ii,:));
    fprintf(fileID,"\n");
end
fclose(fileID);

disp('Training PASAD is complete. The output is the projection matrix U.');
