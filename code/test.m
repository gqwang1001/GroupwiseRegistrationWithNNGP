testMat = ones(100, 100);
[coord_x, coord_y] = find(ones(size(testMat)));
covs = exp(-squareform(pdist([coord_x, coord_y])));

tic
R = chol(covs);
toc
% covs_recover = R' * R;
Rsp = R; Rsp(abs(Rsp)<1e-5)=0;
Rsp = sparse(Rsp);
covsSp = Rsp'*Rsp;

r = symrcm(covsSp);
figure; spy(covsSp);
figure; spy(covsSp(r,r));

tic
U = chol(covsSp(r, r));
toc

tic
U = chol(covsSp);
toc


spI = eye(length(coord_x));
spI = ones(length(coord_x),1);
tic
Invs1 = R\spI;
toc

opts.UT = true;
tic
Invs1SP = linsolve(R, spI, opts);
toc


tic
Iv1 = Invs1* Invs1';
toc

tic
Iv1 = Invs1SP(1:1600,:) * Invs1SP(1:1600,:)';
toc

tic
Iv2 = inv(covs);
toc

%%
close all

nN = 20;
X = magic(10);
[coord_x, coord_y] = find(ones(size(X)));
covs = exp(-squareform(pdist([coord_x, coord_y])));

coord = [coord_x, coord_y];
Ncoord = size(coord,1);
coord1 = coord/1.5+1/5;
coord2 = coord/3;
S1 = exp(-squareform(pdist(coord1)));

K1 = exp(-pdist2(coord1, coord));
K2 = exp(-pdist2(coord2, coord));

InvCov = inv(covs);
TrueCov = K1*InvCov*K2';

nnx = nnsearch_2d(coord, nN);
nnDict = nnDictionary_2d(coord, 10, nN);

nnIDXs1 = nns_2d_square(nnDict, coord, coord1, 1, nN);
nnIDXs2 = nns_2d_square(nnDict, coord, coord2, 1, nN);

condX = CondLatentX_2d(X(:), coord, 1, nnx, true,0);
covInvNNGP = (eye(Ncoord)-condX.A)'*condX.Dinv*(eye(Ncoord)-condX.A);
% figure;
% subplot(1,2,1);imagesc(covInvNNGP); colormap jet; colorbar;
% subplot(1,2,2); imagesc(InvCov); colormap jet; colorbar;
covnngp = inv(covInvNNGP);
% 
% % case 1 Cov_(20,20)
% 
% knn1 = exp(-nnx.Dists);
% c1 = exp(-nnIDXs1.Dists);
% c2 = exp(-nnIDXs2.Dists);
% 
% cov1 = [];
% for i1=1:Ncoord
%     for i2=1:Ncoord
%         cov1(i1, i2) = c1(i1,:)*InvCov(nnIDXs1.Idxs(i1,:),nnIDXs2.Idxs(i2,:)) * c2(i2, :)';
%     end
% end
% ComparisonFig(cov1, TrueCov)
% case 2
b1cov = Cov_NN_Transfer_2d(X(:),coord, nnIDXs1, 1, nN, 1e-10);
b2cov = Cov_NN_Transfer_2d(X(:),coord, nnIDXs2, 1, nN, 1e-10);

b1 = b1cov.Bt;
b2 = b2cov.Bt;

% cov2 = [];
% for i1=1:Ncoord
%     for i2=1:Ncoord
%         cov2(i1, i2) = b1(i1,:)*covnngp(nnIDXs1.Idxs(i1,:),nnIDXs2.Idxs(i2,:))*b2(i2, :)';
%     end
% end
% ComparisonFig(cov2, TrueCov)

cov3 = [];
for i1=1:Ncoord
    for i2=i1:Ncoord
        cov3(i1, i2) = b1(i1,:)*covs(nnIDXs1.Idxs(i1,:),nnIDXs1.Idxs(i2,:))*b1(i2,:)';
        cov3(i2, i1) = cov3(i1,i2);
    end
    cov3(i1,i1) = cov3(i1,i1)+b1cov.Ft(i1);
end

ComparisonFig(cov3,S1);
ComparisonFig(S1-cov3,S1);
ComparisonFig(inv(cov3),inv(S1));
ComparisonFig(inv(cov3)-inv(S1),inv(S1));




%% full marginal approximation
nsubj=2;
Xmat = magic(10);
X = Xmat(:);
testMat = ones(10, 10);
Ncoord = length(testMat(:));
[coord_x, coord_y] = find(ones(size(testMat)));
covsX = exp(-squareform(pdist([coord_x, coord_y])));

InvCov = inv(covsX);
coord = [coord_x, coord_y];

coord1 = coord/2;
coord2 = coord/3;

K1 = exp(-pdist2(coord1, coord));
K2 = exp(-pdist2(coord2, coord));
K12 = K1*InvCov*K2';

S1_1 = K1*InvCov*K1';
figure; imagesc(S1_1); colormap jet; colorbar;
figure; imagesc(S1); colormap jet; colorbar;

figure; imagesc(inv(S1_1)); colormap jet; colorbar;

S1 = exp(-squareform(pdist(coord1)));
S2 = exp(-squareform(pdist(coord2)));

S12 = [S1, K12; K12', S2]+1e-9*eye(nsubj*Ncoord);
invS12 = inv(S12);
figure; imagesc(S12); colormap jet; colorbar;
figure; imagesc(inv(S12)); colormap jet; colorbar;caxis([-7,7]);
figure; imagesc(abs(invS12)>1e-2); colormap jet; colorbar;

nN = 10;
nnx = nnsearch_2d(coord, nN);
condX = CondLatentX_2d(X(:), coord, 1, nnx, true);
covInvNNGP = (eye(Ncoord)-condX.A)'*condX.Dinv*(eye(Ncoord)-condX.A);
figure;
subplot(1,2,1);imagesc(covInvNNGP); colormap jet; colorbar;
subplot(1,2,2); imagesc(InvCov); colormap jet; colorbar;

nnDict = nnDictionary_2d(coord, 10, nN);
nnIDXs1 = nns_2d_square(nnDict, coord, coord1, 1, nN);
nnIDXs2 = nns_2d_square(nnDict, coord, coord2, 1, nN);
nn1 = nnsearch_2d(coord1, nN);
nn2 = nnsearch_2d(coord2, nN);

Cov_NN1 = Cov_NN_Transfer_2d(X, coord, nnIDXs1, 1, nN, 1e-10);
Cov_NN2 = Cov_NN_Transfer_2d(X, coord, nnIDXs2, 1, nN, 1e-10);

[idx,D] = knnsearch(coord2, coord1, 'K', nN);
nnIDXs21.Idxs = idx;
nnIDXs21.Dists = D;

nngp_cov1 = nngpCholesky_marginal(nN, Ncoord, nn1, covsX, nnIDXs1, Cov_NN1, true);
cov1_InvNNGP = (eye(Ncoord)-nngp_cov1.A)'*nngp_cov1.Dinv*(eye(Ncoord)-nngp_cov1.A);

figure;
subplot(1,2,1);imagesc(cov1_InvNNGP); colormap jet; colorbar;
invS1 = inv(S1);
subplot(1,2,2); imagesc(invS1); colormap jet; colorbar;

figure;
subplot(1,2,1);imagesc(inv(cov1_InvNNGP)); colormap jet; colorbar;
subplot(1,2,2); imagesc(S1); colormap jet; colorbar;

nngp_cov1_1 = CondLatentX_2d(X(:), coord1, 1, nn1, true);
cov_1_1InvNNGP = (eye(Ncoord)-nngp_cov1_1.A)'*nngp_cov1_1.Dinv*(eye(Ncoord)-nngp_cov1_1.A);
figure;
subplot(1,2,1);imagesc(cov_1_1InvNNGP); colormap jet; colorbar;
invS1 = inv(S1);
subplot(1,2,2); imagesc(invS1); colormap jet; colorbar;

k32 = Cov_NN1.Bt(3,:)*covsX(nnIDXs1.Idxs(3,:),nnIDXs1.Idxs(2,:))*Cov_NN1.Bt(2,:)';
%% predictive inference
close all
nN = 20;
tau = 1e-3;
nsubj = 2;
nx = 10;
X = magic(nx);
testMat = ones(nx, nx);
Ncoord = length(X(:));
[coord_x, coord_y] = find(ones(size(testMat)));
covs = exp(-squareform(pdist([coord_x, coord_y])));
InvCov = inv(covs);

coord = [coord_x, coord_y];

nnx = nnsearch_2d(coord, nN);
nnDict = nnDictionary_2d(coord, 10, nN);

condX = CondLatentX_2d(X(:), coord, 1, nnx, true,0);
covInvNNGP = (eye(Ncoord)-condX.A)'*condX.Dinv*(eye(Ncoord)-condX.A);
covnngp = inv(covInvNNGP);
coord1 = coord/1.5+1/5;
coord2 = coord/3;
nnDict = nnDictionary_2d(coord, 10, nN);
nnIDXs1 = nns_2d_square(nnDict, coord, coord1, 1, nN);
nnIDXs2 = nns_2d_square(nnDict, coord, coord2, 1, nN);

K1 = exp(-pdist2(coord1, coord));
K2 = exp(-pdist2(coord2, coord));
K12 = K1*InvCov*K2';
S1 = exp(-squareform(pdist(coord1)))+tau*eye(Ncoord);
S2 = exp(-squareform(pdist(coord2)))+tau*eye(Ncoord);
invS1 = inv(S1);

b1cov = Cov_NN_Transfer_2d(X(:),coord, nnIDXs1, 1, nN, 0);
b2cov = Cov_NN_Transfer_2d(X(:),coord, nnIDXs2, 1, nN, 0);

b1 = b1cov.Bt;
b2 = b2cov.Bt;

[idx21,Dist21] = knnsearch(coord1, coord2, 'K', nN);

cov21nngp = zeros(Ncoord);
u21SinvNNGP = zeros(Ncoord);
cov2given1nngp = zeros(Ncoord);
for i=1:Ncoord
    for j=1:nN
        cov21nngp(i, idx21(i,j))=b2(i,:)*covnngp(nnIDXs2.Idxs(i,:),nnIDXs1.Idxs(idx21(i,j),:))*b1(idx21(i,j), :)';
    end
    u21SinvNNGP(i,idx21(i,:))=cov21nngp(i,idx21(i,:))/(S1(idx21(i,:),idx21(i,:)));
end

for i=1:Ncoord
    for j=1:Ncoord
        cov2given1nngp(i,j)=S2(i,j)-u21SinvNNGP(i,idx21(i,:))*S1(idx21(i,:),idx21(j,:))*u21SinvNNGP(j,idx21(j,:))';
    end
end

cov12nngp = cov21nngp';
% figure; 
% subplot(1,2,1); imagesc(cov12nngp); colormap jet; colorbar;title("NNGP");
% subplot(1,2,2); imagesc(K12); colormap jet; colorbar;title("Full GP");

nn1 = nnsearch_2d(coord1, nN);
nngp_cov1_1 = CondLatentX_2d(X(:), coord1, 1, nn1, true, tau);
A1 = nngp_cov1_1.A;
D1inv = nngp_cov1_1.Dinv;

% cov2given1 = S2-K12'*invS1*K12;
% 
% figure;
% subplot(1,2,1);imagesc(cov2given1nngp); colormap jet; colorbar; title("NNGP");
% subplot(1,2,2); imagesc(cov2given1); colormap jet; colorbar;title("Full GP");

% figure;
% subplot(1,2,1);imagesc(inv(cov2given1nngp)); colormap jet; colorbar; title("NNGP");
% subplot(1,2,2); imagesc(inv(cov2given1)); colormap jet; colorbar;title("Full GP");

nn2 = nnsearch_2d(coord2, nN);
A2 = zeros(Ncoord);
D2inv = eye(Ncoord);
D2inv(1,1) = 1/cov2given1nngp(1,1);
for i=2:Ncoord
    nnidx = nn2.Idxs(i,:);
    nidx = nnidx(nnidx~=0);
    A2(i, nidx) = cov2given1nngp(i,nidx)/(cov2given1nngp(nidx,nidx));
    D2inv(i,i) = 1/(cov2given1nngp(i,i)-cov2given1nngp(i,nidx)/(cov2given1nngp(nidx,nidx))*cov2given1nngp(i,nidx)');
end
cov2given1InvNNGP = (eye(Ncoord)-A2)'*D2inv*(eye(Ncoord)-A2);
% cov2given1Inv = inv(cov2given1);

% figure;
% subplot(1,2,1);imagesc(cov2given1InvNNGP); colormap jet; colorbar;title("NNGP");
% subplot(1,2,2);imagesc(inv(cov2given1)); colormap jet; colorbar;title("Full GP");

A21 = (eye(Ncoord)-A2)*u21SinvNNGP;
A2subjs = [A1, zeros(Ncoord); A21, A2];
D2subjsInv = diag([diag(D1inv);diag(D2inv)]);

cov2subjsNNGPInv = (eye(nsubj*Ncoord)-A2subjs)'*D2subjsInv*(eye(nsubj*Ncoord)-A2subjs);
% cov2subjsNNGP = inv(cov2subjsNNGPInv);
% 
S12 = [S1, K12; K12', S2];
invS12 = inv(S12);

% figure;
% subplot(1,2,1);imagesc(cov2subjsNNGPInv); colormap jet; colorbar;title("NNGP");
% subplot(1,2,2);imagesc(invS12); colormap jet; colorbar;title("Full GP");
% % 
% figure;
% subplot(1,2,1);imagesc(cov2subjsNNGP); colormap jet; colorbar;title("NNGP");
% subplot(1,2,2);imagesc(S12); colormap jet; colorbar;title("Full GP");
% 
[Lorg, Dorg] = ldl(S12);
% figure;
% subplot(2,2,1);imagesc(A2subjs); colormap jet; colorbar;title("NNGP");
% subplot(2,2,2);imagesc(Lorg-eye(nsubj*Ncoord)); colormap jet; colorbar;title("Full GP");
% subplot(2,2,3);imagesc(abs(A2subjs)>0); colormap jet; colorbar;title("NNGP");
% subplot(2,2,4);imagesc(abs(Lorg-eye(nsubj*Ncoord))>0); colormap jet; colorbar;title("Full GP");

DorgDiag = diag(Dorg);
DnngpDiag = 1./diag(D2subjsInv);

%% 3 subjects
nsubj = 3;
coord3 = coord/1.5+0.5;

K1 = exp(-pdist2(coord1, coord));
K2 = exp(-pdist2(coord2, coord));
K3 = exp(-pdist2(coord3, coord));
K12 = K1*InvCov*K2';
K13 = K1*InvCov*K3';
K23 = K2*InvCov*K3';
S3=exp(-squareform(pdist(coord3)));
Sx = exp(-squareform(pdist([coord_x, coord_y])));
S12 = [S1, K12; K12', S2];
S123 = [S1, K12, K13; K12', S2, K23; K13', K23', S3];
invS123 = inv(S123);
invS12 = inv(S12);


idx31 = knnsearch(coord1, coord3, 'K', nN);
idx32 = knnsearch(coord2, coord3, 'K', nN);

nnIDXs3 = nns_2d_square(nnDict, coord, coord3, 1, nN);
b3cov = Cov_NN_Transfer_2d(X(:),coord, nnIDXs3, 1, nN, 0);
b3 = b3cov.Bt;

cov31nngp = zeros(Ncoord);
u3SinvNNGP = zeros(Ncoord, 2*Ncoord);
cov32nngp = zeros(Ncoord);

for i=1:Ncoord
    for j=1:nN
        cov31nngp(i, idx31(i,j))=b3(i,:)*Sx(nnIDXs3.Idxs(i,:),nnIDXs1.Idxs(idx31(i,j),:))*b1(idx31(i,j), :)';
        cov32nngp(i, idx32(i,j))=b3(i,:)*Sx(nnIDXs3.Idxs(i,:),nnIDXs2.Idxs(idx32(i,j),:))*b2(idx32(i,j), :)';
    end
    u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])=...
        [cov31nngp(i,idx31(i,:)), cov32nngp(i,idx32(i,:))] /...
        (S12([idx31(i,:),Ncoord+idx32(i,:)],[idx31(i,:),Ncoord+idx32(i,:)]));
end
% 
% tic
cov3given12nngp = nan(Ncoord);
% for i=1:Ncoord
%     for j=1:Ncoord
%         cov3given12nngp(i,j)=S3(i,j)-...
%             u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
%             S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(j,:),Ncoord+idx32(j,:)])*...
%             u3SinvNNGP(j,[idx31(j,:),Ncoord+idx32(j,:)])';
%     end
% end
nn3 = nnsearch_2d(coord3, nN);
% for i=1:Ncoord
%     nnidx = nn3.Idxs(i,:);
%     nidx = nnidx(nnidx~=0);
%     cov3given12nngp(i,i)=S3(i,i)-...
%         u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
%         S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(i,:),Ncoord+idx32(i,:)])*...
%         u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])';
%     for j=1:length(nidx)
%         cov3given12nngp(i,nidx(j))=S3(i,nidx(j))-...
%             u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
%             S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])*...
%             u3SinvNNGP(nidx(j),[idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])';
%         for k=1:length(nidx)
%             if isnan(cov3given12nngp(nidx(k), nidx(j)))
%                 cov3given12nngp(nidx(k),nidx(j))=S3(nidx(k),nidx(j))-...
%                     u3SinvNNGP(nidx(k),[idx31(nidx(k),:),Ncoord+idx32(nidx(k),:)])*...
%                     S12([idx31(nidx(k),:),Ncoord+idx32(nidx(k),:)], [idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])*...
%                     u3SinvNNGP(nidx(j),[idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])';
%             end
%         end
%     end
% end
% toc

% tic 
% S3con = S3 - SPu3SinvNNGP*S12*SPu3SinvNNGP';
% toc

% figure;
% subplot(1,2,1);imagesc(cov3given12nngp); colormap jet; colorbar;title("NNGP");
% subplot(1,2,2);imagesc(S3con); colormap jet; colorbar;title("Full GP");
SPu3SinvNNGP = sparse(u3SinvNNGP);
A3 = zeros(Ncoord);
A31= zeros(Ncoord, 2*Ncoord);
A31(1,:) = SPu3SinvNNGP(1,:);
D3inv = eye(Ncoord);
i=1;
    cov3given12nngp(i,i)=S3(i,i)-...
        u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
        S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(i,:),Ncoord+idx32(i,:)])*...
        u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])';
D3inv(1,1) = 1/cov3given12nngp(1,1);

tic
for i=2:Ncoord
    nnidx = nn3.Idxs(i,:);
    nidx = nnidx(nnidx~=0);
    cov3given12nngp(i,i)=S3(i,i)-...
        u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
        S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(i,:),Ncoord+idx32(i,:)])*...
        u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])';
    for j=1:length(nidx)
        cov3given12nngp(i,nidx(j))=S3(i,nidx(j))-...
            u3SinvNNGP(i,[idx31(i,:),Ncoord+idx32(i,:)])*...
            S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])*...
            u3SinvNNGP(nidx(j),[idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])';
        for k=1:length(nidx)
            if isnan(cov3given12nngp(nidx(k), nidx(j)))
                        cov3given12nngp(nidx(k),nidx(j))=S3(nidx(k),nidx(j))-...
            u3SinvNNGP(nidx(k),[idx31(nidx(k),:),Ncoord+idx32(nidx(k),:)])*...
            S12([idx31(nidx(k),:),Ncoord+idx32(nidx(k),:)], [idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])*...
            u3SinvNNGP(nidx(j),[idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)])';
            end
        end
    end

    A3(i, nidx) = cov3given12nngp(i,nidx)/(cov3given12nngp(nidx,nidx));
    D3inv(i,i) = 1/(cov3given12nngp(i,i)-...
        cov3given12nngp(i,nidx)/...
        (cov3given12nngp(nidx,nidx))*...
        cov3given12nngp(i,nidx)');
%     A31(i,:) = A3(i,nidx)*SPu3SinvNNGP(nidx,:);
end
toc

cov3given12InvNNGP = (eye(Ncoord)-A3)'*D3inv*(eye(Ncoord)-A3);

spI = speye(Ncoord);
spA3 = sparse(A3);
tic
spA31 = (spI-spA3)*SPu3SinvNNGP;
toc
tic
A31 = (eye(Ncoord)-A3)*u3SinvNNGP;
toc

A3subjs = [A2subjs, zeros(2*Ncoord, Ncoord); A31, A3];
D3subjsInv = diag([diag(D2subjsInv);diag(D3inv)]);

cov3subjsNNGPInv = (eye(nsubj*Ncoord)-A3subjs)'*D3subjsInv*(eye(nsubj*Ncoord)-A3subjs);

figure;
subplot(1,2,1);imagesc(cov3subjsNNGPInv); colormap jet; colorbar;title("NNGP");
subplot(1,2,2);imagesc(invS123); colormap jet; colorbar;title("Full GP");

figure;
subplot(1,2,1);imagesc(inv(cov3subjsNNGPInv)); colormap jet; colorbar;title("NNGP");
subplot(1,2,2);imagesc(S123); colormap jet; colorbar;title("Full GP");


figure;
imagesc(S123-inv(cov3subjsNNGPInv)); colormap jet; colorbar;title("Full GP - NNGP");

tic
[Lorg, Dorg] = ldl(S123);
toc

figure;
subplot(2,2,1);imagesc(A3subjs); colormap jet; colorbar;title("NNGP");
subplot(2,2,2);imagesc(Lorg-eye(nsubj*Ncoord)); colormap jet; colorbar;title("Full GP");
subplot(2,2,3);imagesc(abs(A3subjs)>0); colormap jet; colorbar;title("NNGP");
subplot(2,2,4);imagesc(abs(Lorg-eye(nsubj*Ncoord))>0); colormap jet; colorbar;title("Full GP");

%%
coordY = cat(3,coord1,coord2, coord3);
bx = cat(3, b1, b2, b3);
indices1 = cat(3, idx31, idx32);
nnidxs = cat(3, nnIDXs1.Idxs, nnIDXs2.Idxs, nnIDXs3.Idxs);
Y = ones(Ncoord, 3);
tau=.1;

i=1;
j=3;
tic
Stest = CovBlock(2, nN, coordY, coord, 1,1, indices1(i,:,:), indices1(nidx(j),:,:),nnidxs, bx);
toc
tic
Stest = CovBlock_mex(2, nN, coordY, coord, 1,1, indices1(i,:,:), indices1(nidx(j),:,:),nnidxs, bx);
toc

S12test = S12([idx31(i,:),Ncoord+idx32(i,:)], [idx31(nidx(j),:),Ncoord+idx32(nidx(j),:)]);
ComparisonFig(Stest, S12test);

Sx=1*exp(-1*squareform(pdist(coord)));
S1 = 1*exp(-1*squareform(pdist(coord1)));

tic
S2Row = CovBlock_SubjectRow_mex(2,coordY,coord,1,1,tau,nnidxs, bx,Sx);
toc

S12nngp = [S1,S2Row(1:Ncoord,1:Ncoord)'; S2Row];

tic
Aout = MarginalNNGPMuS_v2_mex(3, coord, coordY, Y, 1, 1, tau, nnidxs, nn3.Idxs, nN, indices1, bx,S12nngp, true);
toc

ComparisonFig(Aout.A, [A31,A3]);


%%
close all

nN = 40;
coordLong = [coord1;coord2;coord3];
coordY = cat(3,coord1, coord2, coord3);

tic
nnLong = nnsearch_2d_mex(coordLong, nN);
toc

tic
nnLong1 = knnsearch(coordLong, coordLong, 'K', nN);
toc

nnidxs = [nnIDXs1.Idxs; nnIDXs2.Idxs; nnIDXs3.Idxs];
bx = [ b1; b2; b3];
F = [b1cov.Ft; b2cov.Ft;b3cov.Ft];

Y = ones(size(nnLong.Idxs,1),1);

idxsNN = unique(nnidxs(:));
spSx = zeros(Ncoord);
spSx(idxsNN,idxsNN) = Sx(idxsNN, idxsNN)+1e-9*eye(length(idxsNN));
spSx = sparse(spSx);
bMat = zeros(nsubj*Ncoord, Ncoord);
for i=1:(nsubj*Ncoord)
    bMat(i, nnidxs(i,:))=bx(i,:);
end
spbMat=sparse(bMat);

tic
spL = chol(spSx,'lower');
toc

K1 = exp(-pdist2(coord1, coordX));
K2 = exp(-pdist2(coord2, coordX));
K3 = exp(-pdist2(coord3, coordX));

K123 = [K1;K2;K3];

tic
Syfull = K123/Sx*K123';
toc

tic
spSy=spbMat*Sx*spbMat';
toc



tic
Sy=bMat*Sx*bMat';
toc

tic
cov3 = zeros(nsubj*Ncoord);
for i1=1:(nsubj*Ncoord)
    for i2=i1:(nsubj*Ncoord)
        cov3(i1, i2) = bx(i1,:)*Sx(nnidxs(i1,:),nnidxs(i2,:))*bx(i2,:)';
        cov3(i2, i1) = cov3(i1,i2);
    end
    cov3(i1,i1) = cov3(i1,i1)+F(i1);
end
toc

ComparisonFig(spSy, cov3);

tic
Aoutv3 = MarginalNNGPMuS_v3(Y, 1, tau, nnLong.Idxs, nnidxs, nN, bx, F, Sx, true);
toc

tic
Aoutv4 = MarginalNNGPMuS_v4_mex(Y, 1, tau, nnLong.Idxs, nN, bMat, Sx,F, true);
toc

tic
Aoutv4 = MarginalNNGPMuS_APPR_mex(Y,1, tau, nnLong.Idxs, nN, S123, true, 8);
toc



ComparisonFig(Aoutv3.A, Aoutv4.A);

ComparisonFig(Aoutv4.A, Lorg-eye(nsubj*Ncoord));
ComparisonFig(abs(Aoutv4.A)>0, abs(Lorg-eye(nsubj*Ncoord))>0,false);

cov3subjsNNGPInv = (eye(nsubj*Ncoord)-Aoutv4.A)'*Aoutv4.Dinv*(eye(nsubj*Ncoord)-Aoutv4.A);
ComparisonFig(cov3subjsNNGPInv,invS123);
ComparisonFig(inv(cov3subjsNNGPInv),S123);

figure;
subplot(2,2,1);imagesc(Aoutv4.A); colormap jet; colorbar;title("NNGP");
subplot(2,2,2);imagesc(Lorg-eye(nsubj*Ncoord)); colormap jet; colorbar;title("Full GP");
subplot(2,2,3);imagesc(abs(Aoutv4.A)>0); colormap jet; colorbar;title("NNGP");
subplot(2,2,4);imagesc(abs(Lorg-eye(nsubj*Ncoord))>0); colormap jet; colorbar;title("Full GP");


%% nearest neighbor search
alpha=1;
rho=1;
nrow = 10;
nsubj=3;
testMat = ones(nrow, nrow);
Ncoord = length(testMat(:));
[coord_x, coord_y] = find(ones(size(testMat)));
covsX = exp(-squareform(pdist([coord_x, coord_y])));
coordX = [coord_x, coord_y];
coord1 = coordX/1.5+1/5;
coord2 = coordX/3;
coord3 = coordX/1.5+0.5;
coordY = [];
coordY(:,:,1)=coord1;
coordY(:,:,2)=coord2;
coordY(:,:,3)=coord3;

coordLong = Reshape2d(coordY);
ncoord = size(coordLong, 1);
Y = ones(ncoord, 1);

tic
Sx = alpha*exp(-rho*squareform(pdist(coordX)));
toc
bi = [1,1.5,0.5]';

K1 = bi(1)*exp(-pdist2(coord1, coordX));
K2 = bi(2)*exp(-pdist2(coord2, coordX));
K3 = bi(3)*exp(-pdist2(coord3, coordX));
K12 = K1/Sx*K2';
K13 = K1/Sx*K3';
K23 = K2/Sx*K3';
S1=bi(1)^2*exp(-squareform(pdist(coord1)))+sigma(1)*eye(Ncoord);
S2=bi(2)^2*exp(-squareform(pdist(coord2)))+sigma(2)*eye(Ncoord);
S3=bi(3)^2*exp(-squareform(pdist(coord3)))+sigma(3)*eye(Ncoord);

S123 = [S1, K12, K13; K12', S2, K23; K13', K23', S3];

% 
% tic
% Sy = bMat1(1:Ncoord,:)*Sx*bMat1(1:Ncoord,:)';
% toc

% tic
% Sy = bMatFull*Sx*bMatFull';
% toc

% 
% tic
% SyAppr = Sapproximate_mex(bs,alpha,.1, nnidxs,Sx,Ncoord, 3, 8);
% toc
% 
% 
% tic
% SyAppr = Sapproximate(bMat1,alpha,rho,.1,coordLong, nnidxs,Sx,Ncoord, nsubj, 8);
% toc
% 
% ComparisonFig(SyAppr, Sy);
% % 
% 
% idxsNN = unique(nnidxs(:));
% spSx = zeros(Ncoord);
% spSx(idxsNN,idxsNN) = Sx(idxsNN, idxsNN);
% spSx = sparse(spSx);
% 
% spy(spSx);
% rp = symrcm(spSx);
% 
% spy(spSx(rp, rp))
% 
% % figure; spy(bMat1(:,rp));
% 
% tic
% spSy=bMat1*spSx*bMat1';
% toc
% 
% tic
% ind = 1:Ncoord;
% spSy=bMat1(ind,rp(1:5000))*spSx(rp(1:5000),rp(1:5000))*bMat1(ind,rp(1:5000))';
% toc
% 
% ComparisonFig(spSy, Sy);


% tic
% Aoutv5 = MarginalNNGPMuS_APPR_mex(Y, 1, 0.1,nnidxAPPR, nN, Sy, false, 8);
% Aoutv5.loglik
% toc


sigma = 0.0001*ones(nsubj,1);
K1 = 5;
K2 = 20;
nparts=3;
tic
[nnidxs, bs, ~] = nns_2d_square_parallel(coordX, coordLong, Sx, ncoord, nsubj, K1, 1,1,1, false);
nnidxAPPR = nnsearch_approx(coordLong,K2,ncoord/nparts, nparts, int32(16));
Aoutv1 = MarginalNNGPMuS_APPR_full(Y, coordLong,nsubj,bi,1,1,sigma, nnidxAPPR, K2, Sx, bs, nnidxs, true, 2);
Aoutv1.loglik
toc

tic
[nnidxs, bs, ~] = nns_2d_square_parallel(coordX, coordLong, Sx, ncoord, nsubj, K1, 1,1,1, false);
nnidxAPPR = nnsearch_approx(coordLong,K2,ncoord/nparts, nparts, int32(16));
Aoutv1 = MarginalNNGPMuS_APPR_v1(Y,nsubj,bi,1,0.1*ones(nsubj,1), nnidxAPPR, K2, Sx, bs, nnidxs, true, 2);
Aoutv1.loglik
toc

tic
llk = posteriorLog(Y,coordY, coordX, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma, Sx, 8);
llk
toc

tic
llk1 = posteriorLog_full(Y,coordY, coordX, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma, Sx, 8);
llk1
toc

truellk = log(mvnpdf(Y, 0, S123))

[Lorg, Dorg] = ldl(S123);
ComparisonFig(Aoutv1.A, Lorg-eye(nsubj*Ncoord));
S_NNGPInv = (eye(nsubj*Ncoord)-Aoutv1.A)'*Aoutv1.Dinv*(eye(nsubj*Ncoord)-Aoutv1.A);
SyApprInv = inv(S123);
% ComparisonFig(S_NNGPInv, SyApprInv);
ComparisonFig(inv(S_NNGPInv), S123);

llk1 = log(mvnpdf(Y, Aoutv1.A*Y, inv(Aoutv1.Dinv)))


figure;
subplot(2,2,1);imagesc(Aoutv1.A); colormap jet; colorbar;title("NNGP");
subplot(2,2,2);imagesc(Lorg-eye(nsubj*Ncoord)); colormap jet; colorbar;title("Full GP");
subplot(2,2,3);imagesc(abs(Aoutv1.A)>0); colormap jet; colorbar;title("NNGP");
subplot(2,2,4);imagesc(abs(Lorg-eye(nsubj*Ncoord))>0); colormap jet; colorbar;title("Full GP");




