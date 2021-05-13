function out = MarginalNNGPMuS_APPR_v1(Y, alpha, sigma, nnLong, K,Sx,bs,nnidxs, choleskyMat, parworkers)

% Estimated Marginal NNGP mean and Variance with Approximate marginal
% covariance matrix

% Sy: Covariance matrix of Y.
% K: nearest neighbor size
% nnLong: V*Nsubj-by-K-by nearest neighor info of isubj to 1,...,isubj-1
% Y: V*Nsubj stacked images.
% choleskyMat: logic. whether or not produce cholesky factors
% bMat: V*Nsubj-by-V sparse matrix
% alpha: spatial variance
% sigma: nugget effect.

Ncoord = size(Y,1);
mu = zeros(Ncoord,1);
bcov = zeros(Ncoord, K);
Fdiag = (alpha+sigma)*ones(Ncoord,1);
%
% Sy=bMat*Sx*bMat';
% V = Ncoord/Nsubj;
% for i=1:Nsubj
%     Sy((1:V)+(i-1)*V, (1:V)+(i-1)*V)=alpha*exp(-rho*squareform(pdist(coordY(:,:,i))))+sigma*eye(V);
% end

%i==2
idx = nnLong(2,1);
cov1 = SapprSubCov(bs, alpha, sigma,nnidxs, Sx, 2, idx);
covMat = alpha+sigma;
bcov(2, 1)= cov1/covMat;
Fdiag(2) = alpha+sigma-bcov(2,1)*cov1';
mu(2) = bcov(2,1)*Y(idx);


for i=3:K
    idx = nnLong(i,1:(i-1));
    cov1 = SapprSubCov(bs, alpha, sigma,nnidxs, Sx, i, idx);
    covMat = SapprSubCov(bs, alpha, sigma,nnidxs, Sx, idx, idx);
    bcov(i, 1:(i-1))= cov1/covMat;
    Fdiag(i) = alpha+sigma-bcov(i,1:(i-1))*cov1';
    mu(i) = bcov(i,1:(i-1))*Y(idx);
end

parfor (i=(K+1):Ncoord,parworkers)
    idx = nnLong(i,:);
    cov1 = SapprSubCov(bs, alpha, sigma,nnidxs, Sx, i, idx);
    covMat = SapprSubCov(bs, alpha, sigma,nnidxs, Sx, idx, idx);
    bcov(i,:)=  cov1/covMat;
    Fdiag(i) = alpha+sigma-bcov(i,:)*cov1';
    mu(i) = bcov(i,:)*Y(nnLong(i,:));
end

Aout = nan;
Dinv = nan;

out.bx = bcov;
out.Fs = Fdiag;
out.mu = mu;
out.loglik = logNormalPdf(Y, mu, sqrt(Fdiag));

if choleskyMat
    Aout = zeros(Ncoord);
    Dinv = eye(Ncoord);
    for i = 2:Ncoord
        idxs = nnLong(i,:);
        idx = idxs(idxs~=0);
        nidx = length(idx);
        Aout(i,idx) = bcov(i,1:nidx);
        Dinv(i,i) = 1/Fdiag(i);
    end
end

out.A = Aout;
out.Dinv = Dinv;

end