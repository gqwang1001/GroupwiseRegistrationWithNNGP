function out = MarginalNNGPMuS_APPR_full(Y, coordLong, nsubj, bi, alpha,rho, sigma, nnLong, K, Sx, bs, nnidxs, choleskyMat, parworkers)

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
sigmaMat = repmat(sigma',Ncoord/nsubj,1);
sigmaLong = sigmaMat(:);

biMat = repmat(bi', Ncoord/nsubj,1);
bvec = biMat(:);
Fdiag = alpha*(bvec.^2)+sigmaLong;

%i==2
idx = nnLong(2,1);
% cov1 = SapprSubCov(bs, bvec, alpha, sigmaLong, nnidxs, Sx, 2, idx);
cov1 = SapprSubCov_full(coordLong, Ncoord, nsubj, bs, bvec, alpha, rho, sigmaLong, nnidxs, Sx, 2, idx);

covMat = bvec(idx).^2*alpha+sigmaLong(idx);
bcov(2, 1)= cov1/covMat;
Fdiag(2) = Fdiag(2)-bcov(2,1)*cov1';
mu(2) = bcov(2,1)*Y(idx);


for i=3:K
    idx = nnLong(i,1:(i-1));
    % cov1 = SapprSubCov(bs, bvec, alpha, sigmaLong, nnidxs, Sx, i, idx);
    cov1 = SapprSubCov_full(coordLong, Ncoord, nsubj, bs, bvec, alpha, rho, sigmaLong, nnidxs, Sx, i, idx);
    % covMat = SapprSubCov(bs,bvec, alpha, sigmaLong, nnidxs, Sx, idx, idx);
    covMat = SapprSubCov_full(coordLong, Ncoord, nsubj,bs,bvec, alpha, rho, sigmaLong, nnidxs, Sx, idx, idx);
    bcov(i, 1:(i-1))= cov1/covMat;
    Fdiag(i) = Fdiag(i)-bcov(i,1:(i-1))*cov1';
    mu(i) = bcov(i,1:(i-1))*Y(idx);
end

parfor (i=(K+1):Ncoord,parworkers)
    idx = nnLong(i,:);
    % cov1 = SapprSubCov(bs, bvec, alpha, sigmaLong, nnidxs, Sx, i, idx);
    cov1 = SapprSubCov_full(coordLong, Ncoord, nsubj, bs, bvec, alpha, rho, sigmaLong, nnidxs, Sx, i, idx);
    % covMat = SapprSubCov(bs,bvec, alpha, sigmaLong, nnidxs, Sx, idx, idx);
    covMat = SapprSubCov_full(coordLong, Ncoord, nsubj,bs,bvec, alpha, rho, sigmaLong, nnidxs, Sx, idx, idx);
    
    bcov(i,:)=  cov1/covMat;
    Fdiag(i) = Fdiag(i)-bcov(i,:)*cov1';
    mu(i) = bcov(i,:)*Y(idx);
end

Aout = nan;
Dinv = nan;

out.bx = bcov;
out.Fs = Fdiag;
out.mu = mu;
% out.loglik = logNormalPdf(Y, mu, sqrt(Fdiag));
out.loglik = sum(log(normpdf(Y, mu, sqrt(Fdiag))));

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