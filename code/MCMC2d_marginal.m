function out = MCMC2d_marginal(init, dataIn, niter)

tstart = tic;
nsubj = dataIn.nsubj;
coord = dataIn.coord;
nparts = dataIn.nparts;
parworkers = 8;
K1= dataIn.K1;
K2 = dataIn.K2;
rho = init.rho;
alpha = init.alpha;
Y = dataIn.Y(:);
b = init.b;
sigma = init.sigma;
bnds = dataIn.bnds;
t0 = [init.scaling, init.rot', init.shift];
Sx = alpha*exp(-rho*squareform(pdist(coord)));

coordY=zeros(size(coord,1),2,nsubj);

for subj = 1:nsubj
    coordY(:,:,subj) = Transfer_Coord(init.scaling(subj, :), init.rot(subj), init.shift(subj,:), coord);
end

tMat = zeros(niter,nsubj,size(t0,2));
tMat_adaptLog = zeros(niter,nsubj,size(t0,2)+2);
alphaMat = zeros(niter, 1);
rhoMat = zeros(niter, 1);
acceptMat = zeros(niter, nsubj, size(t0,2)+2);
acceptMatAR = zeros(niter, 2);

ropt = 0.234;
rhat = zeros(nsubj, size(t0,2)+2);
adapt = true;
k=50;
c0=1;
c1=0.8;
logsig = repmat(log(2.4^2/2), nsubj, 5+2);
Sig = ones(5+2,5+2,nsubj);
Sig(4:5, 4:5, :) = ones(2,2,nsubj);

logsigAR = (log(2.4^2)-1)*ones(2,1);
SigAR = eye(2);
rng('default');

for iter = 1:niter
    tstart1 = tic;
    tupdate = transfUpdate_2d_marginal_full_mex(t0, b, Y, Sx, coordY, sigma, alpha, rho, coord, K1, K2, ...
        logsig, Sig,logsigAR, SigAR, nsubj,nparts, bnds, parworkers);
    
    Sx = tupdate.Sx;
    t0 = tupdate.transf;
    b = tupdate.bi;
    alpha = tupdate.alpha;
    rho = tupdate.rho;
    sigma = tupdate.sigma;
    coordY = tupdate.coordY;
    
    tMat(iter,:,:) = tupdate.transf;
    acceptMat(iter,:,:) = tupdate.accept;
     
    tMat_adaptLog(iter,:, 1:2) = log(tupdate.transf(:, 1:2));
    tMat_adaptLog(iter,:, 3) = tupdate.transf(:,3)./180*pi;
    tMat_adaptLog(iter,:,4:5) = tupdate.transf(:,4:5);
    tMat_adaptLog(iter,:,6) = log(tupdate.bi);
    tMat_adaptLog(iter,:,7) = log(tupdate.sigma);
    
    alphaMat(iter) = tupdate.alpha;
    rhoMat(iter) = tupdate.rho;
    
    acceptMatAR(iter,:) = tupdate.accept_alphaRho;
    disp(iter);
    disp(toc(tstart1));
    
    % adaptive MCMC
    if (adapt && rem(iter, k)==0)
        disp(iter);     disp(toc(tstart));
        gamma1 = 1/((floor(iter/k)+1)^c1);
        gamma2 = c0*gamma1;
        for subj = 1:nsubj
            Sig0tHat = cov(squeeze(tMat_adaptLog((iter-k+1):iter,subj,:)));
            for j = 1:(size(t0,2)+2)
                rhat(subj, j) = mean(acceptMat((iter-k+1):iter, subj, j));
                logsig(subj, j) = logsig(subj, j) + gamma2 *(rhat(subj, j)-ropt);
            end
            Sig(:,:,subj) = Sig(:,:,subj)+gamma1*(Sig0tHat-Sig(:,:,subj));
        end
        
        %         SigHatAR = cov(log([alphaMat((iter-k+1):iter),rhoMat((iter-k+1):iter)]));
        rhatAlpha = mean(acceptMatAR((iter-k+1):iter, 1));
        rhatRho = mean(acceptMatAR((iter-k+1):iter, 2));
        
        logsigAR(1) = logsigAR(1)+gamma2*(rhatAlpha-ropt);
        logsigAR(2) = logsigAR(2)+gamma2*(rhatRho-ropt);
        %         SigAR = SigAR+gamma1*(SigHatAR-SigAR);
    end
end

out.Transf = tMat;
out.b = exp(squeeze(tMat_adaptLog(:,:,6)));
% out.X = XMat;
out.sigma = exp(squeeze(tMat_adaptLog(:,:,7)));
out.alpha = alphaMat;
out.rho = rhoMat;

out.time = toc(tstart);
end



