function out = transfUpdate_2d_marginal_full_FULLGP(transf0, bi, Y, Sx, coordY, sigma, alpha, rho, coord, K1,K2, logsig, Sig,logAR, SigAR, nsubj,nparts, bnds, parworkers)

A  = zeros(nsubj,1);
nparams_tf = size(transf0, 2);
% nparams_else = 2;
ncoord = size(Y,1);
accept = zeros(nsubj, nparams_tf);

t1 = transf0;
coordY0 = coordY;
logliks = posteriorLog_full(Y, coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);

for i = 1:nsubj
    
    tOut = t1(i,:);
    
    accept_i = zeros(1, nparams_tf);
    % scaling X
    tfLogProp = randn(1) * sqrt(exp(logsig(i,1))*Sig(1,1,i)) + log(tOut(1));
    t2 = transf0(:,1);
    t2(i) = exp(tfLogProp);
    ScalingX = t2(i);
    %     t2Centered = t2./mean(t2);
    %     ScalingX = t2Centered(i);
    
    if ScalingX>bnds(1,1) && ScalingX<bnds(1,2)
        
        tf1Prop = [exp(tfLogProp), tOut(2:5)];
        coordY(:,:,i) = Transfer_Coord(tf1Prop(1:2), tf1Prop(3), tf1Prop(4:5), coord);
        logliks1 = posteriorLog_full(Y,coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
        A(i) = min(0, logliks1-logliks+log(normpdf(log(tf1Prop(1)), 0, 1))- log(normpdf(log(tOut(1)), 0, 1)));
        
        if (log(rand(1)) < A(i))
            tOut = tf1Prop;
            accept_i(1)  = 1;
            logliks = logliks1;
        end
    end
    
    % scaling Y
    tfLogProp = randn(1) * sqrt(exp(logsig(i,2))*Sig(2,2,i)) + log(tOut(2));
    
    t2 = transf0(:,2);
    t2(i) = exp(tfLogProp);
    %     t2Centered = t2./mean(t2);
    %     ScalingY = t2Centered(i);
    ScalingY = t2(i);
    if ScalingY>bnds(1,1) && ScalingY<bnds(1,2)
        
        tf1Prop = [tOut(1) exp(tfLogProp) tOut(3:5)];
        coordY(:,:,i) = Transfer_Coord(tf1Prop(1:2), tf1Prop(3), tf1Prop(4:5), coord);
        logliks1 = posteriorLog_full(Y,coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
        A(i) = min(0, logliks1-logliks+log(normpdf(log(tf1Prop(2)), 0, 1))- log(normpdf(log(tOut(2)), 0, 1)));
        
        if (log(rand(1)) < A(i))
            tOut = tf1Prop;
            accept_i(2)  = 1;
            logliks = logliks1;
        end
    end
    % rotation
    tfProp = randn(1) * sqrt(exp(logsig(i,3))*Sig(3,3,i)) + tOut(3)/180*pi;
    
    t2 = transf0(:,3);
    t2(i) = tfProp/pi*180;
    %     t2Centered = t2-mean(t2);
    %     Rot = t2Centered(i);
    Rot = t2(i);
    
    if Rot>bnds(2,1) && Rot<bnds(2,2)
        
        tf1Prop = [tOut(1:2), tfProp/pi*180, tOut(4:5)];
        coordY(:,:,i) = Transfer_Coord(tf1Prop(1:2), tf1Prop(3), tf1Prop(4:5), coord);
        logliks1 = posteriorLog_full(Y,coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
        
        A(i) = min(0, logliks1-logliks+log(normpdf(tf1Prop(3)/180*pi,0,.2))- log(normpdf(tOut(3)/180*pi,0,.2)));
        if (log(rand(1)) < A(i))
            tOut = tf1Prop;
            accept_i(3)  = 1;
            logliks = logliks1;
        end
    end
    
    % shifting X
    tfProp = randn(1) * sqrt(exp(logsig(i,4))*Sig(4,4,i)) + tOut(4);
    t2 = transf0(:,4);
    t2(i) = tfProp;
    %     t2Centered = t2-mean(t2);
    %     ShiftX = t2Centered(i);
    ShiftX = t2(i);
    if ShiftX>bnds(3,1) && ShiftX<bnds(3,2)
        
        tf1Prop = [tOut(1:3) tfProp tOut(5)];
        coordY(:,:,i) = Transfer_Coord(tf1Prop(1:2), tf1Prop(3), tf1Prop(4:5), coord);
        logliks1 = posteriorLog_full(Y,coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
        
        A(i) = min(0, logliks1-logliks+ log(normpdf(tf1Prop(4), 0, 3))- log(normpdf(tOut(4), 0, 3)));
        
        if (log(rand(1)) < A(i))
            tOut = tf1Prop;
            accept_i(4)  = 1;
            logliks = logliks1;
        end
    end
    % shifting Y
    tfProp = randn(1) * sqrt(exp(logsig(i,5))*Sig(5,5,i)) + tOut(5);
    t2 = transf0(:,5);
    t2(i) = tfProp;
    %     t2Centered = t2-mean(t2);
    %     ShiftY = t2Centered(i);
    ShiftY = t2(i);
    
    if ShiftY>bnds(3,1) && ShiftY<bnds(3,2)
        
        tf1Prop = [tOut(1:4) tfProp];
        coordY(:,:,i) = Transfer_Coord(tf1Prop(1:2), tf1Prop(3), tf1Prop(4:5), coord);
        logliks1 = posteriorLog_full(Y, coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
        A(i) = min(0, logliks1-logliks+log(normpdf(tf1Prop(5), 0, 3))- log(normpdf(tOut(5), 0, 3)));
        
        if (log(rand(1)) < A(i))
            tOut = tf1Prop;
            accept_i(5) = 1;
            logliks = logliks1;
        end
    end
    
    t1(i,:) = tOut;
    accept(i, :) = accept_i;
    
end


% centerize the transformation parameters
tCentered = [t1(:,1)./ mean(t1(:,1)),...
    t1(:,2)./ mean(t1(:,2)),...
    (t1(:,3) - mean(t1(:,3))),...
    (t1(:,4) - mean(t1(:,4))),...
    (t1(:,5) - mean(t1(:,5)))];

% acceptOut = zeros(nsubj, 5);

% tTreshed_scalingX = tCentered(:,1)<0;
% tTreshed_scalingY = tCentered(:,2)<0;
% tOut = transf0;
coordY = coordY0;

% if sum(tTreshed_scalingX+tTreshed_scalingY) == 0
tOut= tCentered;
acceptOut = accept;
for i=1:nsubj
    coordY(:,:,i) = Transfer_Coord(tOut(i,1:2), tOut(i,3), tOut(i,4:5), coord);
end
% end

logliks = posteriorLog_full(Y, coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
acceptOut_b=zeros(nsubj,1);
% intensity correction
for i=1:nsubj
    b1 = bi;
    bpropLog =  randn(1)*sqrt(exp(logsig(i,nparams_tf+1))*Sig(nparams_tf+1,nparams_tf+1,i))+log(bi(i));
    b1(i) = exp(bpropLog);
    %     bCentered = b1/mean(b1);
    %     if bCentered(i)>bnds(4,1) && bCentered(i)<bnds(4,2)
    %     if b1(i)>bnds(4,1) && b1(i)<bnds(4,2)
    logliks1 = posteriorLog_full(Y, coordY, coord, b1, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
    A(i) = min(0, logliks1-logliks);
    if log(rand(1)) < A(i)
        bi = b1;
        acceptOut_b(i) = 1;
        logliks = logliks1;
    end
    %     end
end
bCentered = bi/mean(bi);

% if bCentered~=bi
logliks = posteriorLog_full(Y, coordY, coord, bCentered, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers);
% end

% subject specific measurement error
acceptOut_Sigma=zeros(nsubj,1);
for i=1:nsubj
    s1 = sigma;
    spropLog = randn(1)*sqrt(exp(logsig(i,nparams_tf+2))*Sig(nparams_tf+2,nparams_tf+2,i))+log(sigma(i));
    s1(i) = exp(spropLog);
    if s1(i)<1 && s1(i)>1e-6
        logliks1 = posteriorLog_full(Y, coordY, coord, bCentered, ncoord, nsubj, nparts, K1, K2, alpha, rho, s1,Sx, parworkers);
        A(i) = min(0, logliks1-logliks-spropLog+log(sigma(i)));
        if log(rand(1)) < A(i)
            sigma = s1;
            acceptOut_Sigma(i)=1;
            logliks = logliks1;
        end
    end
end


% alpha update

accept_alphaRho = zeros(2,1);
lambda = exp(logAR(1));
alpha1 = alpha + lambda * (-1 + 2* rand(1));
p01 = 0; % pdf of rho0 given rho1
p10 = 0; % pdf of rho1 given rho0
if (alpha>bnds(5,1)) && (alpha<bnds(5,2))
    p01 = unifpdf(alpha, alpha1-lambda, alpha1+lambda);
end
if (alpha1>bnds(5,1)) && (alpha1<bnds(5,2))
    p10 = unifpdf(alpha1, alpha-lambda, alpha+lambda);
end
if p10*p01~=0
    Sx1 = Sx/alpha*alpha1;
    logliks1 = posteriorLog_full(Y, coordY, coord, bCentered, ncoord, nsubj, nparts, K1, K2, alpha1, rho, sigma,Sx1, parworkers);
    A = exp(logliks1-logliks)*(alpha/alpha1)*p01/p10;
    if rand(1) < A
        alpha = alpha1;
        accept_alphaRho(1,1) = 1;
        Sx = Sx1;
        logliks = logliks1;
    end
end

% rho update
lambda = exp(logAR(2));
rho1 = rho+ lambda*(-1 + 2* rand(1));
p01 = 0; % pdf of rho0 given rho1
p10 = 0; % pdf of rho1 given rho0

if (rho>bnds(6,1)) && (rho<bnds(6,2))
    p01 = unifpdf(rho, rho1-lambda, rho1+lambda);
end

if (rho1>bnds(6,1)) && (rho1<bnds(6,2))
    p10 = unifpdf(rho1, rho-lambda, rho+lambda);
end

if p10*p01~=0
    Sx1 = alpha*exp(rho1*log(Sx/alpha)/rho);
    logliks1 = posteriorLog_full(Y, coordY, coord, bCentered, ncoord, nsubj, nparts, K1, K2, alpha, rho1, sigma,Sx1, parworkers);
    A = exp(logliks1-logliks)*p01/p10;
    if rand(1) < A
        rho = rho1;
        accept_alphaRho(2,1) = 1;
        Sx = Sx1;
        %     logliks = logliks1;
    end
end

out.transf = tOut;
out.accept = [acceptOut, acceptOut_b, acceptOut_Sigma];
out.bi = bCentered;
out.sigma = sigma;
out.accept_alphaRho = accept_alphaRho;
out.alpha = alpha;
out.rho = rho;
out.Sx = Sx;
end