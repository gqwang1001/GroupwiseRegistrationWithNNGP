function out = FullMarginalLL(Y, coordY, coord, alpha, rho, sigma, b, Sx, nsubj, ncoord)

Ncrd = ncoord/nsubj;

Sy = eye(ncoord);
CovYX = zeros(ncoord, Ncrd);

for i=1:nsubj
    icoord = (Ncrd*(i-1)+1):(Ncrd*i);
    CovYX(icoord,:) = b(i)*alpha*exp(-rho*pdist2(coordY(:,:,i), coord));
end

for i=1:nsubj
    icoord = (Ncrd*(i-1)+1):(Ncrd*i);
    Sy(icoord,icoord)=b(i)^2*alpha*exp(-rho*squareform(pdist(coordY(:,:,i))))+sigma(i)*eye(Ncrd);
    for j=(i+1):nsubj
        jcoord = (Ncrd*(j-1)+1):(Ncrd*j);
        Sy(icoord, jcoord) = CovYX(icoord,:)/Sx*CovYX(ijcoord,:)';
    Sy(jcoord, icoord) = Sy(icoord, jcoord);
    end
end

out = log(mvnpdf(Y,0,Sy));

end