function cov = SapprSubCov_full(coordLong, ncoord, nsubj, bs, bvec, alpha, rho, sigma, nnidxs, Sx, idx1, idx2)

m1= length(idx1);
m2= length(idx2);
cov = zeros(m1, m2);
if idx1==idx2
    %     cov = zeros(m1);
    for i1=1:m1
        isubj1 = floor((idx1(i1)-1)/ncoord*nsubj)+1;
        for i2=(i1+1):m2
            isubj2 = floor((idx2(i2)-1)/ncoord*nsubj)+1;
            if isubj1==isubj2
                cov(i1,i2) = bvec(idx1(i1))^2*alpha*exp(-rho*sqrt(sum((coordLong(idx1(i1),:)-coordLong(idx2(i2),:)).^2)));
            else
                cov(i1,i2) = bvec(idx1(i1))*bs(idx1(i1),:)*Sx(nnidxs(idx1(i1),:),nnidxs(idx2(i2),:))*bs(idx2(i2),:)'*bvec(idx2(i2));
            end
            cov(i2,i1) = cov(i1,i2);
        end
        cov(i1,i1) = bvec(idx1(i1))^2*alpha+sigma(idx1(i1));
    end
else
    %     cov = zeros(1,m2);
    isubj1 = floor((idx1(1)-1)/ncoord*nsubj)+1;
    for i2=1:m2
        isubj2 = floor((idx2(i2)-1)/ncoord*nsubj)+1;
        if isubj1==isubj2
            cov(1,i2) = bvec(idx1(1))^2*alpha*exp(-rho*sqrt(sum((coordLong(idx1(1),:)-coordLong(idx2(i2),:)).^2)));
        else
            cov(1,i2) = bvec(idx1(1))*bs(idx1(1),:)*Sx(nnidxs(idx1(1),:),nnidxs(idx2(i2),:))*bs(idx2(i2),:)'*bvec(idx2(i2));
        end
    end
end

end

