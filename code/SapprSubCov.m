function cov = SapprSubCov(bs, bvec, alpha, sigma, nnidxs, Sx, idx1, idx2)

m1= length(idx1);
m2= length(idx2);
cov = zeros(m1, m2);
if idx1==idx2
%     cov = zeros(m1);
    for i1=1:m1
        for i2=(i1+1):m2
            cov(i1,i2) = bvec(idx1(i1))*bs(idx1(i1),:)*Sx(nnidxs(idx1(i1),:),nnidxs(idx2(i2),:))*bs(idx2(i2),:)'*bvec(idx2(i2));
            cov(i2,i1) = cov(i1,i2);
        end
        cov(i1,i1) = bvec(idx1(i1))^2*alpha+sigma(idx1(i1));
    end
    
else
%     cov = zeros(1,m2);
    for i2=1:m2
        cov(1, i2) = bvec(idx1(1))* bs(idx1(1),:)*Sx(nnidxs(idx1(1),:),nnidxs(idx2(i2),:))*bs(idx2(i2),:)'*bvec(idx2(i2));
    end
end

end

