function cov = Sapproximate(bs,alpha,rho,sigma,coordLong,nnidxs,Sx,Ncoord, nsubj, parworkers)


V = nsubj*Ncoord;
cov = 0.5*(alpha+sigma)*eye(nsubj*Ncoord);

% parfor (i1=1:V, parworkers)
%     for i2=1:V
%         if i1>i2
%             cov(i1, i2) = bs(i1,:)*Sx(nnidxs(i1,:),nnidxs(i2,:))*bs(i2,:)';
%         end
%     end
% end
%
% cov = cov + cov';

for i=1:nsubj
    for j=1:nsubj
        if i==j
            cov((1:Ncoord)+(i-1)*Ncoord, (1:Ncoord)+(i-1)*Ncoord)=...
                alpha*exp(-rho*squareform(pdist(coordLong((1:Ncoord)+(i-1)*Ncoord,:))))+sigma*eye(Ncoord);
        elseif i>j
            cov((1:Ncoord)+(i-1)*Ncoord, (1:Ncoord)+(j-1)*Ncoord)=...
                bs((1:Ncoord)+(i-1)*Ncoord,:)*Sx*bs((1:Ncoord)+(i-1)*Ncoord,:)';
            %             cov((1:Ncoord)+(i-1)*Ncoord, (1:Ncoord)+(j-1)*Ncoord)=...
            %                 Sapproximate_mex(bs((1:Ncoord)+(i-1)*Ncoord,:),alpha, .1, nnidxs((1:Ncoord)+(i-1)*Ncoord,:),Sx,Ncoord,1,parworkers);
        end
        cov((1:Ncoord)+(i-1)*Ncoord, (1:Ncoord)+(j-1)*Ncoord)=...
            cov((1:Ncoord)+(j-1)*Ncoord, (1:Ncoord)+(i-1)*Ncoord);
    end
end

end

