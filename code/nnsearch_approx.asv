function outIdx= nnsearch_approx(coordLong, K, ncoord, nsubj, parworkers)
%Directed Nearest Neighbor: fast algorithm to find the directed nearest
%neighbor

% coordLong: V-by-2 coordinates
% K: number of nn
% nsubj: number of subjects
% ncoord: number of coordinates per subjects
% parworkers: number of core utilized for parrallel computing

outIdx3d = zeros(ncoord, K, nsubj);
parfor (j =1:nsubj,parworkers)
    nnidxsBtwsubj = zeros(ncoord, K);
    if j>1
        nnidxsBtwsubj=knnsearch(coordLong(1:((j-1)*ncoord),:), coordLong((1:ncoord)+(j-1)*ncoord,:), 'K', K);
    end
    for i=1:ncoord
        if i<=K && j==1
            outIdx3d(i,:,j) = [1:(i-1), zeros(1,K-i+1)];
        elseif j==1
            dists = calculateEuclideanDist(coordLong(i,:), coordLong(1:(i-1),:) , i-1);
            [~,I]=mink(dists, K);
            outIdx3d(i,:,j) = I;
        else
            ids = [nnidxsBtwsubj(i,:),(1:(i-1))+(j-1)*ncoord];
            dists = calculateEuclideanDist(coordLong(i+(j-1)*ncoord,:), coordLong(ids,:), i-1+K);
            [~,I] = mink(dists,K);
            outIdx3d(i,:,j) = ids(I);
        end
    end
end
outIdx3d(1,1,1) = 1;

% reshape 3d array
C = permute(outIdx3d,[1 3 2]);
outIdx = reshape(C,[],K,1);

end