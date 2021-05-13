function K = calculateEuclideanDist(P,Q,N)
% Vectorized method to compute pairwise Euclidean distance
% Returns K(i,j) = sqrt((P(i,:) - Q(j,:))'*(P(i,:) - Q(j,:)))

    K = zeros(N,1,'double');
    for q=1:2
        K = K+(P(:,q)-Q(:,q)).^2;
    end
end