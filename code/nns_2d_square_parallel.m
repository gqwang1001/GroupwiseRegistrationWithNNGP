% function [nnIdx, Bt] = nns_2d_square_parallel(nnDict,coord, coordTrans, incre, K, rho,  parworkers)
function [nnIdx, Bt, BtMat] = nns_2d_square_parallel(coord, coordTrans, n_coord, nsubj, K, rho, parworkers)

% n_coord = size(coordTrans, 1);
Bt = zeros(n_coord, K);
% Ft = zeros(n_coord, 1);
% [minco, maxco] = bounds(nnDict.newcoords);
% ImgSz = [(maxco(1)-minco(1))/incre+1, (maxco(2)-minco(2))/incre+1];
% nnIdx = zeros(n_coord, K);
% % nnDists = zeros(n_coord, K);
% round_coordTrans = round(coordTrans);
% 
% parfor (i=1:n_coord, parworkers)
%     
%     if round_coordTrans(i,1) < minco(1)
%         
%         if round_coordTrans(i,2) < minco(2)
%             idx = sub2ind(ImgSz, 1, 1);
%         elseif round_coordTrans(i,2) <= maxco(2)4
%         if round_coordTrans(i,2) < minco(2)
%             idx = sub2ind(ImgSz, floor((round_coordTrans(i,1)-minco(1))/incre)+1, 1);
%         elseif round_coordTrans(i,2) <= maxco(2)
%             idx = sub2ind(ImgSz, floor((round_coordTrans(i,1)-minco(1))/incre)+1, floor((round_coordTrans(i,2)-minco(2))/incre)+1);
%         else
%             idx = sub2ind(ImgSz, floor((round_coordTrans(i,1)-minco(1))/incre)+1, floor((maxco(2)-minco(2))/incre)+1);
%         end
%         
%     else
%         
%         if round_coordTrans(i,2) < minco(2)
%             idx = sub2ind(ImgSz, floor((maxco(1)-minco(1))/incre)+1, 1);
%         elseif round_coordTrans(i,2) <= maxco(2)
%             idx = sub2ind(ImgSz, floor((maxco(1)-minco(1))/incre)+1, floor((round_coordTrans(i,2)-minco(2))/incre)+1);
%         else
%             idx = sub2ind(ImgSz, floor((maxco(1)-minco(1))/incre)+1, floor((maxco(2)-minco(2))/incre)+1);
%         end
%     end
%     % Summarize nn info and tranfermation map
%     nnIdx(i,:) = nnDict.idx(idx,:);
%     nnDists= sqrt(sum(([coordTrans(i,1).*ones(K,1), coordTrans(i,2).*ones(K,1)] - coord(nnIdx(i,:),:)).^2, 2));
%     Bt(i,:) = exp(-rho*nnDists')/exp(-rho*squareform(pdist(coord(nnIdx(i,:), :))));
%     
% end

[nnIdx, nnDist] = knnsearch(coord, coordTrans, 'K', K);
covMat = exp(-rho*nnDist);

parfor (i=1:n_coord, parworkers)
    Bt(i,:) = covMat(i,:)/exp(-rho*squareform(pdist(coord(nnIdx(i,:),:))));
%     Ft(i) = alpha+sigma-alpha*Bt(i,:)*covMat(i,:)';
end

% BtMat = zeros(n_coord, n_coord/nsubj);
% for i=1:n_coord
%      BtMat(i, nnIdx(i,:)) = Bt(i,:);
% end

BtMat=sparse(repmat((1:n_coord)', K,1), nnIdx(:), Bt(:), n_coord, n_coord/nsubj);


end