function array2d = Reshape2d(array3d)
C = permute(array3d,[1 3 2]);
array2d = reshape(C,[],size(array3d, 2),1);
end

