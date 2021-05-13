function ComparisonFig(MatNNGP,Mat, constrast)

if nargin<3
    constrast=true;
end

if constrast
    figure;
    subplot(1,3,1);imagesc(MatNNGP); colormap jet; colorbar;title("NNGP");
    subplot(1,3,2);imagesc(Mat); colormap jet; colorbar;title("Full GP");
    subplot(1,3,3);imagesc(MatNNGP-Mat); colormap jet; colorbar;title("NNGP-Full GP");
else
    figure;
    subplot(1,2,1);imagesc(MatNNGP); colormap jet; colorbar;title("NNGP");
    subplot(1,2,2);imagesc(Mat); colormap jet; colorbar;title("Full GP");
end

end

