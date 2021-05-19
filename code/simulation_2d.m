% load MINST data
[XTrain,YTrain] = digitTrain4DArrayData;

%%
Imgs = [];
imgsInput = [];

for i = 1:3
    temp = XTrain(:,:,:,i+3);
    Imgs(:,:,i) = temp;
    imgsInput(:,i) = temp(:);
end
[coord_x, coord_y] = find(ones(size(temp)));

dat = struct();
dat.XInput = imgsInput;
dat.coord = [coord_x - mean(coord_x), coord_y - mean(coord_y)];
dat.XInput = imgsInput;
dat.Nsubj = size(imgsInput, 2);

figure; 
for i =1:3
    subplot(1,3, i);
    imagesc(XTrain(:,:,:, i+3)); colormap jet; colorbar;
end

%% simulate 1d curves
% 
% close all;
% coords = ((1:100)' - mean(1:100))/10;
% shifts = [0, 2, -2];
% scalings = [1, 1.2, 0.8];
% betas = [1, 0.8, 1.2];
% Xinputs = [];
% 
% rng default;
% for i = 1:3
% %     dat = [dat cos((coords + shifts(i)) * scalings(i)) + 0.1*randn(length(coords), 1)];
%       Xinputs = [Xinputs betas(i)*exp(-((coords + shifts(i)) * scalings(i)).^2)+ 0.05*randn(length(coords), 1)];
% end
% dat.XInput = Xinputs;
% dat.coord = [coords,ones(length(coords),1)];

%% %
% % addpath("D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\code\Matlab\MCMC\2DSupportFunctions\");
% % 
% % X = XTrain(:,:,:,2); figure; imagesc(X); colormap jet; colorbar;
% % 
% % rotations 
% % thetas = [0, 30, -30, 0];
% % scalings = [1 1; 1.2, 1.2; 0.8, 0.8; 1,1];
% % translations = [0,0; 3,3; -3,-3; 0,0]; % fitted value are opposite
% % Imgs = [];
% % imgsInput = [];
% % for i = 1:3
% %     tempIMG = createAffineTransformedImages(scalings(i,:), thetas(i), translations(i, :), X) + ...
% %         randn(size(X)) * 5e-2;
% %     imgsInput(:,i) = tempIMG(:);
% %     Imgs(:,:,i) = tempIMG;
% % end
% % 
% % [xcoords, ycoords] = find(~isnan(X));
% % coords = [xcoords-mean(xcoords), ycoords-mean(xcoords)];
% % 
% % dat = struct();
% % dat.XInput = imgsInput;
% % dat.coord = coords;
% % dat.Nsubj = size(imgsInput, 2);
% % 
% % figure; 
% % for i =1:3
% %     subplot(1,3, i);
% %     imagesc(Imgs(:,:,i)); colormap jet; colorbar;
% % end

%%
% pwd("D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\code\Matlab");
addpath("D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\code\Matlab\MCMC\supportFunctions\");
addpath("D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\code\Matlab\MCMC\2DSupportFunctions\");

init.b = [1,1,1]';
init.rot = [0, 0, 0];
init.shift = [0,0; 0,0; 0,0];
init.scaling = [1 1; 1,1; 1,1];
init.X = dat.XInput(:,1);
init.sigma = [1e-5, 1e-5, 1e-5]';
init.alpha = .05;
init.rho = 1;
bnds = [0.5,1.5;-45,45;-5,5;.7, 1.3;0,1;0.01,5]; %scaling, rotation, shifting, b, alpha, rho

% dataInput
dataIn.Y = dat.XInput;
dataIn.coord = dat.coord;
dataIn.nsubj = 3;
dataIn.nparts = 3;
dataIn.K1 = 10;
dataIn.K2 = 10;
dataIn.K = 10;
dataIn.nnX = nnsearch_2d(dataIn.coord , dataIn.K);
dataIn.nnDict = nnDictionary_2d(dataIn.coord, 10, dataIn.K);
dataIn.bnds = bnds;
niter = 4e4;
% out2d = MCMC2d(init, dataIn, niter);
out2d = MCMC2d_marginal(init, dataIn, niter);

save D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\results\marginal_2d_fit_digit_7.mat out2d;
%%
init.sigma = [1e-3, 1e-3, 1e-3];
init.b = [1,1,1];
out2d = MCMC2d(init, dataIn, niter);
save D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\results\2d_fit.mat out2d;
%% trace plot
nburnin = 1;
nend = 4e4; 
figure('position', [100,100,1000, 1000]);
% out2d=out;
for i = 1:3

    subplot(8, 3, i);
    plot(nburnin:nend, out2d.Transf(nburnin:nend, i, 1));
    title(['Scaling X', num2str(i)]);
    
    subplot(8, 3, i+3);
    plot(nburnin:nend, out2d.Transf(nburnin:nend, i, 2));
    title(['Scaling Y', num2str(i)]);    
    
    subplot(8, 3, i+6);
    plot(nburnin:nend, out2d.Transf(nburnin:nend, i, 4));
    title(['Shift X', num2str(i)]);
    
    subplot(8, 3, i+9);
    plot(nburnin:nend, out2d.Transf(nburnin:nend, i, 5));
    title(['Shift Y', num2str(i)]);    
    
    subplot(8, 3, i+12);
    plot(nburnin:nend, out2d.Transf(nburnin:nend, i, 3));
    title(['Rotation', num2str(i)]);
    
    subplot(8, 3, i+15);
    plot(nburnin:nend,out2d.b(nburnin:nend, i));
    title(['Beta', num2str(i)]);
    
    subplot(8, 3, i+18);
    plot(nburnin:nend,out2d.sigma(nburnin:nend, i));
    title(['Sigma', num2str(i)]);
    
end

subplot(8, 3, 22);
plot(nburnin:nend,out2d.alpha(nburnin:nend));
title("Alpha");
subplot(8, 3, 23);
plot(nburnin:nend,out2d.rho(nburnin:nend));
title("Rho");
saveas(gcf, "D:\Dropbox\projects\SpatialProjectGIT\SpatialAnalysis\GroupwiseReg\results\plots\simulation_marginal_mcmc_matlab_2d_traceplot_digit_7.png");

%%
LatentImg = Imgs(:,:,1);
LatentImg(:) = mean(out2d.X(nburnin:end, :));
figure; imagesc(LatentImg); colormap jet; colorbar;

%%

