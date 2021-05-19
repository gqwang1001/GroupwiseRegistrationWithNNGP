function loglik = posteriorLog(Y, coordY, coord, bi, ncoord, nsubj, nparts, K1, K2, alpha, rho, sigma,Sx, parworkers)

        coordLong = Reshape2d(coordY);
        [nnidxs, bs, ~] = nns_2d_square_parallel(coord, coordLong,Sx, ncoord, nsubj, K1,alpha, rho, parworkers, false);
        nnidxAPPR = nnsearch_approx(coordLong, K2, ncoord/nparts, nparts, parworkers);
        NNGPout = MarginalNNGPMuS_APPR_v1(Y, nsubj, bi, alpha, sigma, nnidxAPPR, K2, Sx, bs, nnidxs, false, parworkers);
        loglik = NNGPout.loglik;
        
end

