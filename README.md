# RANSAC
This repo contains a Matlab implementation of RANSAC and associated functions including homogenous least squares for fitting RANSAC and minimizing error in all dimensions. Useful for edge finding in imagery and other computer vision problems.

To run using homogeneous LS:
```f = @fitFunc;
[m, inliers] = RANSAC(X, f, k, d, w);
```
Where `X` is the design matrix of your data with a column of ones in the last column
To plot result:
```
% Plot that line on top of the scatter plot
siForm = convertToSlopeInterceptForm(m)
% RMSE from RANSAC best line:
RMSE = getRMSError(inliers, m)
plot2dFit(X, siForm, 'RANSAC homogenous LS fit');
```

## Disclaimer
`getRMSError` and plotting only functional on 2d data.