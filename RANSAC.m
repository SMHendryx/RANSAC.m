%Copyright 2018 Sean M. Hendryx
%
%Licensed under the Apache License, Version 2.0 (the "License");
%you may not use this file except in compliance with the License.
%You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS,
%WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%See the License for the specific language governing permissions and
%limitations under the License.

% RANSAC implemenation in Matlab for returning hypothesized model and inliers.

function m = homogenousLS(U)
    % Given the design matrix U, returns the parameters of the homogenous
    % least squares solution IN STANDARD FORM: ax + by = d
    % ONES IN THE LAST COLUMN
    Y = transpose(U) * U;
    [V,D] = eig(Y);
    %produces a diagonal matrix D of eigenvalues and a full matrix V whose columns are the corresponding eigenvectors 
    m = V(:,1);
end

function m = fitFunc(X)
    m = homogenousLS(X);
end

function [m, inliers] = RANSAC(X, fitFunc, k, d, inlierRatio)
    % Returns the best fit straight line parameters to the data in X using RANSAC
    % the minimum number of columns is inferred from the data X
    % param X: the data (with a column of ones in the last column) in
    % an nxdims matrix 
    % k: maximum number of iterations
    % d: threshold distance between the fit model and the points being
    % tested for being a potential inlier
    % inlierRatio: the expected ratio of inliers in the total dataset
    % which determines if we have found a good model
    
    % get dims (which is the mininmum number of points that will be used to describe model, m):
    dims = size(X, 2) - 1;
    % get num points:
    n = size(X, 1);
    % instantiate best-fit model:
    m = zeros(n, 1);
    %instantiate error:
    smallestError = Inf;
    numPointsForCandidateModel = n * inlierRatio;
    
    i = 1;
    while i <= k
        % Get random sample of dims points and those points not selected:
        [maybeInliers, notInMaybeInliers] = getRanSampleOfDimsPoints(X, dims, n);
        distances = getDistancesFromPointsToLine(maybeInliers(1,:), maybeInliers(2,:), notInMaybeInliers);
        
        %get inliers:
        inliersOfMaybeModel = getInliers(notInMaybeInliers, distances, d);
        inliersOfMaybeModel = cat(1, inliersOfMaybeModel, maybeInliers);
        
        % do we have enough points to be a candidate model?
        if size(inliersOfMaybeModel, 1) >= numPointsForCandidateModel
            betterModel = feval(fitFunc, inliersOfMaybeModel);
            error = getRMSError(inliersOfMaybeModel, betterModel);
            if error < smallestError
                m = betterModel;
                inliers = inliersOfMaybeModel;
                smallestError = error;
            end
        end
        
        
        i = i + 1;
    end
    
    
end

function siForm = convertToSlopeInterceptForm(StandardForm)
    % convert coefficients to form y = mx + b
    a = StandardForm(1);
    b = StandardForm(2);
    c = StandardForm(3);
    disp('Slope-Intercept form:')
    m = (-1 * a)/b
    b_1 = (-1 *c)/b
    siForm = [m; b_1];
end

function RMSError = getRMSError(points, model)
    % Not generalized to arbitrary number of dimensions.
    a = model(1);
    b = model(2);
    c = model(3);
    
    orthogErrors = zeros(size(points, 1), 1);
    for i  = 1:size(points, 1)
        d = abs(a * points(i, 1) + b*points(2) + c)/sqrt(a^2 + b^2);
        orthogErrors(i) = d;
    end
    RMSError = rms(orthogErrors);
end

function plot2dFit(X, m, nameOfLine)
    % X is the input 2d design matrix of observations:
    % y in first column, ones in last column
    % explanatory variable in second, and observed y in third.
    % m are the model parameters, last element of m is the intercept
    yHat = X(:,2:end)*m;
    x = X(:,2);
    y = X(:,1);
    figure();
    scatter(x, y);
    hold on
    plot(x,yHat);
    legend('Data',nameOfLine);
    hold off
end

