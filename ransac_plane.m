function disparity_map = ransac_plane(disparity_map,cluster,maxDistance,varargin)
%% Input parser
P = inputParser;
% Plot oder nicht
P.addOptional('do_plot', false, @islogical); 
% den Input lesen
P.parse(varargin{:});
do_plot = P.Results.do_plot;  

%% transform the cluster into 3D point
ind = find(cluster);
[v,u] = ind2sub(size(cluster),ind);
w = disparity_map(ind);
point3D = [u,v,w];

%%  set the parameter
% Set the maximum point-to-plane distance for plane fitting.
ransacParams.maxDistance = maxDistance;
ransacParams.maxNumTrials = 1000;
ransacParams.confidence = 99;
% Use three points to fit a plane
ransacParams.sampleSize = 3;

%% detect the plane and compute its geometric model parameter
[inliers,bestModelParams] = msac(point3D,ransacParams);

%% extract the plane
plane = point3D(inliers,:);

%% fit the plane to the cluster
z = (-bestModelParams(1)*u-bestModelParams(2)*v-bestModelParams(4))/bestModelParams(3);
disparity_map(ind)=z;

%% plot the result
if do_plot
%     figure('Name','ransac Plane'),pcshow(plane);
end
end

%%
function [inliers,bestModelParams] = msac(allPoints, params, varargin)
confidence = params.confidence;
sampleSize = params.sampleSize;
maxDistance = params.maxDistance;
threshold = cast(maxDistance, 'like', allPoints);
numPts    = size(allPoints,1);
idxTrial  = 1;
numTrials = int32(params.maxNumTrials);
maxDis    = cast(threshold * numPts, 'like', allPoints);
bestDis   = maxDis;
bestModelParams = zeros(0, 'like', allPoints);
maxSkipTrials = params.maxNumTrials * 10;
skipTrials = 0;    
bestInliers = false(numPts, 1);

while idxTrial <= numTrials && skipTrials < maxSkipTrials
    % Random selection without replacement
    indices = randperm(numPts, sampleSize);
    
    % Compute a model from samples
    samplePoints = allPoints(indices, :);
    modelParams = fitFunc(samplePoints, varargin{:});
    
    % Validate the model 
    isValidModel = checkFunc(modelParams, varargin{:});
    
    if isValidModel
        % Evaluate model with truncated loss
        [model, dis, accDis] = evaluateModel( modelParams, ...
            allPoints, threshold, varargin{:});

        % Update the best model found so far
        if accDis < bestDis
            bestDis = accDis;
            bestInliers = dis < threshold;
            bestModelParams = model;
            inlierNum = cast(sum(dis < threshold), 'like', allPoints);
            num = computeLoopNumber(sampleSize,confidence, numPts, inlierNum);
            numTrials = min(numTrials, num);
        end
        
        idxTrial = idxTrial + 1;
    else
        skipTrials = skipTrials + 1;
    end
end
isFound = checkFunc(bestModelParams, varargin{:}) && ~isempty(bestInliers) && sum(bestInliers(:)) >= sampleSize;
if isFound
    if isfield(params, 'recomputeModelFromInliers') && ...
            params.recomputeModelFromInliers
        modelParams = fitFunc(allPoints(bestInliers, :, :), varargin{:});        
        [bestModelParams, dis] = evaluateModel(modelParams, ...
            allPoints, threshold, varargin{:});
        isValidModel = checkFunc(bestModelParams, varargin{:});
        inliers = (dis < threshold);
        if ~isValidModel || ~any(inliers)
            disp('cc')
            isFound = false;
            inliers = false(size(allPoints, 1), 1);
            return;
        end
    else
        inliers = bestInliers;
    end
else
    inliers = false(size(allPoints, 1), 1);
end
end

%% evaluate the Model
function [modelOut, distances, sumDistances] = evaluateModel(modelIn,allPoints, threshold, varargin)
dis = evalFunc(modelIn, allPoints);
dis(dis > threshold) = threshold;
accDis = sum(dis);
if iscell(modelIn)
    [sumDistances, minIdx] = min(accDis);
    distances = dis(:, minIdx);
    modelOut = modelIn{minIdx(1)};
else
    distances = dis;
    modelOut = modelIn;
    sumDistances = accDis;
end
end

%% compute the loop number
function N = computeLoopNumber(sampleSize, confidence, pointNum, inlierNum)
pointNum = cast(pointNum, 'like', inlierNum);
inlierProbability = (inlierNum/pointNum)^sampleSize;

if inlierProbability < eps(class(inlierNum))
    N = intmax('int32');
else
    conf = cast(0.01, 'like', inlierNum) * confidence;
    one  = ones(1,    'like', inlierNum);
    num  = log10(one - conf);
    den  = log10(one - inlierProbability);
    N    = int32(ceil(num/den));
end 
end

%% Plane equation: ax + by + cz + d = 0;
function model = fitFunc(points, varargin)
a = points(2, :) - points(1, :);
b = points(3, :) - points(1, :);
% Cross product
normal = [a(2).*b(3)-a(3).*b(2), ...
          a(3).*b(1)-a(1).*b(3), ...
          a(1).*b(2)-a(2).*b(1)];
denom = sum(normal.^2);
if denom < eps(class(points))
    model = [];
else
    normal = normal / sqrt(denom);
    d = -points(1,:) * normal';
    model = [normal, d];
end
end
 
%% Calculate the distance from the point to the plane.
% D = (a*x + b*y + c*z + d)/sqrt(a^2+b^2+c^2). Denominator is 1 because the normal is a unit vector here.
function dis = evalFunc(model, points)
dis = abs(points * model(1:3)' + model(4));
end

%% Validate the plane coefficients
function isValid = checkFunc(model)
isValid = (numel(model) == 4 & all(isfinite(model)));      
end
