function [I1Rect,I2Rect] = image_rectification(I1,I2,t1,t2,varargin)
    %% Input parser
    P = inputParser;
    % Plot oder nicht
    P.addOptional('do_plot', false, @islogical); 
    % Plot oder nicht
    P.addOptional('outputview', 'valid'); 
    % den Input lesen
    P.parse(varargin{:});
    outputView = P.Results.outputview;
    do_plot = P.Results.do_plot;   
    
    %% Compute the transformed location of image corners.
    outPts = zeros(8, 2);

    numRows = size(I1, 1);
    numCols = size(I1, 2);
    inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
    [outPts(1:4,1),outPts(1:4,2)] = point_forward(inPts(:,1),inPts(:,2),t1);

    numRows = size(I2, 1);
    numCols = size(I2, 2);
    inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
    [outPts(5:8,1),outPts(5:8,2)] = point_forward(inPts(:,1),inPts(:,2),t2);

    xSort   = sort(outPts(:,1));
    ySort   = sort(outPts(:,2));
    xLim = zeros(1, 2);
    yLim = zeros(1, 2);

    if strcmp(outputView, 'valid')
        % the output images are cropped to the size of the largest common rectangle containing valid pixels.   
        xLim(1) = ceil(xSort(4)) - 0.5;
        xLim(2) = floor(xSort(5)) + 0.5;
        yLim(1) = ceil(ySort(4)) - 0.5;
        yLim(2) = floor(ySort(5)) + 0.5;
        width   = xLim(2) - xLim(1) - 1;
        height  = yLim(2) - yLim(1) - 1;
        R_out = imref2d([height, width], xLim, yLim);

    elseif strcmp(outputView, 'full')
        % the output images include all pixels from the original images.
        xLim(1) = ceil(xSort(1)) - 0.5;
        xLim(2) = floor(xSort(8)) + 0.5;
        yLim(1) = ceil(ySort(1)) - 0.5;
        yLim(2) = floor(ySort(8)) + 0.5;
        width   = xLim(2) - xLim(1) - 1;
        height  = yLim(2) - yLim(1) - 1;
        R_out = imref2d([height, width], xLim, yLim);
    end

    %% Transform the images.
    tform1 = projective2d(t1);
    tform2 = projective2d(t2);
    I1Rect = imwarp(I1, tform1,  'OutputView', R_out);
    I2Rect = imwarp(I2, tform2,  'OutputView', R_out);
    
    %% plot the result
    if do_plot    
        figure('Name','image_rectification');
        if numel(size(I1))<3
            colormap('gray');
            imagesc([I1Rect,I2Rect]);
        else
            imshow([I1Rect,I2Rect]);
        end
    end
end


%% warp the image
function [outputImage] = image_warp(I,t,R_out)    
    %% getSourceMapping;   
    R_in = imref2d(size(I));
    
     % Form plaid grid of intrinsic points in output image.
    [dstXIntrinsic,dstYIntrinsic] = meshgrid(1:R_out.ImageSize(2),1:R_out.ImageSize(1));

    % Define affine transformation that maps from intrinsic system of
    % output image to world system of output image.
    Sx = R_out.PixelExtentInWorldX;
    Sy = R_out.PixelExtentInWorldY;
    Tx = R_out.XWorldLimits(1)-R_out.PixelExtentInWorldX*(R_out.XIntrinsicLimits(1));
    Ty = R_out.YWorldLimits(1)-R_out.PixelExtentInWorldY*(R_out.YIntrinsicLimits(1));
    tIntrinsictoWorldOutput = [Sx 0 0; 0 Sy 0; Tx Ty 1];

    % Define affine transformation that maps from world system of
    % input image to intrinsic system of input image.
    Sx = 1/R_in.PixelExtentInWorldX;
    Sy = 1/R_in.PixelExtentInWorldY;
    Tx = (R_in.XIntrinsicLimits(1))-1/R_in.PixelExtentInWorldX*R_in.XWorldLimits(1);
    Ty = (R_in.YIntrinsicLimits(1))-1/R_in.PixelExtentInWorldY*R_in.YWorldLimits(1);
    tWorldToIntrinsicInput = [Sx 0 0; 0 Sy 0; Tx Ty 1];

    % Form transformation to go from output intrinsic to input intrinsic
    tComp = tIntrinsictoWorldOutput / t * tWorldToIntrinsicInput;

    % Find the transform to go from input intrinsic to output intrinsic
    tformComposite = invert(projective2d(tComp));
    [Xq,Yq] = tformComposite.transformPointsInverse(dstXIntrinsic,dstYIntrinsic);    
    
    %% interpolation
    X = 1:size(I,2); 
    Y = 1:size(I,1);
    I = I.';
    F = griddedInterpolant({X,Y},I,'linear','none');
    
    %% output the image
    outputImage = F(Xq,Yq);
%     outputImage(~isnan(outputImage)) = 256;
end

function [out1,out2] = point_forward(u,v,T)

    T = double(T);
    x = imlincomb(T(1,1),u, T(2,1),v,  T(3,1), class(T));
    y = imlincomb(T(1,2),u, T(2,2),v,  T(3,2), class(T));
    z = imlincomb(T(1,3),u, T(2,3),v,  T(3,3), class(T)); 

    out1 = x ./ z;
    out2 = y ./ z;

end
