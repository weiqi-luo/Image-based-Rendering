% Detect lines in a grayscale image
%  BUG REPORT:
%  This is an alpha version. Please send your bug reports, comments and
%  suggestions to pengtao@glue.umd.edu . Thanks.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: alpha       Revision: Dec. 05, 2005


%  Reserved for extension
%  accumres:	(Optional, default is [4, 1], minimum is [2, 0.5])
%               The desired resolutions in 'rho' and 'theta'. Is a
%               2-element vector in following format:
%               [resolution in 'rho' (pixels),
%               resolution in 'theta' (degrees)]

function [lines1,lines2,wall_desparity] = line_segmentation(I1,I2,choose_img,varargin)
    %% Input parser
    P = inputParser;
    P.addOptional('do_plot', true, @islogical);
    P.parse(varargin{:});
    do_plot = P.Results.do_plot;
    
    %% prepare the image
    I1 = rgb_to_gray(I1);
    I2 = rgb_to_gray(I2);
    fltr4img = [1 2 3 2 1; 2 3 4 3 2; 3 4 6 4 3; 2 3 4 3 2; 1 2 3 2 1];
    fltr4img = fltr4img / sum(fltr4img(:));
    imgfltrd1 = filter2( fltr4img , I1 );
    imgfltrd2 = filter2( fltr4img , I2 );
    
    if choose_img
        % L1
        grdthres = 4;
        detsens = 0.57;
        x1 = 1500;
        y1 = 1200;
        x2 = 1500;
        y2 = 1200;
        [~, ~, ~, ~, lineseg] = Hough_Grd(imgfltrd1, grdthres, detsens);
        [lines1,I1] = line_find_1(I1, x1,y1,x2,y2); 
        % R1
        grdthres = 1.2;
        detsens = 0.3;
        x1 = 744;
        y1 = 1151;
        x2 = 2200;
        y2 = 1158;
        [~, ~, ~, ~, lineseg] = Hough_Grd(imgfltrd2, grdthres, detsens);
        [lines2,I2] = line_find_1(I2, x1,y1,x2,y2);        
    else
        % L2
        grdthres = 3.5;
        detsens = 0.35;
        x1 = 2141;
        y1 = 887;
        [~, ~, ~, ~, lineseg] = Hough_Grd(imgfltrd1, grdthres, detsens);
        [lines1,I1] = line_find_2(I1, x1,y1);
         % R2
        grdthres = 3;
        detsens = 0.02;
        x1 = 2594;
        y1 = 890;
        [~, ~, ~, ~, lineseg] = Hough_Grd(imgfltrd2, grdthres, detsens);
        [lines2,I2] = line_find_2(I2, x1,y1);
    end    
           
    %% plot
    if do_plot
        figure('Name','line segmentation i1'); imagesc(I1); colormap('gray'); axis image;
        DrawLines_2Ends(lineseg);
        figure('Name','line segmentation i2'); imagesc(I2); colormap('gray'); axis image;
        DrawLines_2Ends(lineseg);
%         figure,imshow(I);
%         figure,imshow(lines);
    end
    
    %% find the lines
    function [lines,I] = line_find_1(I,x1,y1,x2,y2)
        l = zeros(1,4);
        x = lineseg(:,1);
        y = lineseg(:,3);
        lines = zeros(size(I));
        % line 1
        lineseg_ = lineseg(x<x1 & y<y1,:);
        [~,ind] = sort(lineseg_(:,3));
        lineseg_ = lineseg_(ind,:);
        l(1) = lineseg_(1,1);
        l(2) = lineseg_(end,2);     
        l(3) = lineseg_(1,3);
        l(4) = round(lineseg_(end,4));

        line = zeros(2,l(4));
        line(1,:) = 1:l(4);
        line(2,:) = (l(2)-l(1))/(l(4)-l(3))*(line(1,:)-l(3))+l(1);
        line=round(line);

        for ii=-20:20
        ind = sub2ind(size(I),line(1,:),line(2,:)+ii);
        lines(ind) = 1;
        I(ind) = 255;
        end

        % line 2
        lineseg_ = lineseg(x>x2 & y<y2,:);
        [~,ind] = sort(lineseg_(:,3));
        lineseg_ = lineseg_(ind,:);
        l(1) = lineseg_(1,1);
        l(2) = lineseg_(end,2);     
        l(3) = lineseg_(1,3);
        l(4) = round(lineseg_(end,4));

        line = zeros(2,l(4));
        line(1,:) = 1:l(4);
        line(2,:) = (l(2)-l(1))/(l(4)-l(3))*(line(1,:)-l(3))+l(1);
        line=round(line);

        for ii=-20:20
        ind = sub2ind(size(I),line(1,:),line(2,:)+ii);
        lines(ind) = 1;
        I(ind) = 255;
        end   
        lines = logical(lines); 
    end


    function [lines,I] = line_find_2(I, x1,y1)
        l = zeros(1,4);
        x = lineseg(:,1);
        y = lineseg(:,3);
        lines = zeros(size(I)); 
        
        lineseg_ = lineseg(x>x1 & y<y1 ,:);
        [~,ind] = sort(lineseg_(:,3));
        lineseg_ = lineseg_(ind,:);
        l(1) = lineseg_(1,1);
        l(2) = lineseg_(end,2);     
        l(3) = lineseg_(1,3);
        l(4) = round(lineseg_(end,4));

        line = zeros(2,l(4));
        line(1,:) = 1:l(4);
        line(2,:) = (l(2)-l(1))/(l(4)-l(3))*(line(1,:)-l(3))+l(1);
        line=round(line);

        for ii=-20:20
        ind = sub2ind(size(I),line(1,:),line(2,:)+ii);
        lines(ind) = 1;
        I(ind) = 255;
        end
        lines = logical(lines); 
    end
    
end
 


%% semgment the line
function [accum, axis_rho, axis_theta, varargout] = ...
    Hough_Grd(img, varargin)

% Validation of arguments
if ndims(img) ~= 2 || ~isnumeric(img),
    error('Hough_Grd: ''img'' has to be 2 dimensional');
end
if ~all(size(img) >= 16),
    error('Hough_Grd: ''img'' has to be larger than 16-by-16');
end

% Parameters (default values)
prm_grdthres = 8;
prm_accumres = [4, 1];
prm_detsens = 0.08;

func_compu_lineprm = true;
prm_fltraccum = true;

% Validation of arguments
vap_grdthres = 1;
if nargin > vap_grdthres,
    if isnumeric(varargin{vap_grdthres}) && ...
            varargin{vap_grdthres}(1) >= 0,
        prm_grdthres = varargin{vap_grdthres}(1);
    else
        error(['Hough_Grd: ''grdthres'' has to be ', ...
            'a non-negative number']);
    end
end
%{
vap_accumres = 3;
if nargin > vap_accumres,
    if numel(varargin{vap_accumres}) == 2 && ...
            isnumeric(varargin{vap_accumres}) && ...
            ( varargin{vap_accumres}(1) >= 2 && ...
            varargin{vap_accumres}(2) >= 0.5 ),
        prm_accumres = varargin{vap_accumres};
    else
        error(['Hough_Grd: ''accumres'' has to be a two-element ', ...
            'vector and no smaller than [2, 0.5]']);
    end
end
%}
vap_detsens = 2;
if nargin > vap_detsens
    if isnumeric(varargin{vap_detsens}) && ...
            varargin{vap_detsens}(1) > 0 && varargin{vap_detsens}(1) < 1,
        prm_detsens = varargin{vap_detsens};
    else
        error('Hough_Grd: ''detsens'' has to be in the range (0, 1)');
    end
end

func_compu_lineprm = ( nargout > 3 );

% Size of the accumulation array
imgsize = size(img);
coef_rhorng = [ -imgsize(2), sqrt(sum(imgsize.^2)) ];
coef_thetarng = [-pi/18, pi+pi/18];

prm_accumsize = [ ...
    round( (coef_rhorng(2)-coef_rhorng(1)) * (2/prm_accumres(1)) ) , ...
    round( (coef_thetarng(2)-coef_thetarng(1)) * (180/pi) * ...
    (2/prm_accumres(2)) ) ];

% Default filter for the accumulation array
prm_acmfltr_R = 4;
prm_acmfltr_w = [1 2 4 8];

fltr4accum = ones(2 * prm_acmfltr_R - 1) * prm_acmfltr_w(1);
for k = 2 : prm_acmfltr_R,
    fltr4accum(k:(2*prm_acmfltr_R-k), k:(2*prm_acmfltr_R-k)) = ...
        prm_acmfltr_w(k);
end
fltr4accum = fltr4accum / sum(fltr4accum(:));

% Parameters for the algorithm using repeated
% thresholding and segmentation
prm_lp_dthres = min([ 0.1, (0.1+prm_detsens)/2 ]);
prm_lp_maxthres = 0.8;

% Parameters for the algorithm using local maximum filter
prm_useaoi = false;
prm_aoiminsize = [8, 8];

prm_fltrLM_R = 4;
prm_fltrLM_s = 1.3;
prm_fltrLM_r = ceil( prm_fltrLM_R * 0.6 );
prm_fltrLM_npix = 6;

% Reserved parameters
dbg_on = false;      % debug information
dbg_bfigno = 5;
if nargout > 5,  dbg_on = true;  end


%%%%%%%% Building accumulation array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the gradient and the magnitude of gradient
img = double(img);
[grdx, grdy] = gradient(img);
grdmag = sqrt(grdx.^2 + grdy.^2);

% Clear the margins of the gradient field
prm_grdfldmgn = 4;
grdmag([1:prm_grdfldmgn, (end-prm_grdfldmgn+1):end], :) = 0;
grdmag(:, [1:prm_grdfldmgn, (end-prm_grdfldmgn+1):end]) = 0;

% Get the linear indices, as well as the subscripts, of the pixels
% whose gradient magnitudes are larger than the given threshold
grdmasklin = find(grdmag > prm_grdthres);
[grdmask_IdxI, grdmask_IdxJ] = ind2sub(imgsize, grdmasklin);

% Compute (the line parameter) 'theta' for all voting pixels
grdphs_vot = atan2( grdy(grdmasklin), grdx(grdmasklin) );

% -- Regulate 'grdphs_vot' from [-pi, pi] to [0, pi]
grdphs_vot = grdphs_vot + pi * (grdphs_vot < 0);

% Compute the 'theta'-subscript (to the accumulation array)
% for all voting pixels
coef_subtheta = prm_accumsize(2) / ...
    (coef_thetarng(2) - coef_thetarng(1)) * (1 - 1e-6);
sub_theta = ceil( (grdphs_vot - coef_thetarng(1)) * coef_subtheta );

% 'theta' vector for the accumulation array
axis_theta = (coef_thetarng(1) + 0.5 / coef_subtheta) : ...
    (1 / coef_subtheta) : coef_thetarng(2);

% Compute the 'rho' values for all voting pixels
% rho = (J - 0.5) * cos(theta) + (I - 0.5) * sin(theta)
rho_vot = (grdmask_IdxJ - 0.5) .* cos(grdphs_vot) + ...
    (grdmask_IdxI - 0.5) .* sin(grdphs_vot);

% Compute the 'rho'-subscript (to the accumulation array)
% for all voting pixels
coef_subrho = prm_accumsize(1) / ...
    (coef_rhorng(2) - coef_rhorng(1)) * (1 - 1e-6);
sub_rho = ceil( (rho_vot - coef_rhorng(1)) * coef_subrho );

% 'rho' vector for the accumulation array
axis_rho = (coef_rhorng(1) + 0.5 / coef_subrho) : ...
    (1 / coef_subrho) : coef_rhorng(2);

% Build the accumulation array, using gradient magnitude as weight
accum = accumarray( sub2ind(prm_accumsize, sub_rho, sub_theta), ...
    grdmag(grdmasklin) );
accum = [ accum ; ...
    zeros(prm_accumsize(1) * prm_accumsize(2) - numel(accum), 1) ];
accum = reshape( accum, prm_accumsize );


%%%%%%%% Locating peaks in the accumulation array %%%%%%%%%%%%%%%%%%%

% Stop if no need to estimate the parameters of the lines
if ~func_compu_lineprm,
    return;
end

% Smooth the accumulation array
if prm_fltraccum,
    accum = filter2( fltr4accum, accum );
end

% Find the maximum value in the accumulation array
accum_max = max(accum(:));


%------- Algorithm 1: Local maximum filter (begin) -------------
% Build the local maximum filter
fltr4LM = zeros(2 * prm_fltrLM_R + 1);

[mesh4fLM_x, mesh4fLM_y] = meshgrid(-prm_fltrLM_R : prm_fltrLM_R);
mesh4fLM_r = sqrt( mesh4fLM_x.^2 + mesh4fLM_y.^2 );
fltr4LM_mask = ...
	( mesh4fLM_r > prm_fltrLM_r & mesh4fLM_r <= prm_fltrLM_R );
fltr4LM = fltr4LM - ...
	fltr4LM_mask * (prm_fltrLM_s / sum(fltr4LM_mask(:)));

if prm_fltrLM_R >= 4,
	fltr4LM_mask = ( mesh4fLM_r < (prm_fltrLM_r - 1) );
else
	fltr4LM_mask = ( mesh4fLM_r < prm_fltrLM_r );
end
fltr4LM = fltr4LM + fltr4LM_mask / sum(fltr4LM_mask(:));

% Select a number of Areas-Of-Interest from the accumulation array
if prm_useaoi,
    % Thresholding and segmentation
    accummask = ( accum > (accum_max * prm_detsens) );
    [accumlabel, accum_nRgn] = bwlabel( accummask, 8 );

    % Select AOIs from segmented regions
    accumAOI = ones(0,4);
    for k = 1 : accum_nRgn,
        accumrgn_lin = find( accumlabel == k );
        [accumrgn_IdxI, accumrgn_IdxJ] = ...
            ind2sub( size(accumlabel), accumrgn_lin );
        rgn_top = min( accumrgn_IdxI );
        rgn_bottom = max( accumrgn_IdxI );
        rgn_left = min( accumrgn_IdxJ );
        rgn_right = max( accumrgn_IdxJ );        
        % The AOIs selected must satisfy a minimum size
        if ( (rgn_bottom - rgn_top + 1) >= prm_aoiminsize(1) || ...
                (rgn_right - rgn_left + 1) >= prm_aoiminsize(2) ),
            accumAOI = [ accumAOI; ...
                max([ 1, rgn_top - prm_fltrLM_R ]), ...
                min([ size(accum,1), rgn_bottom + prm_fltrLM_R ]), ...
                max([ 1, rgn_left - prm_fltrLM_R ]), ...
                min([ size(accum,2), rgn_right + prm_fltrLM_R ]) ];
        end
    end
else
    % Whole accumulation array as the one AOI
    accumAOI = [1, size(accum,1), 1, size(accum,2)];
end

% **** Debug code (begin)
if dbg_on && prm_useaoi,
    dbg_accumLM = zeros(size(accum));
end
% **** Debug code (end)

% For each of the AOIs selected, locate the local maxima
lineprm = zeros(0,2);
for k = 1 : size(accumAOI, 1),
    aoi = accumAOI(k,:);    % just for referencing convenience

    % Apply the local maxima filter
    candLM = conv2( accum(aoi(1):aoi(2), aoi(3):aoi(4)) , ...
        fltr4LM , 'same' );

	% Thresholding of 'candLM' & 'accum'
    if prm_useaoi,
        candLM_mask = ( candLM > (max(candLM(:))*prm_detsens) & ...
            accummask(aoi(1):aoi(2),aoi(3):aoi(4)) );
    else
        candLM_mask = ( candLM > (max(candLM(:))*prm_detsens) & ...
            accum(aoi(1):aoi(2),aoi(3):aoi(4)) > (accum_max*prm_detsens) );
    end

    % Clear the margins of 'candLM_mask'
    candLM_mask([1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end], :) = 0;
    candLM_mask(:, [1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end]) = 0;

    % **** Debug code (begin)
    if dbg_on && prm_useaoi,
        dbg_accumLM(aoi(1):aoi(2), aoi(3):aoi(4)) = ...
            dbg_accumLM(aoi(1):aoi(2), aoi(3):aoi(4)) + candLM;
    end
    % **** Debug code (end)

    % Group the local maxima candidates by adjacency, compute the
    % centroid position for each group and take that as the center
    % of one circle detected
    [candLM_label, candLM_nRgn] = bwlabel( candLM_mask, 8 );

    for ilabel = 1 : candLM_nRgn,
        % Indices (to current AOI) of the pixels in the group
        candgrp_masklin = find( candLM_label == ilabel );
        [candgrp_IdxI, candgrp_IdxJ] = ...
            ind2sub( size(candLM_label) , candgrp_masklin );

        % Indices (to 'accum') of the pixels in the group
        candgrp_IdxI = candgrp_IdxI + ( aoi(1) - 1 );
        candgrp_IdxJ = candgrp_IdxJ + ( aoi(3) - 1 );
        candgrp_idx2acm = ...
            sub2ind( size(accum) , candgrp_IdxI , candgrp_IdxJ );

        % Minimum number of qulified pixels in the group
        if sum(numel(candgrp_masklin)) < prm_fltrLM_npix,
            continue;
        end

        % Compute the centroid position
        candgrp_acmsum = sum( accum(candgrp_idx2acm) );
        cc_rho = sum( candgrp_IdxI .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
        cc_theta = sum( candgrp_IdxJ .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
        lineprm = [lineprm; cc_rho, cc_theta];
    end
end

% **** Debug code (begin)
if dbg_on,
    figure(dbg_bfigno);
    if prm_useaoi,
        imagesc(dbg_accumLM); axis xy;
    else
        imagesc(candLM); axis xy;
    end
end
% **** Debug code (end)
%------- Algorithm 1: Local maximum filter (end) ---------------


%------- Algo 2: Repeated thresholding & segmentation (begin) --
%{
% Locate the peaks (in pixel coordinates) in the accumulation array
if dbg_on,
    [lineprm, dbg_label] = RecursSegment( accum , ...
        accum_max * [prm_detsens, prm_lp_dthres, prm_lp_maxthres] );
else
    lineprm = RecursSegment( accum , ...
        accum_max * [prm_detsens, prm_lp_dthres, prm_lp_maxthres] );
end

% **** Debug code (begin)
if dbg_on,
    figure(dbg_bfigno);
    imagesc(dbg_label); axis xy; axis equal;
    hold on;
    plot(lineprm(:,2), lineprm(:,1), ...
        'w+', 'LineWidth', 2, 'MarkerSize', 6);
    hold off;
end
% **** Debug code (end)
%}
%------- Algo 2: Repeated thresholding & segmentation (end) ----


% Convert 'lineprm' from pixel coordinates to (rho, theta)
lineprm = [ (lineprm(:,1) - 0.5)/coef_subrho + coef_rhorng(1), ...
	(lineprm(:,2) - 0.5)/coef_subtheta + coef_thetarng(1) ];

% Output 'lineprm'
varargout{1} = lineprm;
if nargout <= 4,
    return;
end


%%%%%%%% Locating the line segments in the raw image %%%%%%%%%%%%%%%%

% Parameters for locating the line segments
prm_ls_tTol = [-pi/7.5, pi/7.5];
prm_ls_pTol = [-3.5, 3.5];
prm_ls_minlen = 10;

% Locate the two ends for all lines detected
rgnmask = logical(zeros(imgsize));
lineseg = zeros(0, 4);

for k = 1 : size(lineprm, 1),
    % Compute the 'rho' values for all pixels voted,
    % assuming they contribute to the current line
    rho2_vot = (grdmask_IdxJ - 0.5) * cos( lineprm(k,2) ) + ...
        (grdmask_IdxI - 0.5) .* sin( lineprm(k,2) );

    % Find the pixels that belong to the line
    PixOnLn_votmask = ( grdphs_vot > (lineprm(k,2) + prm_ls_tTol(1)) & ...
        grdphs_vot < (lineprm(k,2) + prm_ls_tTol(2)) & ...
        rho2_vot > (lineprm(k,1) + prm_ls_pTol(1)) & ...
        rho2_vot < (lineprm(k,1) + prm_ls_pTol(2)) );

	PixOnLn_lin = grdmasklin( find(PixOnLn_votmask) );
    if isempty(PixOnLn_lin),
        continue;
    end
    [PixOnLn_IdxI, PixOnLn_IdxJ] = ind2sub(imgsize, PixOnLn_lin);

    % Find the axis-aligned bounding box for these pixels
    bndbox = [ min(PixOnLn_IdxI), max(PixOnLn_IdxI), ...
        min(PixOnLn_IdxJ), max(PixOnLn_IdxJ) ];

    % Seperate the line segments by adjacencies
    rgnmask( bndbox(1):bndbox(2), bndbox(3):bndbox(4) ) = 0;
    rgnmask( PixOnLn_lin ) = 1;
    [rgnlabel, nrgn] = bwlabel( ...
        rgnmask(bndbox(1):bndbox(2), bndbox(3):bndbox(4)), 8 );

    for k_rgn = 1 : nrgn,
        % Get the linear indices of pixels that belong to one segment
        seg_lin = find( rgnlabel == k_rgn );
        [seg_IdxI, seg_IdxJ] = ind2sub( size(rgnlabel), seg_lin );

        % Ignore if the line segment is too short
        segspan = [ min(seg_IdxI) + bndbox(1) - 1, ...
            max(seg_IdxI) + bndbox(1) - 1, ...
            min(seg_IdxJ) + bndbox(3) - 1, ...
            max(seg_IdxJ) + bndbox(3) - 1 ];
        if (segspan(2) - segspan(1) + 1) < prm_ls_minlen && ...
                (segspan(4) - segspan(3) + 1) < prm_ls_minlen,
            continue;
        end

        % Get the [x1 x2 y1 y2] structure for the line SEGMENT
        if lineprm(k,2) > pi/4 && lineprm(k,2) < 3*pi/4,
            % ls_base_x = 0;
            ls_base_y = lineprm(k,1) / sin(lineprm(k,2));
            lineseg = [ lineseg; segspan(3) - 0.5, segspan(4) - 0.5, ...
                ls_base_y - (segspan(3) - 0.5) / tan(lineprm(k,2)), ...
                ls_base_y - (segspan(4) - 0.5) / tan(lineprm(k,2)) ];
        else
            % ls_base_y = 0;
            ls_base_x = lineprm(k,1) / cos(lineprm(k,2));
            lineseg = [ lineseg; ...
                ls_base_x - (segspan(1) - 0.5) * tan(lineprm(k,2)), ...
                ls_base_x - (segspan(2) - 0.5) * tan(lineprm(k,2)), ...
                segspan(1) - 0.5, segspan(2) - 0.5 ];
        end
    end
end

% Output 'lineseg'
varargout{2} = lineseg;
if nargout > 5,
    varargout{3} = dbg_label;
end
end


%% Recursive thresholding and segmentation
function [lineprm, varargout] = RecursSegment(accum, struct_thres)
% 'struct_thres' contains [thres, deltathres, maxthres]
% 'lineprm' is in pixel coordinate (w.r.t. 'accum')

% Parameters
prm_as_minpixn = 3;
prm_as_maxsize = [12, 12];

% Thresholding
accummask = ( accum > struct_thres(1) );

% Segmentation and locating the centroids of individual regions
[accumlabel, accum_nRgn] = bwlabel( accummask, 8 );

% Segmentation label (for debug purpose)
func_seglbl = ( nargout > 1 );
if func_seglbl,
    seglbl_lblshft = 4;
    seglabel = accumlabel;
end

lineprm = zeros(0, 2);
for k = 1 : accum_nRgn,
    % Linear indices of the pixels in one connected component
	acmrgn_lin = find( accumlabel == k );
    if numel(acmrgn_lin) < prm_as_minpixn,
        continue;
    end
    % Subscripts of the pixels in one connected component
	[acmrgn_IdxPho, acmrgn_IdxTheta] = ...
        ind2sub( size(accumlabel), acmrgn_lin );

	% Further segmentation if the connected region is too big, or
    % computing the centroid of the region
    % -- Axis-aligned bounding box for the region
    bndbox = [ min(acmrgn_IdxPho), max(acmrgn_IdxPho), ...
        min(acmrgn_IdxTheta), max(acmrgn_IdxTheta) ];

    % -- Further segmentation
    bAddCentrdOfRgn = true;
	if ( (bndbox(2) - bndbox(1) + 1) > prm_as_maxsize(1) || ...
            (bndbox(4) - bndbox(3) + 1) > prm_as_maxsize(2) ) && ...
            (struct_thres(1) + struct_thres(2)) <= struct_thres(3),
        if func_seglbl,
            [lineprm_sub, seglabel_sub] = RecursSegment( ...
                accum(bndbox(1):bndbox(2), bndbox(3):bndbox(4)), ...
                [struct_thres(1) + struct_thres(2), struct_thres(2:3)] );
            seglabel(bndbox(1):bndbox(2), bndbox(3):bndbox(4)) = ...
                seglabel(bndbox(1):bndbox(2), bndbox(3):bndbox(4)) + ...
                seglabel_sub + (seglabel_sub > 0) * seglbl_lblshft;
        else
            lineprm_sub = RecursSegment( ...
                accum(bndbox(1):bndbox(2), bndbox(3):bndbox(4)), ...
                [struct_thres(1) + struct_thres(2), struct_thres(2:3)] );
        end
        if ~isempty(lineprm_sub),
            lineprm = [lineprm; lineprm_sub(:,1) + bndbox(1) - 1, ...
                lineprm_sub(:,2) + bndbox(3) - 1 ];
            bAddCentrdOfRgn = false;
        end
    end

    % -- Computing the centroid of the whole region
    if bAddCentrdOfRgn,
        acmrgn_acmsum = sum( accum(acmrgn_lin) );
        lp_IdxPho = sum( acmrgn_IdxPho .* accum(acmrgn_lin) ) / ...
            acmrgn_acmsum;
    	lp_IdxTheta = sum( acmrgn_IdxTheta .* accum(acmrgn_lin) ) / ...
            acmrgn_acmsum;
        lineprm = [ lineprm; lp_IdxPho, lp_IdxTheta ];
	end
end

% Output the segmentation label
if func_seglbl
    varargout{1} = seglabel;
end
end

%% Draw line segments, parameterized as (x1, x2, y1, y2), on graph
function DrawLines_2Ends(lineseg, varargin)
%
%  DrawLines_2Ends(lineseg, varargin)
%  A simple function for drawing line segments on graph. Made as an
%  auxiliary tool for function '[...] = Hough_Grd(...)'.
%
%  INPUT: (lineseg, properties)
%  lineseg:     Parameters (x1, x2, y1, y2) of line segments to draw.
%               Is a Ns-by-4 matrix with each row contains the parameters
%               (x1, x2, y1, y2) that define the two ends of a line
%               segment. The output 'lineseg' from the function
%               '[...] = Hough_Grd(...)' can be put here directly.
%  properties:  (Optional)
%               A string of line drawing properties. Will be transferred
%               to function 'plot' without modification for line drawing.
%
%  OUTPUT: None
%
%  BUG REPORT:
%  Please send your bug reports, comments and suggestions to
%  pengtao@glue.umd.edu . Thanks.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: alpha       Revision: Dec. 02, 2005


hold on;
for k = 1 : size(lineseg, 1),
    % The image origin defined in function '[...] = Hough_Grd(...)' is
    % different from what is defined in Matlab, off by (0.5, 0.5).
    if nargin > 1
        plot(lineseg(k,1:2)+0.5, lineseg(k,3:4)+0.5, varargin{1});
    else
        plot(lineseg(k,1:2)+0.5, lineseg(k,3:4)+0.5, 'LineWidth', 2);
    end
end
hold off;
end
