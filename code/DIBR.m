function output_image = DIBR(I,disparity_map,p,Np,down_ratio,varargin)
% This program performs the popular Depth Image Based Rendering to generate
% Stereoscopic views for 3D Video Displays
%
% This program accepts the color video and the depth video in the YUV 4:2:0
% format and convert it to the stereoscopic views, known as Left and Right
% Videos. The output Left and Right will be in the same YUV format. rows is
% the height of the video and cols is the width of the video, while numf is
% the number of frames that you need to convert
%
% This program requires the yuv toolbox which can be downloaded from
% http://sprljan.com/nikola/matlab/yuv.html for free. YUV Toolbox is a
% provided by Nikola Sprljan. 
%
% This program implements the rendering framework that is used in the paper
%
%       " D. De Silva, W. Fernando, and H. Kodikaraarachchi, 
%	  A NEW MODE SELECTION TECHNIQUE FOR CODING DEPTH MAPS OF 3D VIDEO,? 
%	   IEEE International Conference on Acoustics, Speech and Signal 
%	   Processing, ICASSP 2010. pp. 686-689.  Mar. 2010."
%
%       which has the following features:
%   
%       1. Depth maps are interpreted according to the MPEG Informative
%          Recommendations in MPEG Doc. N8038.
%       2. Left image is rendered from Right most pixel to the Left most 
%          pixel and Right Image is rendered vice versa. This is done to
%          make sure that no background pixels would appear as foreground.
%       3. Disocclusions are filled with Background pixel extrapolation,
%          however with some small modifications. Disocclusions
%          are filled in the opposite direction to rendering. 
%
% Varuna De Silva, I-Lab, University of Surrey, Guildford, GU2 7XH, UK.
% Version 1.0: May 2010. 
% Please send your comments and additional requirements to the author at
% varunax@gmail.com

%% Input parser
P = inputParser;
% Plot oder nicht
P.addOptional('do_plot', false, @islogical); 
% den Input lesen
P.parse(varargin{:});
do_plot = P.Results.do_plot;  

% linear
if p <=0.5
    xb = 0.5*p;
else
    xb = -0.5*p+0.5;
end
%non_linear
%xb = -(p-0.5)^2+0.25;
k = 10;
xb = k*xb;
disparity_map = round( (disparity_map - min(disparity_map(:)) )/down_ratio);

% disp(min(disparity_map(:)))
% 
% disp(max(disparity_map(:)))

Color = double(I);
if rem(size(Color,1),2)
    Color = [Color;Color(end,:,:)];
    disparity_map = [disparity_map;disparity_map(end,:,:)];
end
if rem(size(Color,2),2)
    Color = [Color,Color(:,end,:)];
    disparity_map = [disparity_map,disparity_map(:,end,:)];
end
[L, R] = cd2lr2(Color,disparity_map,xb,Np);
l = yuv2rgb(L.luma,L.chroma1,L.chroma2);
r = yuv2rgb(R.luma,R.chroma1,R.chroma2);
if p<=0.5
    output_image = r;
else
    output_image = l;
end
if do_plot
%     figure;imshow(output_image);title(['Virtual view with p equal to ',num2str(p)])
%     figure('Name','DIBR left image'),imshow(l);
%     figure('Name','DIBR right image'),imshow(r);
%      figure('Name','DIBR stereo image'),imshow(stereoAnaglyph(l,r))
end
end

%% generate the new image
function [L, R] = cd2lr2(color, depth, xb,Np)
%% set the parameters
nkfar = max(depth(:));
N = nkfar+1; %Number of depth planes
rows = size(color,1);
cols = size(color,2);

%% initialization
[Y, U, V]  = rgb2yuv(color);
U = imresize(U,0.5);
V = imresize(V,0.5);
% figure, imshow(yuv2rgb(Y,U,V));
L.luma = zeros(rows,cols);
L.chroma1 = zeros(rows/2, cols/2);
L.chroma2 = zeros(rows/2, cols/2);
R.luma = zeros(rows,cols);
R.chroma1 = zeros(rows/2, cols/2);
R.chroma2 = zeros(rows/2, cols/2);
shiftx = zeros(1,N);

%% Calculate the shift table to lookup
for i = 1:N
    shiftx(1,i) = find_shiftMC3(i-1,xb,nkfar,Np);
end
% This value is the half of the maximum shift
% Maximum shift comes at depth = 0;
S = max(shiftx(:));
if rem(S,2)
    S = S+1;
end

%% Left Image    
    mask.luma = zeros(rows,cols);
    mask.chroma1 = zeros(rows/2, cols/2);

    % Luma Components
    for i = 1:rows
        %Note the order of rendering
        for j = cols:-1:1
            shift = shiftx(1,(depth(i,j)+1));
            if(j-shift+S<cols)
                L.luma(i,j+S-shift)= Y(i,j);
                mask.luma(i,j+S-shift) = 1;
            end
        end        
    end   
    NHOOD = [1 1 1 0 0 0];
    se = strel('arbitrary', NHOOD);
    y = imdilate(mask.luma,se);
    NHOOD2 = [0 0 0 1 1 1 1 1 1];
    se2 = strel('arbitrary', NHOOD2);
    x = imerode(y, se2);
    z = x.*mask.luma;
    mask.luma = z;   
    for i = 1:rows
        % Filling of disocclusions with Background Pixel Extrapolation
        for j = 1:cols
            if(mask.luma(i,j)==0)
                if(j-7<1)
                    L.luma(i,j) = Y(i,j);
                else
                    L.luma(i,j) = sum(L.luma(i,j-7:j-4))/4;
                end
            end
        end
    end
    
    for i = 3:rows-2
        % Smooth the filled occlusions
        for j = 3:cols-2
            if(mask.luma(i,j)==0)
                L.luma(i,j) = sum(sum(L.luma(i-2:i+2,j-2:j+2)))/25;
            end
        end
    end
    
    % Chroma Components
    for i = 1:rows/2
        for j = cols/2:-1:1
            shift = shiftx(1,(depth(i*2,j*2)+1)); %P.S. It is not cols/2 here

            if(j-ceil(shift/2)+S/2<cols/2)
                L.chroma1(i,j-ceil(shift/2)+S/2)= U(i,j);
                L.chroma2(i,j-ceil(shift/2)+S/2)= V(i,j);
                mask.chroma1(i,j-ceil(shift/2)+S/2) = 1;
            end
        end        
    end

    for i = 1:rows/2
        % Filling of disocclusions with Background Pixel Extrapolation
        for j = 1:cols/2
            if(mask.chroma1(i,j)==0)
                if(j-2<1)
                    L.chroma1(i,j) = U(i,j);
                    L.chroma2(i,j) = V(i,j);
                else
                    L.chroma1(i,j) = sum(L.chroma1(i,j-2:j-1))/2;
                    L.chroma2(i,j) = sum(L.chroma2(i,j-2:j-1))/2;
                end
            end
        end
    end  
    
%% Right Image
    mask.luma = zeros(rows,cols);
    mask.chroma1 = zeros(rows/2, cols/2);   
    
    % Right Image
    for i = 1:rows
        for j = 1:cols
            shift = shiftx(1,(depth(i,j)+1));
            
            if(j+shift-S>1)
                R.luma(i,j+shift-S)= Y(i,j);
                mask.luma(i,j+shift-S) = 1;
            end
        end
    end
    
    NHOOD = [0 0 0 1 1 1];
    se = strel('arbitrary', NHOOD);
    y = imdilate(mask.luma,se);
    NHOOD2 = [1 1 1 1 1 1 0 0 0];
    se2 = strel('arbitrary', NHOOD2);
    x = imerode(y, se2);
    z = x.*mask.luma;
    mask.luma = z;
    
    for i = 1:rows
        % Filling of disocclusions with Background Pixel Extrapolation
        for j=cols:-1:1
            if(mask.luma(i,j) == 0)
                if(j+7>cols)
                    R.luma(i,j) = Y(i,j);
                else
                    R.luma(i,j) = sum(R.luma(i,j+4:j+7))/4;
                end
            end
        end
    end
    
    for i = 3:rows-2
        % Smooth out the filled occlusions
        for j=cols-2:-1:3
            if(mask.luma(i,j) == 0)
                R.luma(i,j) = sum(sum(R.luma(i-2:i+2,j-2:j+2)))/25;
            end
        end
    end
    
    % Chroma Components
    for i = 1:rows/2
        for j = 1:cols/2
            shift = shiftx(1,(depth(i*2,j*2)+1));
            
            if(j+ceil(shift/2)-S/2>1)
                R.chroma1(i,j+ceil(shift/2)-S/2)= U(i,j);
                R.chroma2(i,j+ceil(shift/2)-S/2)= V(i,j);
                
                mask.chroma1(i,j+ceil(shift/2)-S/2) = 1;
            end
        end       
    end
    
    for i = 1:rows/2
        % Filling of disocclusions with Background Pixel Extrapolation
        for j=cols/2:-1:1
            if(mask.chroma1(i,j) == 0)
                if(j+2>cols/2)
                    R.chroma1(i,j) = U(i,j);
                    R.chroma2(i,j) = V(i,j);
                else
                    R.chroma1(i,j) = sum(R.chroma1(i,j+1:j+2))/2;
                    R.chroma2(i,j) = sum(R.chroma2(i,j+1:j+2))/2;
                end
            end
        end
    end
end

%% find the shift for each depth plane
function h = find_shiftMC3(depth,xb,nkfar,Np)	
	N=nkfar+1; % Number of depth planes
	D = 300; %Viewing Distance usually 300 cm 
	kfar = nkfar/16;
    knear = 0;    
	A= kfar - depth*(knear+kfar)/(N-1);
	h= round(xb*Np*A/(D*2));
	if (h<0)
        %It will never come here due to Assumption 1
		h=0-h; 
    end
end

%% Converts RGB to YUV
function [Y,U,V]=rgb2yuv(frame)
%Version: 3.00, Date: 2007/11/21, author: Nikola Sprljan
% References: 
%  [1] Rec. ITU-R BT.601-6
%  [2] Rec. ITU-R BT.709-5
%  [3] http://en.wikipedia.org/wiki/YCbCr
%  [4] http://www.poynton.com/ColorFAQ.html
%  [5] http://www.mathworks.com/access/helpdesk/help/toolbox/vipblks/ref/colorspaceconversion.html
%  [6] Keith Jack, Video Demystified, Chapter 3, http://www.compression.ru/download/articles/color_space/ch03.pdf
%  [7] http://msdn.microsoft.com/library/default.asp?url=/library/en-us/wceddraw/html/_dxce_converting_between_yuv_and_rgb.asp
T = [0.2126, 0.7152,  0.0722;
    -0.1146, -0.3854, 0.5000;
     0.5000, -0.4542, -0.0468];
yuvoffset = [0,128,128];
R = frame(:,:,1);
G = frame(:,:,2);
B = frame(:,:,3);
R = double(R);
G = double(G);
B = double(B);
Y = T(1,1) * R + T(1,2) * G + T(1,3) * B + yuvoffset(1);
U = T(2,1) * R + T(2,2) * G + T(2,3) * B + yuvoffset(2);
V = T(3,1) * R + T(3,2) * G + T(3,3) * B + yuvoffset(3);
Y = uint8(round(Y));
U = uint8(round(U));
V = uint8(round(V));
end

%% Converts YUV to RGB
function rgb=yuv2rgb(Y,U,V)
%Version: 4.00, Date: 2007/11/18, author: Nikola Sprljan

%create the 3D YUV array
yuv = zeros(size(Y,1),size(Y,2),3);
yuv(:,:,1) = double(Y);
yuv(:,:,2) = imresize(double(U),2,'bicubic');
yuv(:,:,3) = imresize(double(V),2,'bicubic');

%inversion of the transform matrix
rgb2yuvT = [0.2126, 0.7152,  0.0722;
           -0.1146, -0.3854, 0.5000;
            0.5000, -0.4542, -0.0468];
T = inv(rgb2yuvT);
yuvoffset = [0,128,128];
rgb = zeros(size(Y,1),size(Y,2),3);

if (yuvoffset(1) ~= 0)
  yuv(:,:,1) = yuv(:,:,1) - yuvoffset(1);
end
if (yuvoffset(2) ~= 0)
  yuv(:,:,2) = yuv(:,:,2) - yuvoffset(2);
end
if (yuvoffset(3) ~= 0)
  yuv(:,:,3) = yuv(:,:,3) - yuvoffset(3);
end
rgb(:,:,1) = T(1,1) * yuv(:,:,1) + T(1,2) * yuv(:,:,2) + T(1,3) * yuv(:,:,3);
rgb(:,:,2) = T(2,1) * yuv(:,:,1) + T(2,2) * yuv(:,:,2) + T(2,3) * yuv(:,:,3);
rgb(:,:,3) = T(3,1) * yuv(:,:,1) + T(3,2) * yuv(:,:,2) + T(3,3) * yuv(:,:,3);
rgb = uint8(round(rgb));
end

