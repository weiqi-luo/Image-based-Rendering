function reconstruction_3D(I1,I2,disparity_map,varargin)
%% Input parser
P = inputParser;
% Plot oder nicht
P.addOptional('do_plot', true, @islogical);
% den Input lesen
P.parse(varargin{:});
do_plot = P.Results.do_plot;

D = disparity_map - min(disparity_map(:))+1;
[m,n]=size(D);
x=1:m;
y=1:n;
u=repmat(y,m,1);
v=repmat(x',1,n);
cam0=[3661.6 0 1500; 0 3663.9 987.51; 0 0 1];
cam1=[3661.6 0 1500; 0 3663.9 987.51; 0 0 1];
Tx=100;
doffs = 1500;
fx=cam0(1,1)*0.3;
fy=cam0(2,2)*0.3;
Ox=cam0(1,3)*0.3;
Oy=cam0(2,3)*0.3;
Z=fx*Tx./(D+doffs);

X=(u-Ox).*Z/fx;
Y=(v-Oy).*Z/fy;


I1_resize = imresize(I1,size(Z),'bilinear');

if do_plot
    x = double(I1_resize(:,:,1))/255;
    y = double(I1_resize(:,:,2))/255;
    z = double(I1_resize(:,:,3))/255;
    C(:,1) = x(:);
    C(:,2) = y(:);
    C(:,3) = z(:);
    figure;
    scatter3(X(:), Y(:), Z(:),1.5, C);
    title('3D-model');xlabel('X');ylabel('Y');zlabel('Z');
end
end
