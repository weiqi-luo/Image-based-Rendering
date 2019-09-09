function [t1,t2] = epipolar_rectification( f, pts1, pts2, imageSize )

    %%  Compute the projective transformation for the second camera.

    % Find the epipole
    [u, d, v] = svd(f);
    epipole = u(:, 3);

    % Translate the epipole to put the origin at the center of the image
    t = [1, 0, -imageSize(2)/2; ...
         0, 1, -imageSize(1)/2; ...
         0, 0, 1];
    epipole = t * epipole;

    % Move the epipole to the line at y=0 by rotating the image with the minimum angle.
    % Ensure the homogeneous coordinates have positive sign.
    if epipole(3) < 0
      epipole = -epipole;
    end

    % If the x coordinate is positive, rotate clockwise; otherwise, rotate counter-clockwise.
    if epipole(1) >= 0
      keepOrientation = -1;
    else
      keepOrientation = 1;
    end
    r = [-epipole(1),  -epipole(2),                                     0;
          epipole(2),  -epipole(1),                                     0;
           0,      0,   sqrt(epipole(2)^2+epipole(1)^2) * keepOrientation];
    epipole = r * epipole;

    % Move the epipole to the position at infinit in the x-axis direction
    if epipole(1) ~= 0
      g = [1,            0, 0; 
           0,            1, 0; 
       -epipole(3)/epipole(1),  0, 1];
    else
      g = eye(3);
    end

    % Compute the overall transformation for the second image.
    t2 = g * r * t;
    t([1,2], 3) = -t([1,2], 3);
    t2 = t * t2;
    if(t2(end)) ~= 0
      t2 = t2 / t2(end);
    end

    %% Compute the projective transformation for the first camera.

    % Compute a (partial) synthetic camera matrix for the second camera.
    z = [ 0, -1, 0;
          1,  0, 0;
          0,  0, 1 ];
    d(3,3) = (d(1,1) + d(2,2)) / 2;
    m = u * z * d * v';

    % Compute the transformation for the first camera so the rows in the first
    % camera correspond the rows at the same location in the second camera.
    h1r = t2 * m;

    % Compute the transformation for the first camera so the points have
    % approximately the same column locations in the first and second images.
    numPoints = size(pts1, 1);
    p1 = h1r * [pts1'; ones(1, numPoints)];
    p2 = t2  * [pts2'; ones(1, numPoints)];

    b = p2(1,:) ./ p2(3,:) .* p1(3,:);
    y = p1' \ b';
    h1c = [y(1), y(2), y(3);          
           0,    1,    0;
           0,    0,    1];

    % Compute the overall transformation for the first camera.
    t1 = h1c * h1r;
    if(t1(end)) ~= 0
      t1 = t1 / t1(end);
    end

    t1 = t1';
    t2 = t2';
end