function gray_image = rgb_to_gray(input_image)
    if numel(size(input_image)) > 2
        gray_image = double(input_image);
        gray_image = 0.299 * gray_image(:,:,1) + 0.587 * gray_image(:,:,2) +  0.114 * gray_image(:,:,3);
        gray_image = uint8(gray_image);
    else 
        gray_image = input_image;
    end
end