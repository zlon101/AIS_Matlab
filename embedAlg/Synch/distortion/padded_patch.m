function y = padded_patch( image, i, j )

y = zeros(3);

if i == 1
    if j == 1
        y(1:2,1:2) = image(i,j);
        y(1:2,3)   = image(i,j+1);
        y(3,1:2)   = image(i+1,j);
        y(3,3)     = image(i+1,j+1);
    elseif j == size(image,2)
        y(1:2,1)   = image(i,j-1);
        y(1:2,2:3) = image(i,j);
        y(3,1)     = image(i+1,j-1);
        y(3,2:3)   = image(i+1,j);
    else
        y(1:2,1)   = image(i,j-1);
        y(1:2,2)   = image(i,j);
        y(1:2,3)   = image(i,j+1);
        y(3,1)     = image(i+1,j-1);
        y(3,2)     = image(i+1,j);
        y(3,3)     = image(i+1,j+1);
    end
elseif i == size(image,1)
    if j == 1
        y(2:3,1:2) = image(i,j);
        y(1,1:2)   = image(i-1,j);
        y(1,3)     = image(i-1,j+1);
        y(2:3,3)   = image(i,j+1);
    elseif j == size(image,2)
        y(2:3,2:3) = image(i,j);
        y(1,2:3)   = image(i-1,j);
        y(1,1)     = image(i-1,j-1);
        y(2:3,1)   = image(i,j-1);
    else
        y(2:3,1)   = image(i,j-1);
        y(2:3,2)   = image(i,j);
        y(2:3,3)   = image(i,j+1);
        y(1,1)     = image(i-1,j-1);
        y(1,2)     = image(i-1,j);
        y(1,3)     = image(i-1,j+1);
    end
else
    if j == 1
        y(1,1:2)   = image(i-1,j);
        y(2,1:2)   = image(i,j);
        y(3,1:2)   = image(i+1,j);
        y(1,3)     = image(i-1,j+1);
        y(2,3)     = image(i,j+1);
        y(3,3)     = image(i+1,j+1);
    elseif j == size(image,2)
        y(1,2:3)   = image(i-1,j);
        y(2,2:3)   = image(i,j);
        y(3,2:3)   = image(i+1,j);
        y(1,1)     = image(i-1,j-1);
        y(2,1)     = image(i,j-1);
        y(3,1)     = image(i+1,j-1);
    else
        y          = image(i-1:i+1,j-1:j+1);
    end
end

end