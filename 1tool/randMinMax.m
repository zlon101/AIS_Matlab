function arr = randMinMax(row, col, vmin, vmax)
arr = rand(row,col);
arr = arr.*(vmax-vmin)+vmin;
arr = round(arr,3);
end