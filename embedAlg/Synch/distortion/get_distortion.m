function y = get_distortion( cover, stego, rho )

s = double(stego) - double(cover);
y = 0;

for i = 1:size(cover, 1)
    for j = 1:size(cover, 2)
        s_patch   = padded_patch( s, i, j );
        if sum( s_patch(:) )           
            rho_patch = padded_patch( rho, i, j );           
            y = y + D( s_patch, rho_patch );
        end
    end
end

end