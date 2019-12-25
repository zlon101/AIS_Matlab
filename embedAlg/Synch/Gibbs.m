function [stego_all, distortion_all, costs_all, lambda_all, grid] = Gibbs( cover, cover_rho, stego, sweeps, payload )

assert( isa( cover, 'uint8' ) );
assert( all(size(cover)>=[1 1]) );
assert( isa( stego, 'uint8' ) );
assert( all(size(stego)>=[1 1]) );
assert( isa( cover_rho, 'double' ) );
assert( all(size(cover_rho)>=[1 1]) );
assert( isa( sweeps, 'double' ) );
assert( all(size(sweeps)==[1 1]) );
assert( isa( payload, 'double' ) );
assert( all(size(payload)==[1 1]) );

cover = double( cover );
stego = double( stego );

wetCost = 10^10;

rhoM1_all = zeros( size(cover,1), size(cover,2), sweeps+1 );
rho0_all = rhoM1_all;
rhoP1_all = rhoM1_all;

% sweep 0
stego_all = zeros( size(cover,1), size(cover,2), sweeps+1, 'uint8' );
stego_all(:,:,1) = uint8( stego );
distortion_all = zeros( 1, sweeps+1 );
distortion_all(1) = get_distortion( cover, stego, cover_rho );

% grids
grid = ones( size( cover ) );
grid( 1:2:end, 1:2:end ) = 2;
grid( 2:2:end, 2:2:end ) = 2;
grids = max( grid(:) );

lambda_all = zeros( sweeps+1, grids );

% sweep loop
for s_i = 1:sweeps
    rM1 = zeros( size( cover ) );
    r0 = zeros( size( cover ) );
    rP1 = zeros( size( cover ) );
    % grids loop
    for g_i = 1:grids
        % locate pixels
        bool_indices = grid == g_i;
        [rows, cols] = find( bool_indices );
        rows = rows'; cols = cols';
        % remove changes
        stego( bool_indices ) = cover( bool_indices );
        
        % new costs loop
        z = zeros( size( cols ) );
        rhoM1 = z;
        rho0  = z;
        rhoP1 = z;
        for c_i = 1:size( cols, 2 )
            cover_patch = padded_patch( cover, rows(c_i), cols(c_i) );
            stego_patch = padded_patch( stego, rows(c_i), cols(c_i) );
            rho_patch   = padded_patch( cover_rho, rows(c_i), cols(c_i) );
            
            s = stego_patch - cover_patch;
            impulse = [ 0 0 0; 0 1 0; 0 0 0];
            
            if cover_patch(2,2) > 0
                rhoM1( c_i ) = D( s - impulse, rho_patch );
            else
                rhoM1( c_i ) = wetCost;
            end
            rho0( c_i )  = D( s, rho_patch );
            if cover_patch(2,2) < 255
                rhoP1( c_i ) = D( s + impulse, rho_patch );
            else
                rhoP1( c_i ) = wetCost;
            end
        end
        
        % new changes        
        [pattern, l] = EmbeddingSimulator( z, rhoM1, rho0, rhoP1, payload );
        lambda_all( s_i+1, g_i ) = l; 
        stego( bool_indices ) = cover( bool_indices ) + pattern';
        
        rM1( bool_indices ) = rhoM1;
        r0( bool_indices )  = rho0;
        rP1( bool_indices ) = rhoP1;
    end
    stego_all(:,:,s_i+1) = uint8( stego );
    distortion_all(s_i+1) = get_distortion( cover, stego, cover_rho );
    [rhoM1_all(:,:,s_i+1), rho0_all(:,:,s_i+1), rhoP1_all(:,:,s_i+1)] = deal( rM1, r0, rP1 );
end

costs_all.cover_rho = cover_rho;
costs_all.rhoM1 = rhoM1_all;
costs_all.rho0 = rho0_all;
costs_all.rhoP1 = rhoP1_all;
end
