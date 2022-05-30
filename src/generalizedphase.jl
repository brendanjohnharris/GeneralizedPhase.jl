rewrap(xp) = @. xp - 2*π*floor((xp-π)/(2*π)) - 2*π

interp1(x, y; k=3, bc="extrapolate", kw...) = Spline1D(x, y; k, bc, kw...) # 3rd order spline instead of hermite
interp1(y; kw...) = interp1(axes(x)..., y; kw...)

naninterp(y; kw...) = (idxs = .!isnan.(y); (interp1(findall(idxs), y[idxs]; kw...), idxs))
naninterp!(y; kw...) = ((f, idxs) = naninterp(y; kw...); y[.!idxs] .= f.(findall(.!idxs)); ())


function generalized_phase(x::AbstractVector, fs, lp)
    nwin = 3;
    @assert ndims(x) == 3 "Three-dimensional input required"
    [rows,cols,npts] = size(x);
    dt = 1/fs;

    # Analytic signal representation (single-sided Fourier approach, cf. Marple 1999). This is done for each channel
    x = reshape( x, [rows*cols,npts] )'; # Flatten spatially
    xo = fft( x, npts, 1 )
    h = zeros( npts, ~isempty(x) ); # Fourier transform each channel

    if ( npts > 0 ) && ( mod(npts,2) == 0 )
        h( [1 npts/2+1] ) = 1;
        h(2:npts/2) = 2;
    else
        h(1) = 1;
        h(2:(npts+1)/2) = 2;
    end
    xo = ifft( xo.*h(:,ones(1,size(x,2))), [], 1 ); xo = xo';
    xo = reshape( xo, [rows,cols,npts] );
    ph = angle( xo ); md = abs( xo );

    % calculate IF
    wt = zeros(size(xo));
    wt(:,:,1:end-1) = angle( xo(:,:,2:end) .* conj( xo(:,:,1:end-1) ) ) ./ (2*pi*dt);

    % rectify rotation
    sign_if = sign( nanmean(wt(:)) );
    if ( sign_if == -1 )
        modulus = abs(xo); ang = sign_if .* angle(xo); % rectify rotation
        xo = modulus .* exp( 1i .* ang ); ph = angle( xo ); md = abs( xo );
        wt(:,:,1:end-1) = ... % re-calculate IF
            angle( xo(:,:,2:end) .* conj( xo(:,:,1:end-1) ) ) ./ (2*pi*dt);
    end

    for ii = 1:rows
        for jj = 1:cols

            % check if nan channel
            if all( isnan(ph(ii,jj,:)) ), continue; end

            % find negative frequency epochs (i.e. less than LP cutoff)
            idx = ( squeeze(wt(ii,jj,:)) < lp ); idx(1) = 0; [L,G] = bwlabel( idx );
            for kk = 1:G
                idxs = find( L == kk );
                idx( idxs(1):( idxs(1) + ((idxs(end)-idxs(1))*nwin) ) ) = true;
            end

            % "stitch over" negative frequency epochs
            p = squeeze( ph(ii,jj,:) ); p(idx) = NaN;
            if all( isnan(p) ), continue; end % check if all NaNs
            p = unwrap(p)
            naninterp!(p)
            p = rewrap(p)
            ph(ii,jj,:) = p(1:size(ph,3));

        end
    end

    % save output datacube
    xgp = md .* exp( 1i .* ph );
