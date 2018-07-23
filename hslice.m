function out = hslice(ndarray, d, i)
    sz = size(ndarray);
    pfx = sz(1:d-1);    % dimensions before d
    sfx = sz(d+1:end);  % dimensions after d

    tmp = reshape(ndarray, prod(pfx), sz(d), prod(sfx));
    out = reshape(tmp(:, i, :), [pfx length(i) sfx]);
end