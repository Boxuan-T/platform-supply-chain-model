function validate_parameters(p)
positive = {'aL', 'aH', 'delta', 'v', 'K', 'kL', 'kH', ...
    'CH', 'cM', 'cR', 'sP', 'nuP', 'nuM', 'nuR', 'Dref'};
nonnegative = {'qM', 'qR', 'qC', 'wM', 'wR', 'CL', 'sM', 'sR'};
for k = 1:numel(positive)
    validateattributes(p.(positive{k}), {'double'}, ...
        {'real', 'finite', 'scalar', 'positive'});
end
for k = 1:numel(nonnegative)
    validateattributes(p.(nonnegative{k}), {'double'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'});
end
assert(p.aH > p.aL && p.kH > p.kL && p.CH > p.CL);
assert(p.qM > 0 && p.qR > 0);
assert(p.wM + p.wR <= 1);
assert(abs(p.sP + p.sM + p.sR - 1) < 1e-12);
assert(any(strcmp(p.mode, {'full', 'fixed', 'instant'})));
assert(any(strcmp(p.revenue, {'baseline', 'exp', 'hill'})));
end
