function results = run_tests()
root = fileparts(mfilename('fullpath'));
results = runtests(fullfile(root, 'tests'));
assert(all([results.Passed]));
end
