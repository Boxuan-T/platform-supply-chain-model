function out = run_case(task, outputDirectory, solver)
if nargin < 3
    solver = 'ode45';
end
solver = validatestring(solver, {'ode45', 'ode15s'});
psc.validate_parameters(task.p);
validateattributes(task.U, {'double'}, {'real', 'finite', '2d', 'nrows', 4, 'nonempty'});
assert(all(task.U(1:3, :) > 0 & task.U(1:3, :) < 1, 'all'));
assert(all(task.U(4, :) >= 0));
assert(~isempty(regexp(task.key, '^[a-zA-Z0-9_-]+$', 'once')));

root = fileparts(mfilename('fullpath'));
sources = {'model', 'simulate_batch', 'equilibria', 'cooperation_metrics', 'run_case'};
stamp = '';
for k = 1:numel(sources)
    stamp = [stamp, fileread(fullfile(root, [sources{k}, '.m']))];
end
checkpointDirectory = fullfile(outputDirectory, 'checkpoints');
if ~isfolder(checkpointDirectory)
    mkdir(checkpointDirectory);
end
file = fullfile(checkpointDirectory, [task.key, '_', solver, '.mat']);
if isfile(file)
    old = load(file, 'task', 'stamp', 'out');
    if isequaln(old.task, task) && strcmp(old.stamp, stamp)
        out = old.out;
        return;
    end
end

n = size(task.U, 2);
R = zeros(n, 17);
trace = struct('t', [], 'u', []);
for first = 1:125:n
    ids = first:min(n, first + 124);
    if strcmp(solver, 'ode15s')
        [rows, path] = psc.simulate_batch(task.U(:, ids), task.p, solver, 1e-10, 1e-12);
    else
        [rows, path] = psc.simulate_batch(task.U(:, ids), task.p);
    end
    R(ids, :) = rows;
    if n == 1
        trace = path;
    end
end
out = struct('R', R, 'trace', trace, 'cooperation', psc.cooperation_metrics(task.p));
temporary = [tempname(checkpointDirectory), '.mat'];
cleanup = onCleanup(@() remove_temporary(temporary));
save(temporary, 'task', 'stamp', 'out');
movefile(temporary, file, 'f');
end

function remove_temporary(file)
if isfile(file)
    delete(file);
end
end
