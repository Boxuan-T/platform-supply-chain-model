function tests = test_model
tests = functiontests(localfunctions);
end

function testJacobianAndAccounting(testCase)
rng(2606219);
p = psc.reference_parameters();
worst = 0;
for mode = {'full', 'fixed', 'instant'}
    for revenue = {'baseline', 'exp', 'hill'}
        p.mode = mode{1};
        p.revenue = revenue{1};
        for k = 1:4
            u = [0.1+0.8*rand(3, 1); 10*rand];
            [~, ~, payoffs, jacobian] = psc.model(u, p);
            n = size(jacobian, 1);
            numerical = zeros(n);
            for j = 1:n
                h = 1e-5*(1+abs(u(j)));
                plus = u;
                minus = u;
                plus(j) = plus(j)+h;
                minus(j) = minus(j)-h;
                a = psc.model(plus, p);
                b = psc.model(minus, p);
                numerical(:, j) = (a(1:n)-b(1:n))/(2*h);
            end
            worst = max(worst, max(abs(jacobian-numerical), [], 'all'));
            q = p;
            q.sP = 0.2;
            q.sM = 0.5;
            q.sR = 0.3;
            [~, ~, redistributed] = psc.model(u, q);
            verifyEqual(testCase, sum(payoffs), sum(redistributed), 'AbsTol', 1e-12);
        end
    end
end
verifyLessThan(testCase, worst, 1e-7);
end

function testEquilibriumThresholds(testCase)
p = psc.reference_parameters();
[states, alpha] = psc.equilibria(p);
verifyLessThan(testCase, max(abs(psc.model(states, p)), [], 'all'), 1e-12);
verifyTrue(testCase, alpha(1) < 0 && alpha(5) > 0 && alpha(8) < 0);
metrics = psc.cooperation_metrics(p);
verifyEqual(testCase, metrics.D, 9, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.ell, [13/54, 13/162, 13/162], 'AbsTol', 1e-12);
end

function testSolverAndBatchAgreement(testCase)
p = psc.reference_parameters();
states = psc.strategy_grid([0.2, 0.5, 0.8], 4.5);
a = psc.simulate_batch(states, p);
b = psc.simulate_batch(states, p, 'ode15s', 1e-10, 1e-12);
verifyEqual(testCase, a(:, 1), b(:, 1));
verifyLessThan(testCase, max(abs(a(:, 8)-b(:, 8))), 1e-5);
for j = [1, 14, 27]
    one = psc.simulate_batch(states(:, j), p);
    verifyEqual(testCase, one(1), a(j, 1));
    verifyEqual(testCase, one(8), a(j, 8), 'AbsTol', 1e-5);
end
end

function testMatchedDataAndAsymmetry(testCase)
tasks = psc.experiment_tasks('matched_data');
rows = cell(1, 2);
for k = 1:2
    rows{k} = psc.simulate_batch(tasks(k).U, tasks(k).p);
end
verifyEqual(testCase, rows{1}(:, 1), rows{2}(:, 1));
verifyEqual(testCase, sum(rows{1}(:, 1) == 8), 74);
verifyEqual(testCase, mean(rows{1}(:, 8)), 4.336620619735057, 'AbsTol', 1e-6);
verifyEqual(testCase, mean(rows{2}(:, 8)), 4.439635349604046, 'AbsTol', 1e-6);
tasks = psc.experiment_tasks('asymmetry');
verifyEqual(testCase, numel(tasks), 1395);
verifyEqual(testCase, [tasks(1).p.cM, tasks(466).p.cM, tasks(931).p.qR], [1, 1.5, 1.5]);
end

function testAllocationTies(testCase)
directory = tempname;
mkdir(directory);
cleanup = onCleanup(@() rmdir(directory, 's'));
tasks = psc.experiment_tasks('allocation');
parameters = [tasks.p];
shares = [[parameters.sP]', [parameters.sM]', [parameters.sR]'];
ids = abs(shares(:, 1)-1/3) < 1e-12 & shares(:, 2) >= 0.3-1e-12 & shares(:, 2) <= 11/30+1e-12;
tasks = tasks(ids);
outputs = cell(1, numel(tasks));
for k = 1:numel(tasks)
    outputs{k} = psc.run_case(tasks(k), directory);
end
result = psc.export_results(tasks, outputs, directory, 'allocation');
verifyEqual(testCase, numel(result.maxima.basinKeys), 3);
verifyEqual(testCase, result.maxima.maxBasin, 0.896, 'AbsTol', 1e-12);
verifyEqual(testCase, result.maxima.returnShares, [1/3, 1/3, 1/3], 'AbsTol', 1e-12);
verifyEqual(testCase, result.maxima.maxReturn, 6.7179935167, 'AbsTol', 1e-6);
end

function testExperimentParameters(testCase)
names = {'reference', 'baseline', 'initial_data', 'basin_sections', ...
    'accumulation', 'allocation', 'mechanisms', 'sensitivity', ...
    'adaptation', 'matched_data', 'asymmetry'};
for k = 1:numel(names)
    tasks = psc.experiment_tasks(names{k});
    verifyEqual(testCase, numel(unique({tasks.key})), numel(tasks));
    for j = 1:numel(tasks)
        psc.validate_parameters(tasks(j).p);
    end
end
end

function testExportAndCheckpoint(testCase)
directory = tempname;
mkdir(directory);
cleanup = onCleanup(@() rmdir(directory, 's'));
tasks = psc.experiment_tasks('baseline');
out = psc.run_case(tasks(1), directory);
cached = psc.run_case(tasks(1), directory);
verifyEqual(testCase, out, cached);
result = psc.export_results(tasks(1), {out}, directory, 'baseline');
verifyEqual(testCase, result.statistics.basin, 1);
verifyEqual(testCase, result.statistics.meanW, 8.1364, 'AbsTol', 5e-5);
table = readtable(fullfile(directory, 'baseline_trajectories.csv'));
verifyEqual(testCase, height(table), 1);
tasks(1).p.delta = 1.5;
changed = psc.run_case(tasks(1), directory);
verifyEqual(testCase, changed.R(1), 1);
end
