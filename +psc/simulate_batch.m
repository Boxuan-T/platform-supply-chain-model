function [R, trace] = simulate_batch(U, p, solver, rtol, atol)
if nargin<3
    solver = 'ode45';
end
if nargin<4
    rtol = 1e-8;
end
if nargin<5
    atol = 1e-10;
end
N = size(U, 2);
assert(all(U(1:3, :)>0 & U(1:3, :)<1, 'all'));
[E, alpha] = psc.equilibria(p);
scale = [1;1;1;1+p.Dref];
L = log(U(1:3, :)./(1-U(1:3, :)));
Y = [L;U(4, :);zeros(3, N)];
R = nan(N, 17);
R(:, 15:17) = 0;
active = 1:N;
prevU = [];
prevT = [];
trace = struct('t', [], 'u', []);
opts = odeset('RelTol', rtol, 'AbsTol', atol);
for block = 1:10
    t0 = (block-1)*500;
    t1 = block*500;
    n = numel(active);
    if strcmp(solver, 'ode15s')
        pattern = kron(speye(n), sparse(ones(7)));
        opts = odeset(opts, 'JPattern', pattern);
    end
    fn = @(t, v) logit_rhs(t, v, p, n);
    ts = (t0:.5:t1)';
    [ts, ys] = feval(solver, fn, ts, Y(:), opts);
    if numel(ts)~=1001 || ts(end)~=t1 || any(~isfinite(ys), 'all')
        error('psc:IntegrationHorizon');
    end
    Z = reshape(ys', 7, n, []);
    P = 1./(1+exp(-Z(1:3, :, :)));
    D = Z(4, :, :);
    if strcmp(p.mode, 'fixed')
        D(:) = p.Dref;
    end
    if strcmp(p.mode, 'instant')
        D = (p.aL+(p.aH-p.aL)*P(1, :, :)).*(p.qM*P(2, :, :)+p.qR*P(3, :, :)+p.qC*P(2, :, :).*P(3, :, :))/p.delta;
    end
    A = cat(1, P, D);
    if any(D < -1e-8, 'all')
        error('psc:NegativeData');
    end
    if block==1
        av = squeeze(Z(5:7, :, 201))/100;
        if n==1
            av = reshape(av, 3, 1);
        end
        R(active, 9:11) = av';
        R(active, 8) = sum(av, 1)';
    end
    if N==1
        trace.t = [trace.t;ts(1+(block>1):end)];
        trace.u = cat(3, trace.u, A(:, :, 1+(block>1):end));
    end
    if ~isempty(prevU)
        Atest = cat(3, prevU, A(:, :, 2:end));
        ttest = [prevT;ts(2:end)];
    else
        Atest = A;
        ttest = ts;
    end
    flat = reshape(Atest, 4, []);
    du = psc.model(flat, p);
    residual = reshape(max(abs(du)./scale, [], 1), n, []);
    done = false(1, n);
    for j = 1:n
        u = A(:, j, end);
        [~, idx] = min(max(abs(E-u)./scale, [], 1));
        eq = E(:, idx);
        al = alpha(idx);
        distEnd = max(abs(eq-u)./scale);
        if distEnd>1e-6 && residual(j, end)<1e-5
            [eq, al, ok] = mixed_root(u, p);
            if ok
                idx = 9;
            else
                eq = E(:, idx);
                al = alpha(idx);
            end
        end
        dist = squeeze(max(abs(Atest(:, j, :)-eq)./scale, [], 1))';
        good = (dist<=1e-6)&(residual(j, :)<=1e-7);
        runs = conv(double(good), ones(1, 21), 'valid');
        first = find(runs>=21, 1);
        stable = al < -1e-8;
        if stable && ~isempty(first) && all(good(end-20:end))
            R(active(j), 1) = idx;
            R(active(j), 2) = ttest(first);
            done(j) = true;
        elseif block==10
            R(active(j), 1) = 0;
            if all(good(end-20:end))
                if al>1e-8
                    R(active(j), 1) = -2;
                else
                    R(active(j), 1) = -1;
                end
            end
            done(j) = true;
        end
        if done(j)
            R(active(j), 3:6) = u';
            R(active(j), 7) = al;
            R(active(j), 12) = residual(j, end);
            R(active(j), 13) = t1;
        end
    end
    if all(done)
        break;
    end
    keep = ~done;
    prevU = A(:, keep, end-19:end);
    prevT = ts(end-19:end);
    Y = Z(:, keep, end);
    active = active(keep);
end
R(:, 14) = abs(R(:, 8)-sum(R(:, 9:11), 2));
end

function dv = logit_rhs(~, v, p, n)
Y = reshape(v, 7, n);
s = 1./(1+exp(-Y(1:3, :)));
U = [s;Y(4, :)];
[du, F, pay] = psc.model(U, p);
dL = [p.nuP;p.nuM;p.nuR].*F;
if ~strcmp(p.mode, 'full')
    du(4, :) = 0;
end
dv = reshape([dL;du(4, :);pay], [], 1);
end

function [eq, alpha, ok] = mixed_root(u, p)
eq = u;
free = find(u(1:3)>1e-5 & u(1:3)<1-1e-5);
fixed = setdiff(1:3, free);
eq(fixed) = round(eq(fixed));
ok = false;
alpha = nan;
if isempty(free)
    return;
end
ids = [free(:);4];
if ~strcmp(p.mode, 'full')
    ids = free(:);
end
for it = 1:40
    f = rootfun(eq, p, free);
    if norm(f, inf)<1e-10
        ok = true;
        break;
    end
    H = zeros(numel(f), numel(ids));
    for k = 1:numel(ids)
        h = 1e-5*(1+abs(eq(ids(k))));
        ep = eq;
        em = eq;
        ep(ids(k)) = ep(ids(k))+h;
        em(ids(k)) = em(ids(k))-h;
        H(:, k) = (rootfun(ep, p, free)-rootfun(em, p, free))/(2*h);
    end
    if rcond(H)<1e-13
        break;
    end
    step = H\f;
    accepted = false;
    for damp = 2.^-(0:15)
        en = eq;
        en(ids) = eq(ids)-damp*step;
        if all(en(1:3)>=0 & en(1:3)<=1) && en(4)>=0 && norm(rootfun(en, p, free), inf)<norm(f, inf)
            eq = en;
            accepted = true;
            break;
        end
    end
    if ~accepted
        break;
    end
end
if strcmp(p.mode, 'fixed')
    eq(4) = p.Dref;
end
if strcmp(p.mode, 'instant')
    eq(4) = (p.aL+(p.aH-p.aL)*eq(1))*(p.qM*eq(2)+p.qR*eq(3)+p.qC*eq(2)*eq(3))/p.delta;
end
if ok
    [f, ~, ~, J] = psc.model(eq, p);
    ok = max(abs(f))<1e-8;
    alpha = max(real(eig(J)));
end
end
function f = rootfun(eq, p, free)
[du, F] = psc.model(eq, p);
f = F(free);
if strcmp(p.mode, 'full')
    f = [f;du(4)];
end
end
