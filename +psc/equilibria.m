function [E, alpha] = equilibria(p)
E = zeros(4, 8);
alpha = zeros(1, 8);
for i = 1:8
    s = double(bitget(i-1, 3:-1:1))';
    D = (p.aL+(p.aH-p.aL)*s(1))*(p.qM*s(2)+p.qR*s(3)+p.qC*s(2)*s(3))/p.delta;
    if strcmp(p.mode, 'fixed')
        D = p.Dref;
    end
    E(:, i) = [s;D];
    [~, ~, ~, J] = psc.model(E(:, i), p);
    alpha(i) = max(real(eig(J)));
end
end
