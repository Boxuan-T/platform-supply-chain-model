function [du, F, pay, J] = model(u, p)
x = u(1, :);
y = u(2, :);
z = u(3, :);
D = u(4, :);
a = p.aL+(p.aH-p.aL)*x;
Q = p.qM*y+p.qR*z+p.qC*y.*z;
if strcmp(p.mode, 'fixed')
    D = 0*D+p.Dref;
end
if strcmp(p.mode, 'instant')
    D = a.*Q/p.delta;
end
switch p.revenue
    case 'baseline'
        V = p.v*D./(p.K+D);
        Vp = p.v*p.K./(p.K+D).^2;
    case 'exp'
        V = p.v*(-expm1(-D/p.K));
        Vp = p.v/p.K*exp(-D/p.K);
    case 'hill'
        V = p.v*D.^2./(p.K^2+D.^2);
        Vp = 2*p.v*p.K^2*D./(p.K^2+D.^2).^2;
end
WM = p.sP+p.sM;
WR = p.sP+p.sR;
PM = p.sP*p.wM/WM;
PR = p.sP*p.wR/WR;
MM = p.sM*p.wM/WM;
RR = p.sR*p.wR/WR;
AP = p.sP*y.*z+PM*y.*(1-z)+PR*(1-y).*z;
AM = MM*(1-z)+p.sM*z;
AR = RR*(1-y)+p.sR*y;
dk = p.kH-p.kL;
dc = p.CH-p.CL;
k = p.kL+dk*x;
F = [dk*V.*AP-dc; k.*V.*AM-p.cM; k.*V.*AR-p.cR];
r = [p.nuP;p.nuM;p.nuR].*u(1:3, :).*(1-u(1:3, :));
du = [r.*F; a.*Q-p.delta*D];
if strcmp(p.mode, 'fixed')
    du(4, :) = 0;
end
if strcmp(p.mode, 'instant')
    grad = [(p.aH-p.aL)*Q; a.*(p.qM+p.qC*z); a.*(p.qR+p.qC*y)]/p.delta;
    du(4, :) = sum(grad.*du(1:3, :), 1);
end
if nargout>=3
    pay = [k.*V.*AP-(p.CL+dc*x); y.*F(2, :); z.*F(3, :)];
end
if nargout>=4
    assert(size(u, 2)==1, 'psc:JacobianState');
    APy = PM*(1-z)+(p.sP-PR)*z;
    APz = PR*(1-y)+(p.sP-PM)*y;
    J = [p.nuP*(1-2*x)*F(1), r(1)*dk*V*APy, r(1)*dk*V*APz, r(1)*dk*Vp*AP; ...
    r(2)*dk*V*AM, p.nuM*(1-2*y)*F(2), r(2)*k*V*(p.sM-MM), r(2)*k*Vp*AM; ...
    r(3)*dk*V*AR, r(3)*k*V*(p.sR-RR), p.nuR*(1-2*z)*F(3), r(3)*k*Vp*AR; ...
    (p.aH-p.aL)*Q, a*(p.qM+p.qC*z), a*(p.qR+p.qC*y), -p.delta];
    if strcmp(p.mode, 'fixed')
        J = J(1:3, 1:3);
    end
    if strcmp(p.mode, 'instant')
        J = J(1:3, 1:3)+J(1:3, 4)*grad';
    end
end
end
