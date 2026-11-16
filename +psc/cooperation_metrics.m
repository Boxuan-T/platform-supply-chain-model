function c = cooperation_metrics(p)
[E, alpha] = psc.equilibria(p);
D = E(4, 8);
switch p.revenue
    case 'baseline'
        V = p.v*D/(p.K+D);
    case 'exp'
        V = p.v*(-expm1(-D/p.K));
    case 'hill'
        V = p.v*D^2/(p.K^2+D^2);
end
ell = [(p.CH-p.CL)/(p.kH-p.kL), p.cM/p.kH, p.cR/p.kH]/V;
c = struct('D', D, 'V', V, 'ell', ell, 'stable', alpha(8)<-1e-8, ...
'feasible', sum(ell)<1-1e-12, 'area', max(0, 1-sum(ell))^2, ...
'margin', -alpha(8), 'W', p.kH*V-p.CH-p.cM-p.cR, ...
'viable', alpha(8)<-1e-8 && p.kH*V*p.sP-p.CH>=-1e-12);
end
