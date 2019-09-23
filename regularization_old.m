function [fout, Dfout] = regularization_old(T_bed, parameters)
epsilon = parameters.reg.epsilon; 
epsilon_f = parameters.reg.epsilon_f;

% H = @(x) 1/2*(1 + x./sqrt(epsilon^2+x.^2));
% dH = @(x) epsilon^2./(2*(epsilon^2+x.^2).^(3/2));

H = @(x)  1/2*(1 + tanh(x/epsilon));
dH = @(x) (sech(x/epsilon)).^2./(2*epsilon);%(epsilon+epsilon*cosh(2*x/epsilon)).^(-1);

H_ev = H(-T_bed);
dH_ev = dH(-T_bed);
% chi = @(x) 1/2*(x+sqrt(epsilon^2+x.^2) );
% dchi = @(x) H(x);

chi = @(x) x/2 + epsilon/2*log(cosh(x/epsilon));
dchi = @(x) H(x);

chi_ev = chi(-T_bed);
dchi_ev = dchi(-T_bed);

[check_inf, ind_check_inf] = find(isinf(chi_ev(-T_bed>=0)) == 1);

if isempty(check_inf) == 0
    chi_ev(1:ind_check_inf(end)) = chi_ev(ind_check_inf(end)+1);    
end

f = @(x) exp(x/epsilon_f);
df = @(x) 1/epsilon_f*exp(x/epsilon_f);

fout =  f(0) - df(0).*chi_ev +(f(T_bed) -f(0) -df(0).*T_bed).*H_ev;
Dfout = ( df(0).*dchi_ev) + (df(T_bed) -df(0)).*H_ev - (f(T_bed) -f(0) -df(0).*T_bed).*dH_ev;



