function fout = disp_analytic_v2(lin, parameters)

%physical parameters
alpha = parameters.alpha;
gamma = parameters.gamma;
DD = 1; %already accounted for by DF in the numerical regularization

%sliding law
[f_slide_Tbed, dfslide_Tbed] = regularization_old(parameters.T_bed(1), parameters);
parameters.T_bed(1)
F0 = f_slide_Tbed; 
DF = dfslide_Tbed; 

%steady state 
u0 = 3*F0/(gamma + 3*F0);
tau0 = ( F0 *gamma^-1 )^-1 *u0;
Q0 = parameters.nu + alpha*tau0*u0;

%initialization
klist = lin.klist;
disp = zeros(5, length(klist));
res_plot = zeros(5, length(klist));
disp_sqrt = zeros(5, length(klist));
sigma = NaN*zeros(2, length(klist));
sqrt_ice = NaN*zeros(2, length(klist));
sqrt_bed = NaN*zeros(2, length(klist));
for jj = [1:length(klist)/2 length(klist)/2+2:length(klist)]
    k = klist(jj);
    
    %coefficients of the dispersion relationship for ~k mode and T-dep
    %sliding with  the ice thickness scale velocity field (verified against
    %mathematica)
    % disp_or = a0 + a12*( 1i* k*u0 + sigma)*(-sigma^{1/2} - ( 1i* k*u0 + sigma)^(1/2)) +a1*( 1i* k*u0 + sigma);
    if k ==0
        ch_sh_m1 =cosh(k)*sinh(k);
    else
    ch_sh_m1 = (1-tanh(k).^2)./tanh(k); %alterante form for cosh(k)*sinh(k) well behaved in the large k limit;
    end
    a0 = 2*1i*k*Q0*(k.*ch_sh_m1-1);
    a1 = (-2*k.*alpha .*(tau0.*ch_sh_m1  +2*u0*tanh(k)) +alpha*tau0*2);
    a12 = DD*(-2*k*(F0+gamma).*ch_sh_m1 +2*F0*k*(coth(k) +tanh(k)) +2*gamma)*(DF*tau0)^(-1);
    p = 1i*k*u0;
    
%      a0 = 2*1i*k*Q0*(k- cosh(k).*sinh(k));
%      a1 = (-2*k.*alpha .*(tau0  +2*u0*k*sinh(k)^2) +alpha*tau0*sinh(2*k));
%      a12 = DD*(-2*k*(F0+gamma) +2*F0*k*cosh(2*k) +sinh(2*k)*gamma)*(DF*tau0)^(-1);

    
    c0 = a0^2 + 2*a0*a1*p + p^2*(a1^2 -a12^2*p);
    c1 = -2*a12*p*(a0 +a1*p);
    c2 = 2*(a0*a1 + p*(a1^2 - a12^2*p));
    c3 = -2*a12*(a0 + 2*a1*p);
    c4 = a1^2 - a12^2*p;
    c5 = -2*a1*a12;
    
    %compute y = (sigma)^{1/2}
    y = roots([c5;c4;c3;c2;c1;c0]);
    [~,I]= sort(real(y),'descend');
    y = y(I);
    
    disp_sqrt(1:length(y),jj) = y.^2;
    
    sigma2_12_p = (y.^2+p).^(1/2);
    sigma2_12_m = -(y.^2+p).^(1/2);
    
    disp_or = @(y, sigma2_12) a0 + a12*sigma2_12.^2.*(-y - sigma2_12) +a1*sigma2_12.^2;
    
    res_p = disp_or(y, sigma2_12_p);
    res_m = disp_or(y, sigma2_12_m);
    
    index_p = find(abs(res_p)./(max(abs(res_p)))>1e-10);
    index_m = find(abs(res_m)./(max(abs(res_m)))>1e-10);
    
    index_p = sparse(index_p, ones(length(index_p),1),ones(length(index_p),1),5,1);
    sigma2_12_p(index_p == 1) = NaN + NaN*1i;
    index_m = sparse(index_m, ones(length(index_m),1),ones(length(index_m),1),5,1);
    sigma2_12_m(index_m==1) = NaN+ NaN*1i;
    disp(1:length(y),jj)= y(:).^2;
    
    index_bed = find(real(y)<0);
    disp(index_bed,jj) = NaN+1i*NaN;
    sigma2_12_p(index_bed) = NaN+1i*NaN;
    sigma2_12_m(index_bed) = NaN+1i*NaN;
    
    
    index_ice = zeros(5,1);
    for ss = 1:length(y)
        if isnan(sigma2_12_p(ss)) == 1 && real(sigma2_12_m(ss))<0
            index_ice(ss) = 1;
            res_plot(ss,jj) =  NaN+1i*NaN;
            sigma2_12_m(ss) = NaN+1i*NaN;
        elseif isnan(sigma2_12_m(ss)) == 1 && real(sigma2_12_p(ss))<0
            index_ice(ss)=1;
            res_plot(ss,jj) =  NaN+1i*NaN;
            sigma2_12_p(ss) = NaN+1i*NaN;
        elseif real(sigma2_12_p(ss))<0 && real(sigma2_12_m(ss))<0
            index_ice(ss)=1;
            res_plot(ss,jj) =  NaN+1i*NaN;
            sigma2_12_p(ss) = NaN+1i*NaN;
            sigma2_12_m(ss) = NaN+1i*NaN;
        elseif isnan(sigma2_12_p(ss)) == 1 && isnan(sigma2_12_m(ss)) == 1
            index_ice(ss)=1;
            res_plot(ss,jj) =  NaN+1i*NaN;
        else
            index_ice(ss)=0;
            res_plot(ss,jj) =  min(abs(res_p(ss)),abs(res_m(ss)));
            if abs(res_p(ss))<abs(res_m(ss))
                sigma2_12_m(ss) = NaN+1i*NaN;
            else
                sigma2_12_p(ss) = NaN+1i*NaN;
            end
        end
    end
    
    disp(index_ice==1,jj) = NaN+1i*NaN;
    
    %extract non NaN eigenvalues
    I = find(isnan(disp(:,jj))==0);
    if length(I) == 1,
        sigma(1,jj) = disp(I,jj);
        sqrt_ice(1,jj) = max(sigma2_12_m(I), sigma2_12_p(I));
        sqrt_bed(1,jj) = y(I);
    elseif length(I) == 2,
        sigma(:,jj) = disp(I,jj);
        sqrt_ice(1,jj) = max(sigma2_12_m(I(1)), sigma2_12_p(I(1)));
        sqrt_ice(2,jj) = max(sigma2_12_m(I(2)), sigma2_12_p(I(2)));
        sqrt_bed(:,jj) = y(I);
    end
    
end
sigma(:,length(klist)/2+1) = [0;NaN];
fout. sigma = sigma; 
fout. sqrt_ice = sqrt_ice;
fout. sqrt_bed = sqrt_bed;

%compute amplitudes
AQ = zeros(2,length(klist));
Atheta = zeros(2,length(klist));
for kk = 1:2
   AQ(kk,:) = (1i*Q0*klist.*(-2*klist + sinh(2*klist))./(1i*u0*klist + sigma(kk,:)))./sinh(2*klist); 
   Atheta(kk,:) = ((-2*klist.*(F0.*(1- cosh(2*klist)) + gamma) + gamma * sinh(2*klist))./(tau0*DF))./sinh(2*klist);
end
Atheta(isnan(sigma)==1) = NaN+ 1i*NaN;
AQ(isnan(sigma)==1) = NaN+ 1i*NaN;
fout.Atheta = Atheta;
fout.AQ = AQ; 



