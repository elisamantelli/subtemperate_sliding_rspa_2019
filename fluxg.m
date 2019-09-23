function [fout,Dfout] = fluxg(H_g,parameters)
%flux condition at the grounding line and its derivative

fout.flux = parameters.gr.k*parameters.gamma^(-1/3)*H_g^parameters.gr.m;
if nargout == 2
Dfout.dflux_dHg = parameters.gr.k*parameters.gamma^(-1/3)*parameters.gr.m*H_g^(parameters.gr.m-1);   
end
end