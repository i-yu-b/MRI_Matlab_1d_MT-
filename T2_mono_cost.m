function cost = T2_mono_cost(TE, data, parameters)
% cost function for T2 relaxation data
% inputs: TE - time delay, s
%         data - experimental resuls, don't have to be normalized
%         parameters - fit parameters
    a = parameters(1);
    T2f = parameters(2); %T2 observed of free pool
    fit = a*exp(-TE/T2f);
    semilogy(TE,data,'ok',TE,fit,'-r'); 
    title('T2 relaxation data, monoexp fit')
    xlabel('TE, s')
    ylabel('log(M_{free})')
    drawnow;
    cost = fit - data;
end