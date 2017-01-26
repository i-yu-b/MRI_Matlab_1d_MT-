function cost = T1_IR_cost(TI, data, parameters,flag)
% cost function for absorption lineshape assuming superlorentzian
% inputs: TI - delay after Inversion pulse in s
%         data - experimental resuls, don't have to be normalized
%         parameters - fit parameters
%         if flag == 1 - monoexp fit 
%         if flag == 2 - biexp fit
    if flag == 1 
        a = parameters(1);
        b = parameters(2);
        T1f = parameters(3); %T1 observed of free pool
        fit = a-b*exp(-TI/T1f);
        semilogy(TI,abs(a-data),'ok',TI,abs(a-fit),'-r'); 
        title('Inversion-Recovery data monoexp fit')
        xlabel('TI, s')
        ylabel('log(M_{free})')
        drawnow;
    
    end
    
    if flag == 2
        a = parameters(1);
        b = parameters(2);
        T1f = parameters(3); %T1 observed of free pool
        c = parameters(4);
        T1r = parameters(5); %T1 observed of restricted pool
        fit = a-b*exp(-TI/T1f)-c*exp(-TI/T1r);
        semilogy(TI,abs(a-data),'ok',TI,abs(a-fit),'-r'); 
        title('Inversion-Recovery data biexp fit')
        xlabel('TI, s')
        ylabel('log(M_{free})')
        drawnow;
    end
    
    
    cost = fit - data;
end