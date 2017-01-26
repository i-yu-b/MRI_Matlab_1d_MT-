function cost = Mzfree_cost(omega1,offsets,data,T1,T2)
    R1 = 1/T1;
    fit = R1*(1+(2*pi*offsets*T2).^2)./ (R1*(1+(2*pi*offsets*T2).^2)+...
                                                              omega1^2*T2);
    cost = data - fit;
    semilogx(offsets,data,'ok',offsets,fit,'-r');
    title ('Magnetization transfer data. Henkelman approach')
    xlabel('\Delta\omega, kHz')
    ylabel('M_{free}/M_{free_{0}}')
    drawnow;
end