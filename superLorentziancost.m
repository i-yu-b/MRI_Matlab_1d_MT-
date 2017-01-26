function cost = superLorentziancost(parameters, omega1, delta, data, varargin)
% cost function for absorption lineshape assuming superlorentzian
% inputs: omega1 - B1 amplitudes in Hz
%         delta - frequency offsets in Hz
%         data - experimental resuls, Mz/Mz0
%         parameters - [exchange rate, R1restricted, T2restricted, 
%                       R*Mor/R1f, 1/(R1f*T2f)]
R=parameters(1); %exchange rate
R1r=parameters(2); % R1 rate of restricted pool
% Rb=1;
T2r=parameters(3); %T2r relaxation time of restricted pool
coeff1=parameters(4); % R*Mor/R1f, where Mor - Mo value of restricted pool, 
                  %                  R1f - R1 rate of free pool
coeff2=parameters(5); % 1/(R1f*T2f)

% W - saturation rate of the restricted pool
W=zeros([length(delta), length(omega1)]);
for i = 1:length(delta)
    W(i,:) = pi*superlorentzian(delta(i),T2r) * omega1.^2;
end

numerator = R1r*coeff1 + W + R1r + R; 
denominator = coeff1*(R1r+W) + ... 
             (1 + (1./(2*pi*delta)).^2*omega1.^2*coeff2).*(R1r + W + R);
Mzf = numerator./denominator;

semilogx(delta/1000, data,'o',delta/1000, Mzf,'-k');
title ('Magnetization transfer data. Henkelman simulations')
xlabel('\Delta\omega, kHz')
ylabel('M_{free}/M_{free_{0}}')
lgd = legend(cellstr(num2str(round(omega1/(2*pi*1000),2)', '%0.2f')));
lgd.Location = 'northwest';
lgd.Box = 'off';
lgd.Title.String = 'RF amp attentuation, kHz';
drawnow

cost = Mzf-data;
