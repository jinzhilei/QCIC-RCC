function main_basic()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% QCIC model main program (basic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set: time step
t_start = 0;
t_end = 365;
t_vec = t_start:1:t_end;

%% Set: system initialization values
x0=initialzation_basic();   

%% Calculation: QCIC model by funct ion ode45
[t, x] = ode45(@(t,x) QCIC_basic(t,x), t_vec, x0);

end