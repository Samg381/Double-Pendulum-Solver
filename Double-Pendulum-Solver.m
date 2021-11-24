%{
 Double Pendulum Solver
 --------------------------------------------------------
 This script accepts four inputs: a length and mass for both pendulum arms. 
 The length inputs (L1, L2) are guesses; starting points for optomization.
 The script will attempt to match the specified masses with the smallest
 possible arm length to achieve natural frequencies of 2Hz and 5Hz.
 --------------------------------------------------------
%}

L1 = 0.15;  % Guess of Length 1 (Meters)
L2 = 0.1;  % Guess of Length 2 (Meters)
m1 = 0.19; % Mass of Arm 1 (Kg)
m2 = 0.12; % Mass of Arm 2 (Kg)

[wn] = Pendulum(m1,m2,L1,L2);
NF_Hz = wn/(2*pi);  % Convert natural frequency between Hz and Rad/s

length = fminsearch(@(e) Pendulum_Cost(e, m1,m2), [L1, L2]);  % Minimize L1 and L2

[wn,MS] = Pendulum(m1,m2,length(1),length(2));
NF_Hz = wn/(2*pi);  % Convert natural frequency between Hz and Rad/s
fprintf('Optimized Natural Frequency 1: %.3f Hz \n', NF_Hz(1))
fprintf('Optimized Natural Frequency 2: %.3f Hz \n', NF_Hz(2))
disp('Optimized Mode Shapes:')
disp(MS)
L1 = length(1); % Set L1 to the optomized value
L2 = length(2); % Set L1 to the optomized value
fprintf('Optimized Length 1: %.3f Meters \n', length(1))
fprintf('Optimized Length 2: %.3f Meters \n', length(2))

function [e] = Pendulum_Cost(x,M1,M2)
    L1 = x(1); % Unpack L1 from x(1)
    L2 = x(2); % Unpack L2 from x(2)
    [wn] = Pendulum(M1,M2,L1,L2); % Execute pendulum with optomized L1 and L2
    e = ((wn(1)-12.57)^2 + (wn(2)-31.42)^2); % Determine error between wn and the target
end

function [wn,MS] = Pendulum(m1,m2,L1,L2)
    g = 9.81; % Gravitational constant
    M = [m1 0;0 m2]; % Mass matrix
    K = [(m1/L1+m2/L1+m2/L2)*g -m2/L2*g; -m2/L2*g m2/L2*g]; % Stiffness matrix

    syms lambda

    L = solve(det(M*lambda^2+K)==0); % Determinate of Mass/Stiffness matrix

    wn = 0;
    wn(1,1) = min(abs(imag(L))); % Determine lower natural frequency
    wn(2,1) = max(imag(L));      % Determine higher natural frequency

    MAT1 = K+(wn(1)*i)^2*M;
    r1 = -MAT1(1,1)/MAT1(1,2);

    MAT2 = K+(wn(2)*i)^2*M;
    r2 = -MAT2(1,1)/MAT2(1,2);

    MS(:,1) = [1; r1];
    MS(:,2) = [1; r2];
end