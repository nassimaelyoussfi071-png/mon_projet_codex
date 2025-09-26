%COMPUTE_UIO_GAINS_SCRIPT Solve the LMIs of Theorem 1 and display the gains.
%   This script reproduces the workflow of the original snippet while
%   addressing the issues that prevented valid solutions.  It loads the
%   Takagi-Sugeno model from ETAPE1.mat, formulates one LMI per subsystem,
%   solves the feasibility problem and prints the resulting matrices P and
%   K1.
%
%   Run the script from MATLAB or Octave with the `ETAPE1.mat` file in the
%   current directory.

%% Load data
S = load('ETAPE1.mat', 'A1', 'Hi', 'Ti'); %#ok<NASGU>
A1 = S.A1;

%% Electrical parameters that define matrix C
L_s = 0.1094;
L_r = 0.1071;
L_m = 0.1054;
sigma = 1 - (L_m^2 / (L_s * L_r));
ig = 1200;

C = [
    1/(sigma*L_s), 0, -L_m/(sigma*L_s*L_r), 0, 0, 0;
    0, 1/(sigma*L_s), 0, -L_m/(sigma*L_s*L_r), 0, 0;
    0, 0, 0, 0, 1/ig, 0;
    0, 0, 0, 0, 0, 1/ig];

n = size(A1{1}, 1);
p = size(C, 1);

if size(C, 2) ~= n
    error('Dimension mismatch: the columns of C must equal the state dimension.');
end

%% Define LMI variables
setlmis([]);
P = lmivar(1, [n 1]);
X = cell(1, numel(A1));
for i = 1:numel(A1)
    X{i} = lmivar(2, [n p]);
end

%% Enforce P > 0
lmiterm([-1 1 1 P], 1, 1);

%% One LMI per subsystem (equation 21)
for i = 1:numel(A1)
    lmiterm([i 1 1 P], 1, A1{i}, 's');
    lmiterm([i 1 1 X{i}], 1, -C, 's');
end

%% Solve the LMI feasibility problem
LMISYS = getlmis;
[tmin, xopt] = feasp(LMISYS);

if tmin > 0
    error('The LMIs are infeasible (tmin = %g). Check the problem data.', tmin);
end

%% Recover the matrices
P = dec2mat(LMISYS, xopt, P);
K1 = cell(1, numel(A1));
for i = 1:numel(A1)
    X{i} = dec2mat(LMISYS, xopt, X{i});
    K1{i} = P \ X{i};
end

%% Display the results
fprintf('P =\n');
disp(P);
for i = 1:numel(A1)
    fprintf('K1{%d} =\n', i);
    disp(K1{i});
end

%% Optionally save the results
save('ETAPE2.mat', 'K1', 'P');
