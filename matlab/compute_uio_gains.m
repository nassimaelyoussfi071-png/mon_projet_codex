function [P,K1,X,tmin] = compute_uio_gains(dataFile)
%COMPUTE_UIO_GAINS Solve the LMIs of Theorem 1 to obtain UIO gains.
%   [P,K1,X,tmin] = COMPUTE_UIO_GAINS(dataFile) loads the TS system
%   matrices stored in the MAT-file specified by dataFile (default
%   'ETAPE1.mat') and solves the LMIs in equation (21) of Theorem 1.
%   The function returns the positive definite matrix P, the cell array
%   of gain matrices K1 and the auxiliary decision variables X.  The
%   scalar tmin is the feasibility margin returned by FEASP.
%
%   The implementation fixes two issues that frequently lead to invalid
%   solutions when using the original script:
%     1) A dedicated LMI is created for each subsystem i instead of
%        accumulating all terms in a single LMI.  This reflects the
%        requirements of Theorem 1.
%     2) Dimension checks and the use of the backslash operator avoid
%        numerical instabilities such as calling INV(P).
%
%   Example:
%       [P,K1] = compute_uio_gains('ETAPE1.mat');
%
%   See also SETLMIS, LMIVAR, LMITERM, FEASP, DEC2MAT.

if nargin == 0
    dataFile = 'ETAPE1.mat';
end

% Load the data required by Theorem 1.
S = load(dataFile, 'A1', 'Hi', 'Ti'); %#ok<NASGU>
A1 = S.A1;

% Electrical parameters that define the matrix C.
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

setlmis([]);
P = lmivar(1, [n 1]);
X = cell(1, numel(A1));
for i = 1:numel(A1)
    X{i} = lmivar(2, [n p]);
end

% Enforce P > 0.
lmiterm([-1 1 1 P], 1, 1);

% One LMI per subsystem, as required by (21).
for i = 1:numel(A1)
    lmiterm([i 1 1 P], 1, A1{i}, 's');
    lmiterm([i 1 1 X{i}], 1, -C, 's');
end

LMISYS = getlmis;
[tmin, xopt] = feasp(LMISYS);

if tmin > 0
    error('The LMIs are infeasible (tmin = %g). Check the problem data.', tmin);
end

P = dec2mat(LMISYS, xopt, P);
for i = 1:numel(A1)
    X{i} = dec2mat(LMISYS, xopt, X{i});
end

K1 = cell(1, numel(A1));
for i = 1:numel(A1)
    K1{i} = P \ X{i};
end

end
