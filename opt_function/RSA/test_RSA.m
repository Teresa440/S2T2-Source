clear all
close all
addpath E:\Users\Utente\OneDrive\POLITO\TESI\CODICE\SROC\opt_function\DTLZ_TEST_FUNCTIONS\Octave-Matlab
% addpath D:\Utenti\david\OneDrive\POLITO\TESI\CODICE\SROC\opt_function\DTLZ_TEST_FUNCTIONS\Octave-Matlab % portatile 

% M = 4; % number of objectives
% opt_function = @(x) dtlz1(x',M); % test problem
% k = 5; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific

% M = 4; % number of objectives
% opt_function = @(x) dtlz2(x',M); % test problem
% k = 10; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific

M = 4; % number of objectives
opt_function = @(x) dtlz3(x',M); % test problem
k = 10; % problem specific
nvars = (M - 1) + k; % number of variables, problem specific

% M = 4; % number of objectives
% opt_function = @(x) dtlz4(x',M); % test problem
% k = 10; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific

% M = 4; % number of objectives
% opt_function = @(x) dtlz5(x',M); % test problem
% k = 10; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific

% M = 4; % number of objectives
% opt_function = @(x) dtlz6(x',M); % test problem
% k = 10; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific
 
% M = 4; % number of objectives
% opt_function = @(x) dtlz7(x',M); % test problem
% k = 20; % problem specific
% nvars = (M - 1) + k; % number of variables, problem specific

% lower and upper bounds
lb = zeros(1,nvars);
ub = ones(1,nvars);

% tic
% options = optimoptions("gamultiobj","FunctionTolerance",1e-15);
% [x, fval] = gamultiobj(opt_function,nvars,[],[],[],[],lb,ub,[],[],options);
% toc
% figure
% % plot(fval(:,1),fval(:,2),'ok')
% plot3(fval(:,1),fval(:,2),fval(:,3),'ok')

tic
[x, fval] = RSA(opt_function,lb,ub);
toc
% figure
% plot(fval(:,1),fval(:,2),'ok')
% plot3(fval(:,1),fval(:,2),fval(:,3),'ok')
