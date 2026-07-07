% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) MOP2(x);  % Cost Function

nVar = 5; % Number of Decision Variables

VarSize = [1 nVar]; % Size of Decision Variables Matrix

VarMin = -1*ones(VarSize); % Lower Bound of Variables
VarMax = 1*ones(VarSize); % Upper Bound of Variables

%%

[x, fval] = nsga3(CostFunction, VarMin, VarMax, []);
