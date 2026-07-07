%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Definition

CostFunction=@(x) ZDT(x);  % Cost Function

nVar=3;             % Number of Decision Variables

VarSize=[nVar 1];   % Decision Variables Matrix Size

VarMin = 0*ones(VarSize);         % Decision Variables Lower Bound
VarMax = 1*ones(VarSize);         % Decision Variables Upper Bound

[pos_pareto, cost_pareto, pos_all, cost_all, pop] = moead(CostFunction, VarMin, VarMax);

