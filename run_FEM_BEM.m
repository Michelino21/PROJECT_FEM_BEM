clear all; close all; clc;
addpath(genpath('src'));

% === PARAMETRI ===
params.Lx               = 0.05;
params.Ly               = 0.05;
params.delx             = 0.003;
params.dely             = params.delx;
params.margin_BEM_coeff = 3;
params.e0               = 8.85e-12;
params.er               = 1;
params.V1               = +0.5;
params.V2               = -0.5;

% === SALVATAGGIO ===
params.SALVATAGGIO    = 1;
params.ERROR_ANALYSIS = 0;
params.margini_list = [params.margin_BEM_coeff, 50];

% === CHIAMA LA FUNZIONE ===
FEM_BEM(params);