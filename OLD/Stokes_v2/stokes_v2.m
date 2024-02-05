clear all
close all
clc

N=2;

[xiG,wG]=Gnodes(N);
[xiGLL,wGLL]=GLLnodes(N);
xiEG = [-1 xiG 1];

[Dp,Gd,NGp,Cd,Cp,Dd,Gp] = topology(N);

[Aq,Au] = hodge_d1p1(N,xiGLL,xiEG,wG,wGLL);

[Bv,Bu_bc,Bu_in] = som_d1p0(N,xiGLL,xiEG,wG,wGLL,NGp);

[Cl,Cv] = hodge_p1d1(N,xiGLL,xiEG,wG,wGLL);

