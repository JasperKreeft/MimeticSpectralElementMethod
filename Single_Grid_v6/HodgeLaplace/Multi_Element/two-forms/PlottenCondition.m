clear all
% close all
clc

load 'Hconv_N6_c0.3.mat'

HconvRange = 2.^(1:6);
h = 2./HconvRange;

rate = mean((log(ConditionNumber(end))-log(ConditionNumber(2:end-1)))./(log(h(end))-log(h(2:end-1))))

logCondRef = log(10000)+(-2.5)*(log(h)-log(h(1)));
CondRef = exp(logCondRef);

loglog(2./HconvRange,ConditionNumber,'-m')
hold on
loglog(2./HconvRange,CondRef,'-r')
xlabel('h')