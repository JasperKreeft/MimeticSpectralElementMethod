function [F] = force_zeroform(X,Y,FunctionType)

F = exact_solution(X,Y,FunctionType,'force');