function [k11,k12,k21,k22] = constitutive(Nx,Ny)

k11 = ones(Nx,Ny);
k12 = zeros(Nx,Ny);
k21 = zeros(Nx,Ny);
k22 = ones(Nx,Ny);