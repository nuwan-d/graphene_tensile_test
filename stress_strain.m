%% Nuwan Dewapriya
%% 2019/04/10
%% This code extracts data from the LAMMPS output file log.lammps and plot the stress-strain curve.

clear all
close all
clc

%% Extracting stress-strain data from the log.lammps file

[fid] = fopen('log.lammps');

[Stress,count] = fscanf(fid, '%*d %*f %f %*f %f %*f %*f %*f %*f %*f ',[2,inf]);%% etract only Step and E_pair     

Stress = Stress'; %

%% Plotting

figure
plot(Stress(:,1), Stress(:,2),'-or','MarkerSize',2)
xlabel('Strain','fontsize',12)
ylabel('Stress (GPa)','fontsize',12)
grid on
set(gca,'LineWidth',1,'Fontsize',12)
axis square
axis([0 0.25 0 110])
grid off
