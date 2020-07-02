%% Nuwan Dewapriya
%% 2019/04/09
%% This code has three parts.
%% The first part computes the coordinates of carbon atoms in a graphene sheet with specified dimensions on lines 18 & 19.
%% The second part generates .data file for LAMMPS
%% The third part generates .in file for LAMMPS

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1
%% Setting dimensions of the sheet

length=50; % All lengths are in Angstome
width=50;
c_c_bl=1.396; %% c-c bond acording to the airebo potential

%%
unit_x=c_c_bl*2*(1+cos(pi/3)); 
unit_y=c_c_bl*2*cos(pi/6);

u_x=floor(width/unit_x);
u_y=floor(length/unit_y);

%% Positions of the first 4 atoms

x(1,1)=c_c_bl*cos(pi/3);
x(1,2)=x(1,1)+c_c_bl;
x(2,1)=0;
x(2,2)=c_c_bl*(1+2*cos(pi/3));

y(1,1)=0;
y(1,2)=0;
y(2,1)=c_c_bl*cos(pi/6);
y(2,2)=c_c_bl*cos(pi/6);

%% Repetition along the x direction

for i=1: u_x
    x(1,(i*2)+1)=x(1,(i-1)*2+1)+c_c_bl*2*(1+cos(pi/3));
    x(1,(i*2)+2)=x(1,(i-1)*2+2)+c_c_bl*2*(1+cos(pi/3));
    x(2,(i*2)+1)=x(2,(i-1)*2+1)+c_c_bl*2*(1+cos(pi/3));
    x(2,(i*2)+2)=x(2,(i-1)*2+2)+c_c_bl*2*(1+cos(pi/3));
    
    y(1,(i*2)+1)=y(1,(i-1)*2+1);
    y(1,(i*2)+2)=y(1,(i-1)*2+2);
    y(2,(i*2)+1)=y(2,(i-1)*2+1);
    y(2,(i*2)+2)=y(2,(i-1)*2+2);
end

%% Repetition along the y direction

for i=1:u_y

    x((i*2)+1,:)=x((i-1)*2+1,:);
    x((i*2)+2,:)=x((i-1)*2+2,:);

    y((i*2)+1,:)=y((i-1)*2+1,:)+c_c_bl*cos(pi/6)*2;
    y((i*2)+2,:)=y((i-1)*2+2,:)+c_c_bl*cos(pi/6)*2;

end

Num_of_atoms = size(x,1)*size(x,2)

          
%% Getting x and y coordinates in to 2 columns

x=reshape(x,Num_of_atoms,1);
y=reshape(y,Num_of_atoms,1);

coord=zeros(Num_of_atoms,3);
coord(:,1)=x;
coord(:,2)=y;

l_x = max(x);
l_y = max(y);


%% Getting dimensions of the simulation box.

x_max=max(coord(:,1))+c_c_bl/2;
x_min=min(coord(:,1))-c_c_bl/2;
y_max=max(coord(:,2))+c_c_bl*cos(pi/6)/2;
y_min=min(coord(:,2))-c_c_bl*cos(pi/6)/2;

b_x =[x_min, x_min, x_max, x_max, x_min; y_min, y_max, y_max, y_min, y_min]; % Boundaries of the simulation box

Lx = x_max-x_min
Ly = y_max-y_min

%% Plotting the graphene sheet

plot(coord(:,1),coord(:,2),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k') % change this in an appropriate way for zig-zag and arm chair
axis equal  
hold on
plot(b_x(1,:),b_x(2,:),'--','Color','red','LineWidth',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2
%% Preparing the data file for LAMMPS

fid=fopen('grap.data','w');

fprintf(fid,'uniaxial tensile test of graphene\n');
fprintf(fid,'\n');
fprintf(fid,'%d atoms \n',size(coord,1));
fprintf(fid,'\n');
fprintf(fid,'%d atom types \n',1);
fprintf(fid,'\n');
fprintf(fid,'#simulation box \n');
fprintf(fid,'%f %f xlo xhi\n',x_min, x_max);
fprintf(fid,'%f %f ylo yhi\n',y_min, y_max);
fprintf(fid,'%f %f zlo zhi\n',-10.0, 10.0);
fprintf(fid,'\n');
fprintf(fid,'Masses\n');
fprintf(fid,'\n');
fprintf(fid,'%d %f \n',1, 12.010000);
fprintf(fid,'\n');
fprintf(fid,'Atoms\n');
fprintf(fid,'\n');

%% Defining atoms 

number=size(coord,1);
for i=1:number

    fprintf(fid,'%d 1 %f %f %f \n',i,coord(i,1),coord(i,2),rand/10); 
end

fclose(fid);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3
%% Preparing the simulation file for LAMMPS

fid=fopen('grap.in','w');

fprintf(fid,'#uniaxial tensile test of graphene\n');
fprintf(fid,'\n');

fprintf(fid,'##---------------INITIALIZATION-------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'units          metal\n');
fprintf(fid,'dimension 	    3 \n');

fprintf(fid,'boundary       p p f\n');

fprintf(fid,'atom_style 	atomic\n');
fprintf(fid,'newton 		on\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'##---------------ATOM DEFINITION------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'read_data 	grap.data\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'##---------------FORCE FIELDS---------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'pair_style 	airebo 3.0\n');
fprintf(fid,'pair_coeff     * * CH.airebo C\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'##---------------SETTINGS-------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'timestep 	0.0005\n');
fprintf(fid,'variable   ts equal 0.0005\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'##---------------COMPUTES-------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'compute 	1 all stress/atom NULL\n'); % for information, please check  http://lammps.sandia.gov/doc/compute_stress_atom.html
fprintf(fid,'compute    2 all reduce sum c_1[1] c_1[2]\n'); 
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'variable   Lx equal lx\n'); % See https://lammps.sandia.gov/doc/variable.html
fprintf(fid,'variable   Ly equal ly\n');
fprintf(fid,'variable   Lz equal lz\n');
fprintf(fid,'variable   Vol equal vol\n');

fprintf(fid,'variable   thickn equal 3.4\n'); %interlayer spacing is 3.4 A

fprintf(fid,'fix		1 all npt temp 300 300 0.05 x 0 0 0.5 y 0 0 0.5\n');

fprintf(fid,'thermo 	2000\n');

fprintf(fid,'##---------------RELAXATION--------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'run            50000\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'##---------------DEFORMATION--------------------------------------\n');
fprintf(fid,'unfix              1\n');
fprintf(fid,'reset_timestep     0\n');

fprintf(fid,'fix		1 all npt temp 300 300 0.05 x 0 0 0.5\n');

fprintf(fid,'fix        2 all ave/time 1 100 100 c_2[1] c_2[2]\n');% http://lammps.sandia.gov/doc/fix_ave_atom.html
fprintf(fid,'fix        3 all ave/time 1 100 100 v_Lx v_Ly v_Lz v_Vol\n');


fprintf(fid,'variable   srate equal 1.0e9\n');
fprintf(fid,'variable   srate1 equal "v_srate / 1.0e12"\n');


fprintf(fid,'fix		4 all deform 1 y erate ${srate1} units box remap x\n');
    
fprintf(fid,'run            100\n'); % to activate fixes


fprintf(fid,'##---------------THERMO-OUTPUTS--------------------------------------\n');


fprintf(fid,'variable   CorVol equal f_3[4]*v_thickn/(f_3[3])\n'); %corrected volume ignoring large box size

fprintf(fid,'variable   ConvoFac equal 1/1.0e4\n'); 

fprintf(fid,'variable   sigmaxx equal f_2[1]*v_ConvoFac/v_CorVol\n');
fprintf(fid,'variable   sigmayy equal f_2[2]*v_ConvoFac/v_CorVol\n'); 

fprintf(fid,'variable   StrainPerTs equal v_srate1*v_ts\n'); %strain per time step
fprintf(fid,'variable   strain equal v_StrainPerTs*step\n'); %strain 

fprintf(fid,'thermo 	100\n');
fprintf(fid,'thermo_style custom step temp v_strain v_sigmaxx v_sigmayy pe ke lx ly vol \n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'##---------------DEFORMATION--------------------------------------\n');

fprintf(fid,'dump           1 all atom 5000 tensile_test.lammpstrj\n');
fprintf(fid,'run            500000\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fclose(fid);
