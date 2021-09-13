% chapter-6_Question 8
clearvars
clc

rho=0; % density 7800 kg/m^3
G = 0.8e11; % in N/m^2

d = [0.04 0.04];   % diameters of shafts in m

    nele = 2;       % no. of elements
    connect = [1 1 2;   %% First Column is element number
               2 2 3];    % Second & Third Column are Nodes (in sequence)
               
           
     coord = [1 0.000;       % first column node second column x coodinate
              2 0.610;       %provided length convert into equivalent length
              3 0.620];
         
     %mass = [4   5   6   7];     % mass of disc at each node

     %Dia = [0.08   0.1   0.12    0.14];  % dia of disc at node
     
     Ip = [0.03 0.02 0.04]; % in kg-m^2
     
bc = 1;         %1: free-free ; 2:fixed-fixed ; 3:fixed-free ; 4:free-fixed


for i=1:nele+1
    %Ip(i) = mass(i)*Dia(i)^2/8;   % polar moment of inertia of each disc
end

M = zeros(nele+1,nele+1);   %initializing global mass and stiffness matrix
K = zeros(nele+1,nele+1);

for i = 1:nele              %loop for finding elemental mass &
    nd1 = connect(i,2);     %stiffness matrix and assembly
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
    l = x2-x1;
    J = pi*d(i)^4/32;
    kele(:,:,i) = G*J/l*[1 -1;-1 1];
    if(i<nele)
        mele(:,:,i) = [Ip(nd1) 0; 0 0] ;      %disc is included only in left side
    else
        mele(:,:,i) = [Ip(nd1) 0;0 Ip(nd2)] ;   %disc is considered in both side
    end
    %assembly
    vec = [nd1 , nd2];  %global dof vector for assembly
    for ii = 1:2
        for jj = 1:2
            K(vec(ii),vec(jj)) =K(vec(ii),vec(jj))+kele(ii,jj,i);
            M(vec(ii),vec(jj)) =M(vec(ii),vec(jj))+mele(ii,jj,i);
        end
    end
end

% imposing boundary condition

if(bc==1)       %free-free
    Kred = K;
    Mred = M;
elseif(bc==2)      %fixed-fixed
    Kred = K(2:nele,2:nele);
    Mred = M(2:nele,2:nele);
elseif(bc==3)       %fixed-free
    Kred = K(2:nele+1,2:nele+1);
    Mred = M(2:nele+1,2:nele+1);
else                %free-fixed
    Kred = K(1:nele,1:nele);
    Mred = M(1:nele,1:nele);
end

D = Mred\Kred;
[eig_vec,eig_val] = eig(D);


for i=1:nele+1
    wnf(i) = sqrt(eig_val(i,i));      %critical speed
    mode(:,i) = eig_vec(:,i)/max(abs(eig_vec(:,i)));   %normalised eigen vector
end

for p=1:length(wnf)-1
  for i =1:length(wnf)-1      %arranging in asceending order
    if wnf(i)>wnf(i+1)
        temp_wnf = wnf(i);
        wnf(i) = wnf(i+1);
        wnf(i+1) = temp_wnf;
        temp_mode(:,1) = mode(:,i);
        mode(:,i) = mode(:,i+1);
        mode(:,i+1) = temp_mode(:,1);
    end
  end
end


%post processing

disp('Solution is printed to a text file "chapter-6_Question 8_output_FEM.txt"');
disp('mode shape is saved as "Mode_shape_chapter-6_Question 8_FEM.fig"');
fid = fopen('Output_chapter-6_Question 8_FEM.txt','w');
fprintf(fid,'Finite Element method\n\n');
if rho==0
    fprintf(fid,'shaft with multiple disc \n Neglecting shaft mass\n\n\n');
else
    fprintf(fid,'shaft with multiple disc \n Considering shaft mass\n\n\n');
end    
fprintf(fid,'Number of elements = %d\n',nele);
fprintf(fid,'Density of Shaft = %d\n',rho);

if(nele <= 6)
for i=1:nele
    fprintf(fid,'Element-%d\n',i);
    fprintf(fid,'----------\n\n');
    fprintf(fid,'Mass Matrix [M]%d\n',i);
    for ii = 1:2
        for jj=1:2
            fprintf(fid,'%.2f\t',mele(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Stiffness matrix [K]%d\n',i);
    for ii = 1:2
        for jj=1:2
            fprintf(fid,'%.2f\t',kele(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end


fprintf(fid,'Global Mass matrix \n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf(fid,'%.4f \t',M(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');
fprintf(fid,'Global Stiffness matrix \n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf(fid,'%1.2e \t',K(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
if(bc==1)
    fprintf(fid, '"Free-Free" Boundary condition\n');
elseif(bc==2)      
   fprintf(fid, 'Fixed-Fixed Boundary condition\n');
elseif(bc==3)      
   fprintf(fid, 'Fixed-Free Boundary condition\n');
else                
    fprintf(fid, 'Free-Fixed Boundary condition\n');
end
fprintf(fid,'\n');
fprintf(fid,'Natural frequencies:\n');
fprintf(fid,'%.3f \n',wnf);
fprintf(fid,'\n');
fprintf(fid,'Eigen vector matrix \n');
for i=1:length(wnf)
    for j= 1:length(wnf)
        fprintf(fid,'%.4f \t',eig_vec(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');
fprintf(fid,'Normalised eigen vector matrix \n');
for i=1:length(wnf)
    for j= 1:length(wnf)
        fprintf(fid,'%.4f \t',mode(i,j));
    end
    fprintf(fid,'\n');
end

fclose(fid);
end


for i = 1:length(wnf)
    phi_z=[];
    x_z = [];
    for j = 1:nele
        nd1 = connect(j,2);  
        nd2 = connect(j,3);
        x1 = coord(nd1,2);
        x2 = coord(nd2,2);
        l = x2-x1;
        x = x1:l/10:x2;     %deviding each element into small element for smooth mode shape
        phi1 = mode(j,i);
        phi2 = mode(j+1,i);
        for k=1:length(x)
            z = x(k)-x1;
            N1 = 1-(z/l);
            N2 = z/l;
            phi(k)=[N1 N2]*[phi1; phi2];
        end
        phi_z = [phi_z,phi];
        x_z =[x_z,x];
    end
       
    mode_1(:,i)=phi_z; % /max(abs(phi_z));  % normalised
end



h = figure(1);
set(gcf, 'Position', get(0,'Screensize'));

if(size(mode_1,2) < 4)
    nmode=size(mode_1,2);
else
    nmode=4;
end


for i = 1: nmode %length(wnf)
    if(i ==1)
        plot(x_z,mode_1(:,i), '-k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==2)
         plot(x_z,mode_1(:,i), ':k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==3)
         plot(x_z,mode_1(:,i),'-.k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==4)
         plot(x_z,mode_1(:,i), '--k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    else
         plot(x_z,mode_1(:,i), 'DisplayName',['\omega=',num2str(wnf(i))]);
    end
         
    %hold on;
    %plot(x_z,mode_1(:,i),'Displayname',['\omega=',num2str(wnf(i))]);
    hold on;
end
grid on;
xlabel('shaft length(m)','fontsize',16);
ylabel('relative amplitude','fontsize',16);
%title('Mode shape','fontsize',20);
ylim([-1.2 1.2]);
legend('show');
saveas(h,'chapter-6_Question 8_mode_shape_FEM','png');