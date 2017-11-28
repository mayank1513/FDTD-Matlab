if unit_cell == square
    nx=n_grid_pt_per_lattice_const; dx=a/nx;
    ny=nx; dy=dx;
    
    Ind1 = repmat([1:nx]',1,ny+1); %#ok<NBRAK>
    Ind2 = repmat(1:ny+1,nx,1);
    Rods_x = ((Ind1-0.5-nx/2)*dx).^2+((Ind2-1-ny/2)*dy).^2<r^2;
    
    Ind1 = repmat([1:nx+1]',1,ny); %#ok<NBRAK>
    Ind2 = repmat(1:ny,nx+1,1);
    Rods_y = ((Ind1-1-nx/2)*dx).^2+((Ind2-0.5-ny/2)*dy).^2<r^2;
    
    Ind1 = repmat([1:nx]',1,ny); %#ok<NBRAK>
    Ind2 = repmat(1:ny,nx,1);
    Rods_z = ((Ind1-0.5-nx/2)*dx).^2+((Ind2-0.5-ny/2)*dy).^2<r^2;
    % source posiitons
    SourcePos_x = round([nx/12  nx/12     11*nx/12]);
    SourcePos_y = round([ny/12  11*ny/12  11*ny/12]);
    % prob posiitons
    ProbPos_x = round([nx/12 nx/12   nx/2   nx/2]);
    ProbPos_y = round([ny/2  4*ny/5  ny/12  11*ny/12]);
    
    % Vertices of brillion zone
    xVert=[0 pi/a pi/a];      % x coordinates of the vertices of the loop closing irreducible brillion zone in k space
    yVert=[0 0    pi/a];
end

if unit_cell == triangular
    nx=n_grid_pt_per_lattice_const; dx=a/nx;
    ny=round(sqrt(3)*nx/2);    dy=sqrt(3)*a/2/ny;  
    
    Ind1 = repmat([1:nx]',1,ny+1); %#ok<NBRAK>
    Ind2 = repmat(1:ny+1,nx,1);
    Rods_x = (((Ind1-0.5-nx/2)*dx).^2+((Ind2-1)*dy).^2<r^2)|(((Ind1-0.5-nx)*dx).^2+((Ind2-1-ny)*dy).^2<r^2) | (((Ind1-0.5)*dx).^2+((Ind2-1-ny)*dy).^2<r^2);
    
    Ind1 = repmat([1:nx+1]',1,ny); %#ok<NBRAK>
    Ind2 = repmat(1:ny,nx+1,1);
    Rods_y = (((Ind1-1-nx/2)*dx).^2+((Ind2-0.5)*dy).^2<r^2)|(((Ind1-1-nx)*dx).^2+((Ind2-0.5-ny)*dy).^2<r^2) | (((Ind1-1)*dx).^2+((Ind2-0.5-ny)*dy).^2<r^2);
    
    Ind1 = repmat([1:nx]',1,ny); %#ok<NBRAK>
    Ind2 = repmat(1:ny,nx,1);
    Rods_z = (((Ind1-0.5-nx/2)*dx).^2+((Ind2-0.5)*dy).^2<r^2)|(((Ind1-.5-nx)*dx).^2+((Ind2-0.5-ny)*dy).^2<r^2) | (((Ind1-.5)*dx).^2+((Ind2-0.5-ny)*dy).^2<r^2);
    % source posiitons
    SourcePos_x = round([nx/2              nx/12                11*nx/12]);
    SourcePos_y = round([nx/2+(ny-nx/2)/2  ny-nx/2-(ny-nx/2)/3  ny-nx/2-(ny-nx/2)/3]);
    % prob posiitons
    ProbPos_x = round([nx/2                nx/4  3*nx/4]);
    ProbPos_y = round([nx/2+5*(ny-nx/2)/6  ny/2  ny/2]);
    
    % Vertices irreducible brillion zone
    xVert=[0 0                 2*pi/(3*a)];      % x coordinates of the vertices of the loop closing irreducible brillion zone in k space
    yVert=[0 2*pi/(sqrt(3)*a)  2*pi/(sqrt(3)*a)];
end
%}
n_Probes = length(ProbPos_x);
n_Sources = length(SourcePos_x);
% Duration of a time step in seconds
dt = 1/(c*sqrt((1/dx^2)+(1/dy^2)));
dt = courant_factor*dt;
time = (1:number_of_time_steps)*dt;

% computing kx, ky arrays
%{
disp('Computing k-space parameters');
dK=(sqrt((xVert(2)-xVert(1))^2+(yVert(2)-yVert(1))^2)+sqrt((xVert(3)-xVert(2))^2+(yVert(3)-yVert(2))^2)+sqrt((xVert(3)-xVert(1))^2+(yVert(3)-yVert(1))^2))/no_of_k;

kx=[]; ky=[]; 
vertPos(3)=0;
for PathInd=1:3
    kx1=xVert(PathInd);                 ky1=yVert(PathInd);
    kx2=xVert(mod(PathInd,3)+1);        ky2=yVert(mod(PathInd,3)+1);
    if kx1~=kx2
        m=(ky2-ky1)/(kx2-kx1);
        dkx=dK/sqrt(m^2+1);     dky=m*dkx;
        if kx1<kx2
            kx_=kx1:dkx:(kx2-dkx);
            kx=[kx kx_]; %#ok<*AGROW>
        else
            kx_=kx1:-dkx:(kx2+dkx);
            kx=[kx kx_];
        end
        if ky1<ky2
            ky=[ky ky1:dky:(ky2-dky)];
        elseif ky1>ky2
            ky=[ky ky1:-dky:(ky2+dky)];
        else
            ky=[ky repmat(ky1,size(kx_))];
        end
    else
        if ky1<ky2
            ky_=ky1:dK:(ky2-dK);
        else
            ky_=ky1:-dK:(ky2+dK);
        end
        ky=[ky ky_];
        kx=[kx repmat(kx1,size(ky_))];
    end
    vertPos(PathInd)=length(kx)+1;
end
%}
% computing source waveform -- derivative gausian

tau = 1/(2*max_source_freq);
t_0 = 4.5 * tau;
waveform = -(sqrt(2*exp(1))/tau)*(time - t_0).*exp(-((time - t_0)/tau).^2);
