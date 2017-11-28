if unit_cell == square
    % Vertices of brillion zone
    xVert=[0 pi/a pi/a];      % x coordinates of the vertices of the loop closing irreducible brillion zone in k space
    yVert=[0 0    pi/a];
end

if unit_cell == triangular
    % Vertices irreducible brillion zone
    xVert=[0 0                 2*pi/(3*a)];      % x coordinates of the vertices of the loop closing irreducible brillion zone in k space
    yVert=[0 2*pi/(sqrt(3)*a)  2*pi/(sqrt(3)*a)];
end

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
