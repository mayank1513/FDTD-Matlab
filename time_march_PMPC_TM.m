% Computing Coefficients
collision_freq = 0;
Cehx = dt/(eps_0*dy);    Cehy = -dt/(eps_0*dx); 
Cezp = -dt/eps_0;
Cpp = (1-collision_freq*dt/2)/(1+collision_freq*dt/2);
Cpe = eps_0*omega_p^2*dt/(1+collision_freq*dt/2);
Chxe = dt/(mu_0*dy);   
Chye = -dt/(mu_0*dx);   %Ceyp = -dt/mu_0;

BandSq(length(kx),30)=0;    % array to store band sequances
if unit_cell == square
    % computing band structure
    for kInd = 1:length(kx)
        % resetting probes
        ProbValues = zeros(number_of_time_steps, n_Probes);
        % resetting fields and currents
        Hx(nx,ny+1) = 0; %#ok<SAGROW>
        Hy(nx+1,ny) = 0; %#ok<SAGROW>
        Ez(nx,ny) = 0; Ptz = Ez*0; %#ok<SAGROW>

        for time_step = 1:number_of_time_steps
            % update magnetic fields
            Ezhx = [Ez Ez(:,1)*exp(1i*a*ky(kInd))];
            Hx(:,2:end) = Hx(:,2:end) + Chxe*(Ezhx(:,2:end)-Ezhx(:,1:end-1));% + Cexp*Ptz(:,2:end);
            Hx(:,1) = Hx(:,end)*exp(-1i*a*ky(kInd));
            % Updating Hy
            Ezhy = [Ez; Ez(1,:)*exp(1i*a*kx(kInd))];
            Hy(2:end,:) = Hy(2:end,:) + Chye*(Ezhy(2:end,:)-Ezhy(1:end-1,:));% + Ceyp*Pty(2:end,:);
            Hy(1,:) = Hy(end,:)*exp(-1i*a*kx(kInd));
            % update polarization currents
            Ptz(~Rods_z) = Cpp*Ptz(~Rods_z) + Cpe*Ez(~Rods_z);     %Ptz(Rods_x) = 0;
            % update electric field Ez -- applying bloch Boundary condition
            Ez = Ez + Cehx*(Hx(:,2:end)-Hx(:,1:end-1)) + Cehy*(Hy(2:end,:)-Hy(1:end-1,:)) + Cezp*Ptz;
            % record magnetic fields at probe positions
            for ProbInd = 1:n_Probes
                ProbValues(time_step,ProbInd) = Ez(ProbPos_x(ProbInd),ProbPos_y(ProbInd));
            end
            % Updating Sources
            for SourceInd = 1:n_Sources
                Ez(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) = Ez(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) + Cezp*waveform(time_step); %%#ok<SAGROW>
            end
            Ez(Rods_z)=0;
        end
        PostProcess;
    end
end


if unit_cell == triangular
    % computing band structure
    for kInd = 1:length(kx)
        % resetting probes
        ProbValues = zeros(number_of_time_steps, n_Probes);
        % resetting fields and currents
        Hx(nx,ny+1) = 0; %#ok<SAGROW>
        Hy(nx+1,ny) = 0; %#ok<SAGROW>
        Ez(nx,ny) = 0; Ptz = Ez*0; %#ok<SAGROW>

        for time_step = 1:number_of_time_steps
            % update magnetic fields
            % update electric field Hx -- applying bloch Boundary condition
            Ezhx=[Ez [Ez(end/2+1:end,1)*exp(-1i*a*(kx(kInd)-sqrt(3)*ky(kInd))/2); Ez(1:end/2,1)*exp(1i*a*(kx(kInd)+sqrt(3)*ky(kInd))/2)]];
            Hx(:,2:end) = Hx(:,2:end) + Chxe*(Ezhx(:,2:end)-Ezhx(:,1:end-1));% + Cezp*Ptz(:,2:end);
            Hx(:,1)=[Hx(end/2+1:end,end)*exp(-1i*a*(kx(kInd)+sqrt(3)*ky(kInd))/2); Hx(1:end/2,end)*exp(1i*a*(kx(kInd)-sqrt(3)*ky(kInd))/2)];
            % Updating Hy
            Ezhy = [Ez; Ez(1,:)*exp(1i*a*kx(kInd))];
            Hy(2:end,:) = Hy(2:end,:) + Chye*(Ezhy(2:end,:)-Ezhy(1:end-1,:));% + Ceyp*Pty(2:end,:);
            Hy(1,:) = Hy(end,:)*exp(-1i*a*kx(kInd));
            % update polarization currents
            Ptz(~Rods_z) = Cpp*Ptz(~Rods_z) + Cpe*Ez(~Rods_z);
            % Updating Ez
            Ez = Ez + Cehx*(Hx(:,2:end)-Hx(:,1:end-1)) + Cehy*(Hy(2:end,:)-Hy(1:end-1,:)) + Cezp*Ptz;
            % record magnetic fields at probe positions
            for ProbInd = 1:n_Probes
                ProbValues(time_step,ProbInd) = Ez(ProbPos_x(ProbInd),ProbPos_y(ProbInd));
            end
            % Updating Sources
            for SourceInd = 1:n_Sources
                Ez(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) = Ez(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) + Cezp*waveform(time_step); %%#ok<SAGROW>
            end
            Ez(Rods_z)=0; 
        end
        PostProcess;
    end
end
% displaying results