% Computing Coefficients
collision_freq = 0;
Chex = dt/(mu_0*dy);    Chey = -dt/(mu_0*dx);
Cpp = (1-collision_freq*dt/2)/(1+collision_freq*dt/2);
Cpe = eps_0*omega_p^2*dt/(1+collision_freq*dt/2);
Cexh = dt/(eps_0*dy);   Cexp = -dt/eps_0;
Ceyh = -dt/(eps_0*dx);   Ceyp = -dt/eps_0;

BandSq(length(kx),30)=0;    % array to store band sequances
if unit_cell == square
    % computing band structure
    for kInd = 1:length(kx)
        % resetting probes
        ProbValues = zeros(number_of_time_steps, n_Probes);
        % resetting fields and currents
        Ex(nx,ny+1) = 0; Ptx = Ex*0; %#ok<SAGROW>
        Ey(nx+1,ny) = 0; Pty = Ey*0; %#ok<SAGROW>
        Hz(nx,ny) = 0; %#ok<SAGROW>

        for time_step = 1:number_of_time_steps
            % update magnetic fields
            Hz = Hz + Chex*(Ex(:,2:end)-Ex(:,1:end-1)) + Chey*(Ey(2:end,:)-Ey(1:end-1,:));
            % record magnetic fields at probe positions
            for ProbInd = 1:n_Probes
                ProbValues(time_step,ProbInd) = Hz(ProbPos_x(ProbInd),ProbPos_y(ProbInd));
            end
            % update polarization currents
            Ptx(~Rods_x) = Cpp*Ptx(~Rods_x) + Cpe*Ex(~Rods_x);     %Ptx(Rods_x) = 0;
            Pty(~Rods_y) = Cpp*Pty(~Rods_y) + Cpe*Ey(~Rods_y);     %Pty(Rods_y) = 0;
            % update electric field Ex -- applying bloch Boundary condition
            Hzex = [Hz Hz(:,1)*exp(1i*a*ky(kInd))];
            Ex(:,2:end) = Ex(:,2:end) + Cexh*(Hzex(:,2:end)-Hzex(:,1:end-1)) + Cexp*Ptx(:,2:end);
            Ex(:,1) = Ex(:,end)*exp(-1i*a*ky(kInd));
            % Updating Ey
            Hzey = [Hz; Hz(1,:)*exp(1i*a*kx(kInd))];
            Ey(2:end,:) = Ey(2:end,:) + Ceyh*(Hzey(2:end,:)-Hzey(1:end-1,:)) + Ceyp*Pty(2:end,:);
            Ey(1,:) = Ey(end,:)*exp(-1i*a*kx(kInd));
            % Updating Sources
            for SourceInd = 1:n_Sources
                Ex(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) = Ex(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) + Cexp*waveform(time_step); %#ok<SAGROW>
            end
            Ex(Rods_x)=0; Ey(Rods_y)=0;  %#ok<SAGROW>
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
        Ex(nx,ny+1) = 0; Ptx = Ex*0; %#ok<SAGROW>
        Ey(nx+1,ny) = 0; Pty = Ey*0; %#ok<SAGROW>
        Hz(nx,ny) = 0; %#ok<SAGROW>

        for time_step = 1:number_of_time_steps
            % update magnetic fields
            Hz = Hz + Chex*(Ex(:,2:end)-Ex(:,1:end-1)) + Chey*(Ey(2:end,:)-Ey(1:end-1,:));
            % record magnetic fields at probe positions
            for ProbInd = 1:n_Probes
                ProbValues(time_step,ProbInd) = Hz(ProbPos_x(ProbInd),ProbPos_y(ProbInd));
            end
            % update polarization currents
            Ptx(~Rods_x) = Cpp*Ptx(~Rods_x) + Cpe*Ex(~Rods_x);
            Pty(~Rods_y) = Cpp*Pty(~Rods_y) + Cpe*Ey(~Rods_y);
            % update electric field Ex -- applying bloch Boundary condition
            Hzex=[Hz [Hz(end/2+1:end,1)*exp(-1i*a*(kx(kInd)-sqrt(3)*ky(kInd))/2); Hz(1:end/2,1)*exp(1i*a*(kx(kInd)+sqrt(3)*ky(kInd))/2)]];
            Ex(:,2:end) = Ex(:,2:end) + Cexh*(Hzex(:,2:end)-Hzex(:,1:end-1)) + Cexp*Ptx(:,2:end);
            Ex(:,1)=[Ex(end/2+1:end,end)*exp(-1i*a*(kx(kInd)+sqrt(3)*ky(kInd))/2); Ex(1:end/2,end)*exp(1i*a*(kx(kInd)-sqrt(3)*ky(kInd))/2)];
            % Updating Ey
            Hzey = [Hz; Hz(1,:)*exp(1i*a*kx(kInd))];
            Ey(2:end,:) = Ey(2:end,:) + Ceyh*(Hzey(2:end,:)-Hzey(1:end-1,:)) + Ceyp*Pty(2:end,:);
            Ey(1,:) = Ey(end,:)*exp(-1i*a*kx(kInd));
            % Updating Sources
            for SourceInd = 1:n_Sources
                Ex(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) = Ex(SourcePos_x(SourceInd),SourcePos_y(SourceInd)) + Cexp*waveform(time_step); %#ok<SAGROW>
            end
            Ex(Rods_x)=0; Ey(Rods_y)=0;  %#ok<SAGROW>
        end
        PostProcess;
    end
end
% displaying results