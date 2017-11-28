clear; clc; %close all;

number_of_time_steps=2^14;  no_of_k=30;
courant_factor = 0.99;
max_source_freq = 8e10;
n_grid_pt_per_lattice_const=100;
% Constants
eps_0 = 8.854187817e-12; mu_0  = 4*pi*1e-7; c = 1/sqrt(mu_0*eps_0);
% number_of_frequencies = 2^10;   % for post_pocess without fft
syms triangular square TE TM PPC MPC PMPC Magnetized_PMPC % PPC type 1 == plasma rods in vaccum background
%triangular=0; square=1; TE=0; TM=1; PPC=0; MPC=1; PMPC=2; Magnetized_PMPC=3; % PPC type 1 == plasma rods in vaccum background

unit_cell = triangular;
mode = TM;
Case = PMPC; ppcType = 1; %type 1 == plasma rods in vaccum background

    a=1;
    r=0.1; %sqrt(.5*a^2/pi);%0.2;0.2e-2; %
    omega_p_Array = (0.1:0.1:2)*pi*c/a;
    collision_freq = 0;
%% ---------------------------------
r_by_a = [0:0.05:0.15 0.2:0.025:0.5];
BandSq_X(length(r_by_a),length(omega_p_Array),30)=0;
BandSq_Gamma=BandSq_X*0;
BandSq_J=BandSq_X*0;
% figure(2);
    for plasmaFreqInd = 1:length(omega_p_Array)
        omega_p = omega_p_Array(plasmaFreqInd);
        omega_str='000000000';
        omega_str(1:length(num2str(omega_p*1e-9)))=num2str(omega_p*1e-9);
        omega_folder_string = [char(unit_cell) '\' char(mode) '\plasmaFreq_' omega_str '_GHz']; omega_folder_string(omega_folder_string=='.') = '_';
        mkdir(omega_folder_string);
    end
compute_k_space_parameters;
plasmaFreqInd = 1; rInd=13;
    
if exist('status2.mat','file'), load status2; end
plasmaFreqInd0=plasmaFreqInd;
rInd0 = rInd;
for rInd=rInd0:length(r_by_a)
    r = r_by_a(rInd);
    r_str='00000';
    r_str(1:length(num2str(r)))=num2str(r); r_str(r_str=='.') = '_';
    r_folder_string = [char(unit_cell) '\' char(mode) '\r_by_a_' r_str]; 
    mkdir(r_folder_string);
    initialize_parameters;
    for plasmaFreqInd = plasmaFreqInd0:length(omega_p_Array)
        close all
        omega_p = omega_p_Array(plasmaFreqInd);
        omega_str='000000000';
        omega_str(1:length(num2str(omega_p*1e-9)))=num2str(omega_p*1e-9);   omega_str(omega_str=='.') = '_';
        omega_folder_string = [char(unit_cell) '\' char(mode) '\plasmaFreq_' omega_str '_GHz']; 
        %         mkdir(omega_folder_string);
        if mode==TE
            disp('mode = TE')
            if Case == PPC
                time_march_PPC_TE
            elseif Case == PMPC
                time_march_PMPC_TE
            end
        else
            disp('mode = TM')
            time_march_PMPC_TM
        end
        BandSqOriginal = BandSq;
        save([omega_folder_string, '\BandSqOriginal_r_by_a_', r_str],'BandSqOriginal','-MAT');
        save([r_folder_string, '\BandSqOriginal_PlasmaFreq_', omega_str, '_GHz'],'BandSqOriginal','-MAT');
        %% --- Plotting results
        BandSq = BandSqOriginal;
        % process_BandSq_2
        figure()
        if Case == PPC
            plot(0:length(kx)-1,BandSq*c/a*1e-9,'x','LineWidth',1,'color',[0 0 0]);
            axis([0 length(kx)-1 0 80]); axis1=axis;
        else
            plot(0:length(kx)-1,BandSq,'x','LineWidth',1,'color',[0 0 0]);
            axis([0 length(kx)-1 0 4]); axis1=axis;
        end    
        xlabel('Wave vector','FontSize',20);
        ylabel('Normalised Frequency','FontSize',20);
        if mode == TE
            title('Tranverse Electric (TEz) Photonic Band Structure')
        else
            title('Trasverse Magnetic (TMz) Photonic Band Structure')
        end
        %line(0:length(kx)-1,BandSq);
        line([vertPos(1)-1 vertPos(1)-1],[0 axis1(4)])
        line([vertPos(2)-1 vertPos(2)-1],[0 axis1(4)])
        %%
        BandSq_Gamma(rInd,plasmaFreqInd,:)=BandSq(end,:); % the first frequency in the lowest band is generally missing
        BandSq_J(rInd,plasmaFreqInd,:)=BandSq(vertPos(1),:);
        BandSq_X(rInd,plasmaFreqInd,:)=BandSq(vertPos(2),:);
        saveas(gcf,[omega_folder_string '\r_by_a_' r_str '.jpg'],'jpg')
        saveas(gcf,[r_folder_string '\PlasmaFreq_' omega_str '_GHz.jpg'],'jpg')
        save status2
    end
plasmaFreqInd0=1;
end