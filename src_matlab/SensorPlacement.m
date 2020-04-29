%%% © Stefano Maranò (2016)
%%%
%%% This script implements an algorithm for sensor placement presented in
%%% S. Maranò, D. Fäh, and Y. M. Lu,
%%% "Sensor placement for the analysis of seismic surface waves: Sources of
%%% error, design criterion and array design algorithms"
%%% Geophys. J. Int., 2014.
%%%
%%% Available online from GJI https://doi.org/10.1093/gji/ggt489
%%% 

clear all; close all;

run gurobi651/linux64/matlab/gurobi_setup           % TODO set the path to the actual gurobi_setup.m on your machine
if isempty(strfind(path, [pwd, '/utils']))          % adding ./utils folder to path
    path(path, [pwd, '/utils'])
end

%%%
%%% Some options
%%%

OUTPUT_DIR='output';        % Folder where output MAT files are saved

N_sensors_vec=[5, 7, 11];    % This is an array containing the number of sensors of interest
k_max=1;                    % Largest wavenumber of interest. Normalized wavenumber MUST BE k_max=1.
k_min_vec=k_max./[3, 5];    % Array with the smallest wavenumbers of interest eg: k_min_vec=[0.2, 0.25, 0.5]


% parameters defining the possible sensor locations
% these parameters can be overridden in the code.
% Reson for overriding, for example:
% 1) You may wish to define a specific layout of possible sensor locations for a specific number of sensors
% 2) You may wish to define the grid as a function of the value of k_min
N_rings = 15;                                         % number of rings
R_rings = linspace(0,0.6/min(k_min_vec), N_rings);    % the radius of each ring [meters] 

Ns_ring=zeros(1, N_rings);                   % number of possible sensor locations on the corresponding ring
for nn=1:N_rings
    if R_rings(nn) == 0                      % the central sensor (ring with radius 0)
        Ns_ring(nn)=1;
    elseif R_rings(nn)/max(R_rings) < 0.25   % some innermost rings
        Ns_ring(nn)=12;
    elseif R_rings(nn)/max(R_rings) < 0.5    % some intermediate rings
        Ns_ring(nn)=18; % 8*9
    elseif R_rings(nn)/max(R_rings) < 0.75   % some more intermediate rings
        Ns_ring(nn)=24; % 8*5*3
    else
        Ns_ring(nn) = 30;                    % outermost rings
    end
end

% additional linear constraints
CENTRAL_SENSOR = -1;  % constrains the central sensor -1 unconstrained, 0 force absence, 1 force presence of sensor
MOI_RING = true;      % constrains the moment of intertia (MOI) on each ring
MOI_GLOBAL = true;    % constrains the moment of intertia (MOI) on the whole array

% solver parameters
TimeLimit=500;    % after this time [seconds] the solver is halted

%%%
%%% Code begins here
%%%
path(path,[pwd '/utils']);
assert (k_max == 1)                 % we work with normalized wavenumbers/coordinates, k_max must be 1
assert (k_max > max(k_min_vec))     % k_min cannot exceed k_max
assert (isscalar(k_max))            

k_max= 2 * k_max; % ROI is twice as large

for nnn=1:numel(N_sensors_vec)
    for kkk=1:numel(k_min_vec)
        clear model;
        k_max=1;
        k_min=k_min_vec(kkk);
        
        N_sensors=N_sensors_vec(nnn);


        N_k = ceil(2*(k_max-k_min)/k_min); % related to number of samples in spatial freq
        
        % Here parameters specifying the possible sensor locations are
        % defined for a given number of sensors (N_sensors)
        % Based on heuristic, we may want to help the algorithm by
        % specifing how many rings, how many possible sensor locations on
        % each ring, whether the central sensor is constrained or not.
        switch N_sensors
            case 5
                R_rings=linspace(0,0.4/k_min,50);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                CENTRAL_SENSOR=0;
                for nn=1:N_rings
                    Ns_ring(nn)=5;
                end
            case 6
                R_rings=linspace(0,0.3/k_min,12);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                %CENTRAL_SENSOR=1;
                for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=5;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=10;
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=20;
                    else
                        Ns_ring(nn) = 30;
                    end
                end
            case 7
                R_rings=linspace(0,0.4/k_min,12);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=6;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=18;
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=24;
                    else
                        Ns_ring(nn) = 36;
                    end
                end
            case {8}
                R_rings=linspace(0,0.7/k_min,50);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                CENTRAL_SENSOR=1;
                 for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=7  ;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=7; 
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=7;
                    else
                        Ns_ring(nn) = 7;
                    end
                 end
            case {9}
                R_rings=linspace(0,0.4/k_min,20);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                 for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=9  ;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=18; % 8*9
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=21; % 8*5*3
                    else
                        Ns_ring(nn) = 9;
                    end
                 end
            case {10}
                R_rings=linspace(0,0.4/k_min,30);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                 for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=15   ;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=15; % 8*9
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=30; % 8*5*3
                    else
                        Ns_ring(nn) = 45;
                    end
                 end
            case 11
                R_rings=linspace(0,0.4/k_min,20);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=5   ;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=10; % 8*9
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=15; % 8*5*3
                    else
                        Ns_ring(nn) = 30;
                    end
                end
  
            case 14
                R_rings=linspace(0,0.4/k_min,16);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                CENTRAL_SENSOR=0;
                for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=0;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=0; 
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=7; 
                    else
                        Ns_ring(nn) = 21;
                    end
                end
            case 15
                R_rings=linspace(0,0.5/k_min,12);
                N_rings=numel(R_rings);
                Ns_ring=zeros(1, N_rings);
                for nn=1:N_rings
                    if R_rings(nn) == 0
                        Ns_ring(nn)=1;
                    elseif R_rings(nn)/max(R_rings) < 0.25
                        Ns_ring(nn)=12;
                    elseif R_rings(nn)/max(R_rings) < 0.5
                        Ns_ring(nn)=15; 
                    elseif R_rings(nn)/max(R_rings) < 0.75
                        Ns_ring(nn)=21; 
                    else
                        Ns_ring(nn) = 30;
                    end
                end
            otherwise
                disp('using default values for:')
                disp('   number of rings (R_rings)')
                disp('   number of sensors on each ring (Ns_rings)')

        end
        
        Nx=sum(Ns_ring); % total number of possible sensor positions
        
        
        % Create output folder(s)
        tOUTPUT_DIR=sprintf('%s/k_ratio=%.2f',OUTPUT_DIR,0.5*k_max/k_min);
        if ~exist(tOUTPUT_DIR,'dir')
            mkdir(tOUTPUT_DIR)
        end
        
        FILENAME=sprintf('MIP_Ns=%.2d',N_sensors);
        FILENAME_MAT=sprintf('%s.mat',FILENAME);
        
        
        if ~exist('model', 'var')
            
            %% linear constraints Ax=b
            % TODO here one could add linear constraints to account for the
            % presence of obstacles
            A_sensors = ones(1,Nx);
            b_sensors = N_sensors;
            
            % frequency sampling points
            K_vec=linspace(k_min, k_max, N_k);
            clear CostMask_Azimuth;
            clear CostMask_K;
            mm=0;
            for kk=1:numel(K_vec)
                Azimuth_vec=linspace(-pi/2,pi/2,ceil(pi*K_vec(kk)/(K_vec(2)-K_vec(1)))); Azimuth_vec(end)=[];
                for aa=1:numel(Azimuth_vec)
                    mm=mm+1;
                    CostMask_Azimuth(mm)=Azimuth_vec(aa);
                    CostMask_K(mm)=K_vec(kk);
                end
            end
            M=mm;
            
            % possible spatial positions
            clear p;
            xx=1;
            for nn1=1:N_rings
                for nn2=1:Ns_ring(nn1)
                    p(xx,:)=R_rings(nn1)*[cos(nn2*2*pi/Ns_ring(nn1)) sin(nn2*2*pi/Ns_ring(nn1))];
                    xx = xx +1;
                end
            end
            
            % linear constraints
            
            if MOI_RING % sm=0,st=0,Qst=0 per ring
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*cos(nn2*2*pi/Ns_ring(nn1))*sin(nn2*2*pi/Ns_ring(nn1)); % st mean = zero
                        xx = xx +1;
                    end
                end
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*cos(nn2*2*pi/Ns_ring(nn1)); % s mean zero
                        xx = xx +1;
                    end
                end
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*sin(nn2*2*pi/Ns_ring(nn1)); % t mean zero
                        xx = xx +1;
                    end
                end
                % sm=0,st=0,Qst=0 per ring,% +45 deg rotation
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*cos(nn2*2*pi/Ns_ring(nn1)+pi/4)*sin(nn2*2*pi/Ns_ring(nn1)+pi/4); % st mean = zero
                        xx = xx +1;
                    end
                end
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*cos(nn2*2*pi/Ns_ring(nn1)+pi/4); % s mean zero
                        xx = xx +1;
                    end
                end
                xx=1;
                for nn1=1:N_rings
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn2=1:Ns_ring(nn1)
                        A_sensors(end, xx )=R_rings(nn1)*sin(nn2*2*pi/Ns_ring(nn1)+pi/4); % t mean zero
                        xx = xx +1;
                    end
                end
            end
            if MOI_GLOBAL % sm=0,st=0,Qst=0 whole array
                
                % sm=0,st=0,Qst=0 per ring,% +45 deg rotation
                for theta_r=[0 pi/4];
                    
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn1=1:Nx
                        A_sensors(end, nn1 )=(cos(theta_r)*p(nn1,1)+sin(theta_r)*p(nn1,2))*(-sin(theta_r)*p(nn1,1)+cos(theta_r)*p(nn1,2)); % st mean = zero
                    end
                    
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn1=1:Nx
                        A_sensors(end, nn1 )=cos(theta_r)*p(nn1,1)+sin(theta_r)*p(nn1,2); % s mean zero
                    end
                    
                    A_sensors=[A_sensors; zeros(1,Nx)];
                    b_sensors=[b_sensors; 0];
                    for nn1=1:Nx
                        A_sensors(end, nn1 )=-sin(theta_r)*p(nn1,1)+cos(theta_r)*p(nn1,2); % t mean zero
                    end
                end
            end
            
            if CENTRAL_SENSOR == 1 || CENTRAL_SENSOR == 0
                A_sensors=[A_sensors; zeros(1,Nx)];
                b_sensors=[b_sensors; CENTRAL_SENSOR];
                A_sensors(end, 1)=1;
            end
            
            N_lc=size(A_sensors,1);
            
            disp('Building Fourier operator')
            F_operator_c=complex(zeros(M,Nx));
            for mm=1:M
                for nn=1:Nx
                    F_operator_c(mm,nn) = exp(-2*pi*1i*[CostMask_K(mm)*cos(CostMask_Azimuth(mm)), CostMask_K(mm)*sin(CostMask_Azimuth(mm))]*transpose(p(nn,:)));
                end
            end
            
            
            disp('Creating Gurobi model')
            model.obj = [1; zeros(Nx,1)];
            model.modelsense = 'min';
            
            model.A = sparse([
                -ones(M,1), real(F_operator_c)/N_sensors;
                -ones(M,1),  imag(F_operator_c)/N_sensors;
                -ones(M,1),  -real(F_operator_c)/N_sensors;
                -ones(M,1), -imag(F_operator_c)/N_sensors;
                zeros(N_lc,1), A_sensors;
                ]);
            model.rhs = [zeros(4*M,1); b_sensors;];
            model.sense = [repmat('<', [4*M 1]); repmat('=', [N_lc 1]);];
            
            clear F_operator_c;
            
            model.vtype = ['C' repmat('B', [1 Nx])];
            
            
            params.outputflag = 1;
            params.resultfile = '';
            params.MIPFocus = 1;
            params.TimeLimit = TimeLimit; % [s]
            
        end
        
        if exist('x_MIP', 'var')
            disp('Starting gurobi from existing solution')
            model.start = zeros(1+Nx,1);
            model.start(1) = NaN;
            model.start(2:end) = x_MIP;
        end
        
        try
            
            disp('Calling Gurobi')
            result = gurobi(model, params);
            
            fprintf('Obj: %e\n', result.objval);
            
        catch gurobiError
            fprintf('Error reported\n');
            gurobiError
        end
        
        if exist('result','var') && isfield(result,'x') %% solution was found
            
            x_sol = result.x(2:end);
            pos_sol=p(logical(x_sol),:);
            
            
            
            if exist('FILENAME_MAT','var')
                disp('Loading existing file')
                if exist(sprintf('%s/%s',tOUTPUT_DIR,FILENAME_MAT),'file')
                    TMP=load(sprintf('%s/%s',tOUTPUT_DIR,FILENAME_MAT));
                    [TMPbp TMPmax_bp]=plotBeampattern_C(TMP.pos_sol, TMP.k_min, TMP.k_max);
                    close
                    [bp max_bp]=plotBeampattern_C(pos_sol, k_min, k_max);
                    close
                    clear TMP;
                    if TMPmax_bp > max_bp % check if existing file has better solution or not
                        %             clear max_bp
                        %             clear TMPmax_bp
                        disp(sprintf('Saving results on disk (MAT file). %s/%s',tOUTPUT_DIR,FILENAME_MAT) )
                        eval(sprintf('save %s/%s',tOUTPUT_DIR,FILENAME_MAT))
                        
                    end
                else
                    disp(sprintf('Saving results on disk (MAT file). %s/%s',tOUTPUT_DIR,FILENAME_MAT) )
                    eval(sprintf('save %s/%s',tOUTPUT_DIR,FILENAME_MAT))
                    
                end
            end
            
            
        end
        
    end
end

%%%
%%% Some plotting
%%%

figure
subplot(2,2,1)
plot(p(:,1),p(:,2),'o')
hold on
plot(pos_sol(:,1),pos_sol(:,2),'r.')
xlabel('[m]'); ylabel('[m]');
title('Possible sensor positions')
axis square

subplot(2,2,2)
plotArray_C(pos_sol)

subplot(2,2,3)
[bp max_bp]=plotBeampattern_C(pos_sol, k_min, k_max);
plotWavenumberLimits(k_min, k_max)
title(max_bp)
xlabel('[1/m]'); ylabel('[1/m]');
axis square

subplot(2,2,4)
plot(CostMask_K.*cos(CostMask_Azimuth), CostMask_K.*sin(CostMask_Azimuth),'.')
hold on
plot(CostMask_K.*cos(CostMask_Azimuth+pi), CostMask_K.*sin(CostMask_Azimuth+pi),'o')
plotWavenumberLimits(k_min, k_max)
hold off
xlim(1.5*k_max*[-1 1 ])
ylim(1.5*k_max*[-1 1 ])
xlabel('[1/m]'); ylabel('[1/m]');
axis square
title(M)

