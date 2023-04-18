% BSL comsol input code
% The code makes a comsol input file (txt) from a matlab struct (.mat).

% assumes .mat has a table I, T, varaibles(SOC)


%% Interface

% .mat (struct) fullpath
origin_fullpath ='C:\Users\jsong\Downloads\fullfile_Final.mat';
Ra_fullpath = 'C:\Users\jsong\Downloads\cellmodel_I_T_SOC_Ra.txt';
Rc_fullpath = 'C:\Users\jsong\Downloads\cellmodel_I_T_SOC_Rc.txt';

% parameters
A = 1; % scaling factor
Roff = 0; % correction factor







%% Engine

% load cell model 
load(origin_fullpath) % assume variable name 'data_accum'


% Crate vector
Crate_column = cell2mat(data_accum.Crate); % the table column is in 'cell'. convert cell to matrix (double).
Crate_unique = unique(Crate_column);

% Temperature vector
T_column = cell2mat(data_accum.Temperature);
T_unique = unique(T_column);


% make comsol input

Ra_comsol = [];
Rc_comsol = [];

for i = 1:length(Crate_unique)
    for j = 1:length(T_unique)
        
        % i, j - th data to write in .txt
        SOC_vec_ij = data_accum.SOC{(i-1)*length(T_unique)+j};
        Ra_vec_ij = data_accum.Resistance_Anode{(i-1)*length(T_unique)+j};
        Rc_vec_ij = data_accum.Resistance_Cathode{(i-1)*length(T_unique)+j};

        % length of the i, j - vector
        L_ij = size(data_accum.SOC{(i-1)*length(T_unique)+j});
                    % SOC, R_a, R_c, and other simulated data have the same length 
 
        Ra_comsol = [Ra_comsol;...
                    Crate_unique(i)*ones(L_ij), T_unique(j)*ones(L_ij), SOC_vec_ij, Ra_vec_ij];

        Rc_comsol = [Rc_comsol;...
                    Crate_unique(i)*ones(L_ij), T_unique(j)*ones(L_ij), SOC_vec_ij, Rc_vec_ij];
                    % R_a and R_c have the same length 


    end
end


% correction, if needed
Ra_comsol(:,4) = (Ra_comsol(:,4)-Roff)*A;
Rc_comsol(:,4) = (Rc_comsol(:,4)-Roff)*A;



% save
writematrix(Ra_comsol,Ra_fullpath);
writematrix(Ra_comsol,Rc_fullpath);
