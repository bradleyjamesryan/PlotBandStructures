clear; clc; close all;
tic
cd '/home/anon/Documents/GraduateSchool/Data/Group IV/DFT/2020_01_06_SpinPolarizedDefects/HVacancy'
filez = dir('*.xml');
ylimz = [-3 4];
set(0,'defaultAxesFontName', 'Arial');
for randvar = 1:length({filez(:).name})
    axLineWidth = 1.5; plotLineWidth = 1.7; gridTransparency = .15; textSize = 7; LabelSize = 14; titleFontSize = 24; AxisNumberSize = 18; LegFontSize = 14; colorz = [1 0 0; 0 0 1; 0 1 0]; alpha = 1;xlabelz = {'K' '\Gamma' 'M' 'K'};
    if (exist([filez(randvar).name(1:end-4) '_SiOrbital' '.png'],'file')==2) && (exist([filez(randvar).name(1:end-4) '_AllAtoms' '.png'],'file')==2) && (exist([filez(randvar).name(1:end-4) '_SiH' '.png'],'file')==2)
        continue
    else
        %% Misc values
        if  (exist(['data_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
            load(['data_',filez(randvar).name(1:end-4),'.mat'],'data');
        else
            data = xml2struct(filez(randvar).name);
            save(['data_',filez(randvar).name(1:end-4),'.mat'],'data');
        end
        if  (exist(['natom_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['ntype_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['atom_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['atomnames_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['fermienergy_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
            
            load(['natom_',filez(randvar).name(1:end-4),'.mat'],'natom');
            load(['ntype_',filez(randvar).name(1:end-4),'.mat'],'ntype');
            load(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
            load(['atomnames_',filez(randvar).name(1:end-4),'.mat'],'atomnames');
            load(['fermienergy_',filez(randvar).name(1:end-4),'.mat'], 'fermienergy');
        else
            natom = str2num(data.modeling.atominfo.atoms.Text); % Number of atoms in unit cell
            for ii = 1:natom
                atomz(ii) = str2num(data.modeling.atominfo.array{1}.set.rc{ii}.c{2}.Text); % Array of numeric identifiers of atoms
                atomnames_dummy{ii} = [data.modeling.atominfo.array{1}.set.rc{ii}.c{1}.Text];
            end
            
            ntype = length(unique(atomz)); % Number of different typoes of atoms
            
            for ii = 1:ntype
                atom{ii} =  [sum(atomz==ii)];
            end
            atomnames = unique(atomnames_dummy,'stable');
            save(['natom_',filez(randvar).name(1:end-4),'.mat'],'natom');
            save(['ntype_',filez(randvar).name(1:end-4),'.mat'],'ntype');
            save(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
            save(['atomnames_',filez(randvar).name(1:end-4),'.mat'],'atomnames');
        end
        
        %% Loading Plotting Data
        if (exist(['nband_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['fermienergy_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['band_xaxis_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['numkpoints_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['projected_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['partial_up_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)  &&...
                (exist(['partial_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['total_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['total_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['bands_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['recip_basis_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['symm_points_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['band_DLocx_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['projected_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) &&...
                (exist(['bands_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
            
            load(['nband_',filez(randvar).name(1:end-4),'.mat'], 'nband');
            load(['fermienergy_',filez(randvar).name(1:end-4),'.mat'], 'fermienergy');
            load(['band_xaxis_',filez(randvar).name(1:end-4),'.mat'], 'band_xaxis');
            load(['numkpoints_',filez(randvar).name(1:end-4),'.mat'], 'numkpoints');
            load(['projected_',filez(randvar).name(1:end-4),'.mat'], 'projected');
            load(['partial_up_',filez(randvar).name(1:end-4),'.mat'], 'partial_up');
            load(['partial_down_',filez(randvar).name(1:end-4),'.mat'], 'partial_down');
            load(['total_',filez(randvar).name(1:end-4),'.mat'], 'total');
            load(['total_down_',filez(randvar).name(1:end-4),'.mat'], 'total_down');
            load(['bands_',filez(randvar).name(1:end-4),'.mat'], 'bands');
            load(['recip_basis_',filez(randvar).name(1:end-4),'.mat'],'recip_basis');
            load(['symm_points_',filez(randvar).name(1:end-4),'.mat'],'symm_points');
            load(['band_x_',filez(randvar).name(1:end-4),'.mat'], 'band_x');
            load(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat'], 'cbm_vbm');
            load(['projected_down_',filez(randvar).name(1:end-4),'.mat'],'projected_down');
            load(['bands_down_',filez(randvar).name(1:end-4),'.mat'], 'bands_down');
        else
            %% x-axis for band structure
            for ii = 1:length(data.modeling.kpoints.generation.v)
                band_x{ii} = str2num(data.modeling.kpoints.generation.v{ii}.Text);
            end
            numkpoints = str2num(data.modeling.kpoints.generation.i.Text);
            
            % Calculating the recriprocal-lattice-special-points in the basis of the recriprocal lattice
            for ii = 1:3
                recip_basis{ii} =  str2num(data.modeling.structure{1}.crystal.varray{2}.v{ii}.Text);
            end
            
            % A bit of magic here. Even I am not sure what exactly is going
            % on, but (I believe) it works.
            for ii = 1:length(data.modeling.kpoints.generation.v)
                symm_points{ii} = str2num(data.modeling.kpoints.generation.v{ii}.Text);
                for jj = 1:3
                    recip_symm_points(ii,jj) = symm_points{ii}(1)*recip_basis{1}(jj) + symm_points{ii}(2)*recip_basis{2}(jj) + symm_points{ii}(3)*recip_basis{3}(jj);
                end
            end
            
            band_xaxis_old = zeros(numkpoints,1);
            for ii = 1:length(data.modeling.kpoints.generation.v)-1
                dist_betwn_recip_symm_points(ii) = sqrt((recip_symm_points(ii+1,1) - recip_symm_points(ii,1))^2 + (recip_symm_points(ii+1,2) - recip_symm_points(ii,2))^2 + (recip_symm_points(ii+1,3) - recip_symm_points(ii,3))^2)./(numkpoints-1);
                discretized_space_betwn_recip_symm_points(:,ii) = linspace(0,dist_betwn_recip_symm_points(ii).*numkpoints,numkpoints+1);
                band_xaxis_old((ii-1).*numkpoints+1:ii.*numkpoints,1) = discretized_space_betwn_recip_symm_points(1:end-1,ii) + band_xaxis_old(end);
            end
            
            save(['band_x_',filez(randvar).name(1:end-4),'.mat'], 'band_x');
            save(['numkpoints_',filez(randvar).name(1:end-4),'.mat'], 'numkpoints');
            save(['recip_basis_',filez(randvar).name(1:end-4),'.mat'],'recip_basis');
            save(['symm_points_',filez(randvar).name(1:end-4),'.mat'],'symm_points');
            
            %% Projected, DOS, and Bands
            nband = length(data.modeling.calculation.projected.array.set.set{1}.set{1}.set);
            for kpoint = 1:length(data.modeling.calculation.projected.array.set.set{1}.set)
                for band = 1:nband
                    for atm = 1:length(data.modeling.calculation.projected.array.set.set{1}.set{1}.set{1}.r)
                        projected{kpoint}{band}(atm,:) = str2num(data.modeling.calculation.projected.array.set.set{1}.set{kpoint}.set{band}.r{atm}.Text); % Projected: projected = [s py pz px dxy dyz dz2 dxz x2-y2]
                        projected_down{kpoint}{band}(atm,:) = str2num(data.modeling.calculation.projected.array.set.set{2}.set{kpoint}.set{band}.r{atm}.Text);
                    end
                end
                for A = 1:length(data.modeling.calculation.eigenvalues.array.set.set{1}.set{kpoint}.r)
                    bands_old{kpoint}{A} = str2num(data.modeling.calculation.eigenvalues.array.set.set{1}.set{kpoint}.r{A}.Text);% Bands: bands_old = [energy, occypancy]
                    bands{kpoint}(A) = bands_old{kpoint}{A}(1);
                    
                    bands_old_down{kpoint}{A} = str2num(data.modeling.calculation.eigenvalues.array.set.set{2}.set{kpoint}.r{A}.Text);% Bands: bands_old = [energy, occypancy]
                    bands_down{kpoint}(A) = bands_old_down{kpoint}{A}(1);
                end
                band_xaxis{kpoint} = linspace(band_xaxis_old(kpoint),band_xaxis_old(kpoint),length(bands{kpoint}));
            end
            
            save(['band_xaxis_',filez(randvar).name(1:end-4),'.mat'], 'band_xaxis');
            save(['nband_',filez(randvar).name(1:end-4),'.mat'], 'nband');
            save(['projected_',filez(randvar).name(1:end-4),'.mat'], 'projected');
            save(['projected_down_',filez(randvar).name(1:end-4),'.mat'],'projected_down');
            
            %% Calculating VBM, CBM, and Fermi Energy
            fermienergy = str2num(data.modeling.calculation.dos.i.Text);            vbm = -1E10;            cbm = 1E10;
            for ii = 1:length(bands_old)
                for jj = 1:length(bands_old{ii})
                    if bands_old{ii}{jj}(1) > fermienergy && bands_old{ii}{jj}(1) < cbm
                        cbm = bands_old{ii}{jj}(1);
                    elseif bands_old{ii}{jj}(1) < fermienergy && bands_old{ii}{jj}(1) > vbm
                        vbm = bands_old{ii}{jj}(1);
                    end
                end
            end
            cbm_vbm = [cbm vbm];
            fermienergy = fermienergy - cbm_vbm(2);
            
            save(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat'], 'cbm_vbm');
            save(['fermienergy_',filez(randvar).name(1:end-4),'.mat'],'fermienergy');
            
            %% Extracting total and partial DOS
            for atm = 1:length(data.modeling.calculation.dos.partial.array.set.set)
                for A = 1:length(data.modeling.calculation.dos.partial.array.set.set{atm}.set{1}.r)
                    partial_up{atm}(A,:) = str2num(data.modeling.calculation.dos.partial.array.set.set{atm}.set{1}.r{A}.Text); % pDOS: partial = [f2]
                    partial_down{atm}(A,:) = str2num(data.modeling.calculation.dos.partial.array.set.set{atm}.set{2}.r{A}.Text); % pDOS: partial = [f2]
                    
                    total(A,:) = str2num(data.modeling.calculation.dos.total.array.set.set{1}.r{A}.Text);  % Total DOS: total = [energy total integrated]
                    total_down(A,:) = str2num(data.modeling.calculation.dos.total.array.set.set{2}.r{A}.Text);  % Total DOS: total = [energy total integrated]
                end
            end
            
            save(['partial_up_',filez(randvar).name(1:end-4),'.mat'], 'partial_up');
            save(['partial_down_',filez(randvar).name(1:end-4),'.mat'], 'partial_down');
            save(['total_',filez(randvar).name(1:end-4),'.mat'], 'total');
            save(['total_down_',filez(randvar).name(1:end-4),'.mat'], 'total_down');
            save(['bands_',filez(randvar).name(1:end-4),'.mat'], 'bands');
            save(['bands_down_',filez(randvar).name(1:end-4),'.mat'],'bands_down');
            
        end
        
        %         load(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
        xlower = min([band_xaxis{:}]); xupper = max([band_xaxis{:}]);
        
        for ii = 1:length(atomnames)
            if atomnames{ii} == 'Si'; SiLoc = ii; end
            if atomnames{ii} == 'H '; HLoc = ii;  end
        end
        sumz = HLoc + SiLoc;
        if sumz == 3;            DLoc = 3;        end
        if sumz == 4;            DLoc = 2;        end
        if sumz == 5;            DLoc = 1;        end
        
        %% 1/3 Projected DOS for Si(px), Si(py), Si(pz)
        if  (exist(['band_x_plot_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['band_plot_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color1_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color2_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color3_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color1_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color2_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color3_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['actymin_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['actymax_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['DOS_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['DOS_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['band_plot_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['bands_down_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
            
            load(['band_x_plot_',filez(randvar).name(1:end-4),'.mat'], 'band_x_plot');
            load(['band_plot_',filez(randvar).name(1:end-4),'.mat'], 'band_plot');
            load(['color1_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color1');
            load(['color2_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color2');
            load(['color3_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color3');
            load(['color1_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
            load(['color2_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
            load(['color3_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
            load(['actymin_',filez(randvar).name(1:end-4),'.mat'], 'actymin');
            load(['actymax_',filez(randvar).name(1:end-4),'.mat'], 'actymax');
            load(['DOS_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
            load(['DOS_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
            load(['band_plot_down_',filez(randvar).name(1:end-4),'.mat'], 'band_plot_down');
            load(['bands_down_',filez(randvar).name(1:end-4),'.mat'], 'bands_down');
        else
            
            %% Calculating percentages of pDOS (to determine colors for plotting bands)
            % Projected: orbital = [s py pz px dxy dyz dz2 dxz x2-y2]
            % Projected: projected{kpoint}{band}(atom,orbital)
            spdf_vals = {2,3,4}; %orbital = [s py pz px dxy dyz dz2 dxz x2-y2]
            range{1} = 1:atom{1};
            for kp = 1:length(projected)
                for bnd = 1:length(projected{kp})
                    for spdf = 1:3
                        A = sum(sum(projected{kp}{bnd}(:,:)))+eps;
                        atomic_projected{kp}{bnd}{1}(spdf) = sum(sum(projected{kp}{bnd}(range{1},spdf_vals{spdf})));
                        percent_atomic_projected{kp}{bnd}{1}(spdf) = atomic_projected{kp}{bnd}{1}(spdf)./A;
                        
                        A2 = sum(sum(projected_down{kp}{bnd}(:,:)))+eps;
                        atomic_projected_down{kp}{bnd}{1}(spdf) = sum(sum(projected_down{kp}{bnd}(range{1},spdf_vals{spdf})));
                        percent_atomic_projected_down{kp}{bnd}{1}(spdf) = atomic_projected_down{kp}{bnd}{1}(spdf)./A2;
                    end
                end
            end
            
            for kp = 1:length(projected)
                for bnd = 1:length(projected{kp})
                    for atm = 2:length(atom)
                        range{atm} = (sum([atom{1:atm-1}])+1):(sum([atom{1:atm-1}])+atom{atm});
                        for spdf = 1:3
                            A = sum(sum(projected{kp}{bnd}(:,:)))+eps;
                            atomic_projected{kp}{bnd}{atm}(spdf) = sum(sum(projected{kp}{bnd}(range{atm},spdf_vals{spdf})));
                            percent_atomic_projected{kp}{bnd}{atm}(spdf) = atomic_projected{kp}{bnd}{atm}(spdf)./A;
                            
                            A2 = sum(sum(projected_down{kp}{bnd}(:,:)))+eps;
                            atomic_projected_down{kp}{bnd}{atm}(spdf) = sum(sum(projected_down{kp}{bnd}(range{atm},spdf_vals{spdf})));
                            percent_atomic_projected_down{kp}{bnd}{atm}(spdf) = atomic_projected_down{kp}{bnd}{atm}(spdf)./A2;
                        end
                    end
                end
            end
            
            %% Calculating Colors from pDOS percentages
            band_x_plot = cell(length(band_xaxis{1}),1);
            band_plot = cell(length(band_xaxis{1}),1);
            band_plot_down = cell(length(band_xaxis{1}),1);
            isymin(1) = min([bands{1}]); isymax(1) = max([bands{1}]);
            color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
            color1_down = cell(length(band_xaxis{1}),1); color2_down = cell(length(band_xaxis{1}),1); color3_down = cell(length(band_xaxis{1}),1);
            for bnd = 1:length(band_xaxis{1})
                band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{1}(bnd)];
                band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
                band_plot_down{bnd} = [band_plot_down{bnd} bands_down{1}(bnd)];
                color1{bnd} = [color1{bnd} percent_atomic_projected{1}{bnd}{SiLoc}(3)];
                color2{bnd} = [color2{bnd} percent_atomic_projected{1}{bnd}{SiLoc}(1)];
                color3{bnd} = [color3{bnd} percent_atomic_projected{1}{bnd}{SiLoc}(2)];
                
                color1_down{bnd} = [color1_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(3)];
                color2_down{bnd} = [color2_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(1)];
                color3_down{bnd} = [color3_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(2)];
                for kp = 2:length(band_xaxis)
                    band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{kp}(bnd)];
                    band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
                    band_plot_down{bnd} = [band_plot_down{bnd} bands_down{kp}(bnd)];
                    color1{bnd} = [color1{bnd} percent_atomic_projected{kp}{bnd}{SiLoc}(3)];
                    color2{bnd} = [color2{bnd} percent_atomic_projected{kp}{bnd}{SiLoc}(1)];
                    color3{bnd} = [color3{bnd} percent_atomic_projected{kp}{bnd}{SiLoc}(2)];
                    
                    color1_down{bnd} = [color1_down{bnd} percent_atomic_projected_down{kp}{bnd}{SiLoc}(3)];
                    color2_down{bnd} = [color2_down{bnd} percent_atomic_projected_down{kp}{bnd}{SiLoc}(1)];
                    color3_down{bnd} = [color3_down{bnd} percent_atomic_projected_down{kp}{bnd}{SiLoc}(2)];
                    isymin(kp) = min([bands{kp}]); isymax(kp) = max([bands{kp}]);
                end
                %                 Commented section below enables adjustng kpoint path
                %                 a = linspace(0,49*(band_x_plot{bnd}(102)-band_x_plot{bnd}(101)),50);
                %                 b = linspace(a(end),49*(band_x_plot{bnd}(2)-band_x_plot{bnd}(1))+a(end),50);
                %                 c = linspace(b(end),49*(band_x_plot{bnd}(52)-band_x_plot{bnd}(51))+b(end),50);
            end
            actymin = min(isymin);    actymax = max(isymax);
            
            %% Calculating DOS: DOS{atom_number}{spdf_val}
            for ii = 1:ntype
                for jj = 1:length(spdf_vals)
                    A = zeros(length(partial_up{range{ii}(1)}(:,spdf_vals{jj}+1)),1);
                    A2 = zeros(length(partial_down{range{ii}(1)}(:,spdf_vals{jj}+1)),1);
                    for kk = 1:length(range{ii})
                        A = A + partial_up{range{ii}(kk)}(:,spdf_vals{jj}+1)+eps;
                        A2 = A2 + partial_down{range{ii}(kk)}(:,spdf_vals{jj}+1)+eps;
                    end
                    DOS{ii}{jj} = A;
                    DOS_down{ii}{jj} = A2;
                end
            end
            
            save(['band_x_plot_',filez(randvar).name(1:end-4),'.mat'], 'band_x_plot');
            save(['band_plot_',filez(randvar).name(1:end-4),'.mat'], 'band_plot');
            save(['color1_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color1');
            save(['color2_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color2');
            save(['color3_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color3');
            save(['color1_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
            save(['color2_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
            save(['color3_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
            save(['actymin_',filez(randvar).name(1:end-4),'.mat'], 'actymin');
            save(['actymax_',filez(randvar).name(1:end-4),'.mat'], 'actymax');
            save(['DOS_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
            save(['DOS_down_SiOrbital_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
            save(['band_plot_down',filez(randvar).name(1:end-4),'.mat'], 'band_plot_down');
            save(['bandplot_',filez(randvar).name(1:end-4),'.mat'],'band_plot');
            
        end
        xtickz(1) = 0;        count = 1;
        for ii = [1 2]
            xtickz(count+1) = band_x_plot{1}(numkpoints*ii);
            count = count + 1;
        end
        
        xtickz(end+1) = max([band_xaxis{:}]);
        amnt_down = -0.05;         amnt_up = -0.04;
        f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); subplot('Position',[0.12 0.15 0.6 0.82]); hold on;
        xlim([xlower xupper]);    ylim(ylimz)
        for ii = [1 2]
            h = plot([band_x_plot{1}(numkpoints*ii) band_x_plot{1}(numkpoints*ii)],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
        end
        for bnd = 1:nband
            for kp = 2:150
                plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot{bnd}(kp-1)-cbm_vbm(2), band_plot{bnd}(kp)-cbm_vbm(2)],'color',[color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)],'linewidth',plotLineWidth)
                plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot_down{bnd}(kp-1)-cbm_vbm(2), band_plot_down{bnd}(kp)-cbm_vbm(2)],'color',[color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)],'linewidth',plotLineWidth)
            end
            for kp = 3:4:100
                t = text(([band_x_plot{bnd}(kp-2)]), band_plot{bnd}(kp-2)-cbm_vbm(2)-amnt_up,'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
                t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2)-amnt_down,'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
            end
            for kp = 103:8:150
                t = text(([band_x_plot{bnd}(kp-4)]), band_plot{bnd}(kp-4)-cbm_vbm(2)-amnt_up,'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
                t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2)-amnt_down,'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
            end
        end
        
        %         ▲▼△▽˅˄
        plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
        ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');        xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
        set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
        set(gca,'xtick',xtickz,'xticklabel',xlabelz);
        titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
        
        % Plotting DOS
        subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
        p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;% Total
        plot(total(:,1)-cbm_vbm(2),total(:,2), 'color', [0 0 0],'linewidth',1);% Total
        
        d = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{3}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(px)
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{3}, 'color', [1 1 1],'linewidth',1);% Si(px)

        d = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(py)
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1}, 'color', [1 1 1],'linewidth',1);% Si(py)

        d = area(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2)); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(pz)
        plot(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2), 'color', [1 1 1],'linewidth',1);% Si(pz)

        h3 = area(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = alpha;% Si(pz)
        plot(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2), 'color', [0 0 1],'linewidth',1);% Si(pz)

        h2 = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = alpha;% Si(py)
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1}, 'color', [0 1 0],'linewidth',1);% Si(py)

        h1 = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{3}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = alpha;% Si(px)                
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{3}, 'color', [1 0 0],'linewidth',1);% Si(px)
        
        
       
        

        p = area(total(:,1)-cbm_vbm(2),-total_down(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;% Total
        plot(total(:,1)-cbm_vbm(2),-total_down(:,2), 'color', [0 0 0],'linewidth',1);% Total

        d = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{3}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(px)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{3}, 'color', [1 1 1],'linewidth',1);% Si(px)

        d = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(py)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1}, 'color', [1 1 1],'linewidth',1);% Si(py)

        
        d = area(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2)); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(pz)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2), 'color', [1 1 1],'linewidth',1);% Si(pz)

        h3 = area(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = alpha;% Si(pz)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2), 'color', [0 0 1],'linewidth',1);% Si(pz)

        h2 = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = alpha;% Si(py)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1}, 'color', [0 1 0],'linewidth',1);% Si(py)

        h1 = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{3}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = alpha;% Si(px)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{3}, 'color', [1 0 0],'linewidth',1);% Si(px)

        
        
        leg = legend([p h1 h2 h3],{'Total','Si(px)','Si(py)','Si(pz)'},'location','south','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
        view([90 -90]);
        %         ylim([0 1.05.*max(total(:,2))]);        
        xlim(ylimz)
        ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
        set(gca,'xticklabel','');        set(gca,'linewidth',axLineWidth,'XMinorTick','on','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
        plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
        saveas(f,[filez(randvar).name(1:end-4) '_SiOrbital' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_SiOrbital' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_SiOrbital' '.fig']);
        close(f)
        
        %% 2/3 Projected DOS for H(s), Si(s), and Si(p)
        if  (exist(['band_x_plot_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['band_plot_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color1_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color2_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color3_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color1_down_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color2_down_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['color3_down_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['DOS_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
                (exist(['DOS_down_SiH_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
            
            load(['color1_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color1');
            load(['color2_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color2');
            load(['color3_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color3');
            load(['color1_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
            load(['color2_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
            load(['color3_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
            load(['DOS_SiH_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
            load(['DOS_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
        else
            %% Calculating percentages of pDOS (to determine colors for plotting bands)
            % Projected: orbital = [s py pz px dxy dyz dz2 dxz x2-y2]
            % Projected: projected{kpoint}{band}(atom,orbital)
            spdf_vals = {1,2:4,5:9};        range{1} = 1:atom{1};
            for kp = 1:length(projected)
                for bnd = 1:length(projected{kp})
                    for spdf = 1:3
                        A = sum(sum(projected{kp}{bnd}(:,:)))+eps;
                        atomic_projected{kp}{bnd}{1}(spdf) = sum(sum(projected{kp}{bnd}(range{1},spdf_vals{spdf})));
                        percent_atomic_projected{kp}{bnd}{1}(spdf) = atomic_projected{kp}{bnd}{1}(spdf)./A;
                        
                        A2 = sum(sum(projected_down{kp}{bnd}(:,:)))+eps;
                        atomic_projected_down{kp}{bnd}{1}(spdf) = sum(sum(projected_down{kp}{bnd}(range{1},spdf_vals{spdf})));
                        percent_atomic_projected_down{kp}{bnd}{1}(spdf) = atomic_projected_down{kp}{bnd}{1}(spdf)./A2;
                    end
                end
            end
            for kp = 1:length(projected)
                for bnd = 1:length(projected{kp})
                    for atm = 2:length(atom)
                        range{atm} = (sum([atom{1:atm-1}])+1):(sum([atom{1:atm-1}])+atom{atm});
                        for spdf = 1:3
                            A = sum(sum(projected{kp}{bnd}(:,:)))+eps;
                            atomic_projected{kp}{bnd}{atm}(spdf) = sum(sum(projected{kp}{bnd}(range{atm},spdf_vals{spdf})));
                            percent_atomic_projected{kp}{bnd}{atm}(spdf) = atomic_projected{kp}{bnd}{atm}(spdf)./A;
                            
                            A2 = sum(sum(projected_down{kp}{bnd}(:,:)))+eps;
                            atomic_projected_down{kp}{bnd}{atm}(spdf) = sum(sum(projected_down{kp}{bnd}(range{atm},spdf_vals{spdf})));
                            percent_atomic_projected_down{kp}{bnd}{atm}(spdf) = atomic_projected_down{kp}{bnd}{atm}(spdf)./A2;
                        end
                    end
                end
            end
            
            %% Calculating Colors from pDOS percentages
            color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
            color1_down = cell(length(band_xaxis{1}),1); color2_down = cell(length(band_xaxis{1}),1); color3_down = cell(length(band_xaxis{1}),1);
            
            for bnd = 1:length(band_xaxis{1})
                color1{bnd} = [color1{bnd} percent_atomic_projected{1}{bnd}{HLoc}(1)];
                color2{bnd} = [color2{bnd} percent_atomic_projected{1}{bnd}{SiLoc}(1)];
                color3{bnd} = [color3{bnd} percent_atomic_projected{1}{bnd}{SiLoc}(2)];
                
                color1_down{bnd} = [color1_down{bnd} percent_atomic_projected_down{1}{bnd}{HLoc}(1)];
                color2_down{bnd} = [color2_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(1)];
                color3_down{bnd} = [color3_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(2)];
                for kp = 2:length(band_xaxis)
                    color1{bnd} = [color1{bnd} percent_atomic_projected{kp}{bnd}{HLoc}(1)];
                    color2{bnd} = [color2{bnd} percent_atomic_projected{kp}{bnd}{SiLoc}(1)];
                    color3{bnd} = [color3{bnd} percent_atomic_projected{kp}{bnd}{SiLoc}(2)];
                    
                    color1_down{bnd} = [color1_down{bnd} percent_atomic_projected_down{1}{bnd}{HLoc}(1)];
                    color2_down{bnd} = [color2_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(1)];
                    color3_down{bnd} = [color3_down{bnd} percent_atomic_projected_down{1}{bnd}{SiLoc}(2)];
                end
            end
            
            for ii = 1:ntype
                for jj = 1:length(spdf_vals)
                    A = zeros(length(partial_up{range{ii}(1)}(:,spdf_vals{jj}+1)),1)+eps;
                    A2 = zeros(length(partial_down{range{ii}(1)}(:,spdf_vals{jj}+1)),1);
                    for kk = 1:length(range{ii})
                        A = A + partial_up{range{ii}(kk)}(:,spdf_vals{jj}+1);
                        A2 = A2 + partial_down{range{ii}(kk)}(:,spdf_vals{jj}+1)+eps;
                    end
                    DOS{ii}{jj} = A;
                    DOS_down{ii}{jj} = A2;
                end
            end
            save(['color1_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color1');
            save(['color2_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color2');
            save(['color3_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color3');
            save(['color1_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
            save(['color2_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
            save(['color3_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
            save(['DOS_SiH_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
            save(['DOS_down_SiH_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
        end
        
        f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); subplot('Position',[0.12 0.15 0.6 0.82]); hold on;
        for ii = [1 2]
            h = plot([band_x_plot{1}(numkpoints*ii) band_x_plot{1}(numkpoints*ii)],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
        end
        for bnd = 1:nband
            for kp = 2:150
                plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot{bnd}(kp-1)-cbm_vbm(2), band_plot{bnd}(kp)-cbm_vbm(2)],'color',[color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)],'linewidth',plotLineWidth)
                plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot_down{bnd}(kp-1)-cbm_vbm(2), band_plot_down{bnd}(kp)-cbm_vbm(2)],'color',[color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)],'linewidth',plotLineWidth)
            end
            for kp = 3:4:100
                t = text(([band_x_plot{bnd}(kp-2)]), band_plot{bnd}(kp-2)-cbm_vbm(2),'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
                t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2),'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
            end
            for kp = 108:8:150
                t = text(([band_x_plot{bnd}(kp-4)]), band_plot{bnd}(kp-4)-cbm_vbm(2),'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';t.HorizontalAlignment = 'center';
                t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2),'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';t.HorizontalAlignment = 'center';
            end
        end
        
        plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
        xlim([xlower xupper]);    ylim(ylimz)
        ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');        xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
        set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';        set(gca,'xtick',xtickz,'xticklabel',xlabelz);
        titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
        
        % Plotting DOS
        subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
        p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
        plot(total(:,1)-cbm_vbm(2),total(:,2),'color', [0 0 0],'linewidth',1);

        h3 = area(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = alpha;% Si(p)
        plot(partial_up{1}(:,1)-cbm_vbm(2),sum(DOS{SiLoc}{2},2),'color', [0 0 1],'linewidth',1);% Si(p)

        h2 = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = alpha;% Si(s)
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{SiLoc}{1},'color', [0 1 0],'linewidth',1);% Si(s)

        h1 = area(partial_up{1}(:,1)-cbm_vbm(2),DOS{HLoc}{1}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = alpha;% H(s)
        plot(partial_up{1}(:,1)-cbm_vbm(2),DOS{HLoc}{1},'color', [1 0 0],'linewidth',1);% H(s)

        
        
        
        p = area(total(:,1)-cbm_vbm(2),-total_down(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
        plot(total(:,1)-cbm_vbm(2),-total_down(:,2),'color', [0 0 0],'linewidth',1);
        
        h3 = area(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = alpha;% Si(p)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-sum(DOS_down{SiLoc}{2},2),'color', [0 0 1],'linewidth',1);% Si(p)
        
        h2 = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = alpha;% Si(s)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{SiLoc}{1},'color', [0 1 0],'linewidth',1);% Si(s)
        
        h1 = area(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{HLoc}{1}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = alpha;% H(s)
        plot(partial_up{1}(:,1)-cbm_vbm(2),-DOS_down{HLoc}{1},'color', [1 0 0],'linewidth',1);% H(s)
        
        
        leg = legend([p h1 h2 h3],{'Total','H(s)','Si(s)','Si(p)'},'location','south','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
        view([90 -90]);
        %         ylim([0 1.05.*max(total(:,2))]); 
        xlim(ylimz);
        ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
        set(gca,'xticklabel','');        set(gca,'linewidth',axLineWidth,'XMinorTick','on','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
        plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
        saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.fig']);
        close(f)
        
%         %% 3/3 Projected DOS for H, Si, and Dopant
%         % Projected: orbital = [s py pz px dxy dyz dz2 dxz x2-y2]
%         % Projected: projected{kpoint}{band}(atom,orbital)
%         if      (exist(['color1_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['color2_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['color3_SiHDl_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['color1_down_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['color2_down_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['color3_down_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['DOS_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2) && ...
%                 (exist(['DOS_down_SiHD_',filez(randvar).name(1:end-4),'.mat'],'file') == 2)
%             
%             load(['color1_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color1');
%             load(['color2_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color2');
%             load(['color3_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color3');
%             load(['color1_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
%             load(['color2_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
%             load(['color3_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
%             load(['DOS_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
%             load(['DOS_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
%         else
%             %             s = 1; p = 2; d = 3;
%             %             spdf_vals = {1,2,3,4};
%             spdf_vals = {1,2:4,5:9};
%             for kp = 1:length(projected)
%                 for bnd = 1:length(projected{kp})
%                     for atm = 1:length(atom)
%                         if atm == 1
%                             range{atm} = 1:atom{1};
%                         else
%                             range{atm} = (sum([atom{1:atm-1}])+1):(sum([atom{1:atm-1}])+atom{atm});
%                         end
%                         for spdf = 1:3
%                             A = sum(sum(projected{kp}{bnd}(:,:)));
%                             atomic_projected{kp}{bnd}{atm}(spdf) = sum(sum(projected{kp}{bnd}(range{atm},spdf_vals{spdf})));
%                             percent_atomic_projected{kp}{bnd}{atm}(spdf) = atomic_projected{kp}{bnd}{atm}(spdf)./A;
%                             
%                             A2 = sum(sum(projected_down{kp}{bnd}(:,:)));
%                             atomic_projected_down{kp}{bnd}{atm}(spdf) = sum(sum(projected_down{kp}{bnd}(range{atm},spdf_vals{spdf})));
%                             percent_atomic_projected_down{kp}{bnd}{atm}(spdf) = atomic_projected_down{kp}{bnd}{atm}(spdf)./A2;
%                         end
%                     end
%                 end
%             end
%             
%             color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
%             color1_down = cell(length(band_xaxis{1}),1); color2_down = cell(length(band_xaxis{1}),1); color3_down = cell(length(band_xaxis{1}),1);
%             for bnd = 1:length(band_xaxis{1})
%                 band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
%                 band_plot_down{bnd} = [band_plot_down{bnd} bands_down{1}(bnd)];
%                 color1{bnd} = [color1{bnd} sum(percent_atomic_projected{1}{bnd}{DLoc}(:))];
%                 color2{bnd} = [color2{bnd} sum(percent_atomic_projected{2}{bnd}{HLoc}(:))];
%                 color3{bnd} = [color3{bnd} sum(percent_atomic_projected{3}{bnd}{SiLoc}(:))];
%                 
%                 color1_down{bnd} = [color1_down{bnd} sum(percent_atomic_projected_down{1}{bnd}{DLoc}(:))];
%                 color2_down{bnd} = [color2_down{bnd} sum(percent_atomic_projected_down{2}{bnd}{HLoc}(:))];
%                 color3_down{bnd} = [color3_down{bnd} sum(percent_atomic_projected_down{3}{bnd}{SiLoc}(:))];
%                 for kp = 2:length(band_xaxis)
%                     band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
%                     color1{bnd} = [color1{bnd} sum(percent_atomic_projected{kp}{bnd}{DLoc}(:))];
%                     color2{bnd} = [color2{bnd} sum(percent_atomic_projected{kp}{bnd}{HLoc}(:))];
%                     color3{bnd} = [color3{bnd} sum(percent_atomic_projected{kp}{bnd}{SiLoc}(:))];
%                     
%                     band_plot_down{bnd} = [band_plot_down{bnd} bands_down{kp}(bnd)];
%                     color1_down{bnd} = [color1_down{bnd} sum(percent_atomic_projected_down{kp}{bnd}{DLoc}(:))];
%                     color2_down{bnd} = [color2_down{bnd} sum(percent_atomic_projected_down{kp}{bnd}{HLoc}(:))];
%                     color3_down{bnd} = [color3_down{bnd} sum(percent_atomic_projected_down{kp}{bnd}{SiLoc}(:))];
%                 end
%             end
%             
%             A = [];
%             for ii = 1:ntype
%                 for jj = 1:length(spdf_vals)
%                     A = [];
%                     A2 = [];
%                     for kk = 1:length(range{ii})
%                         if kk == 1
%                             A = partial_up{range{ii}(kk)}(:,spdf_vals{jj}+1);
%                             A2 = partial_down{range{ii}(kk)}(:,spdf_vals{jj}+1);
%                         else
%                             A = A + partial_up{range{ii}(kk)}(:,spdf_vals{jj}+1);
%                             A2 = A2 + partial_down{range{ii}(kk)}(:,spdf_vals{jj}+1);
%                         end
%                     end
%                     DOS{ii}{jj} = A;
%                     DOS_down{ii}{jj} = A2;
%                 end
%             end
%             
%             save(['color1_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color1');
%             save(['color2_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color2');
%             save(['color3_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color3');
%             save(['color1_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color1_down');
%             save(['color2_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color2_down');
%             save(['color3_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'color3_down');
%             save(['DOS_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'DOS');
%             save(['DOS_down_SiHD_',filez(randvar).name(1:end-4),'.mat'], 'DOS_down');
%             
%         end
%         
%         amnt_down = -0.05;
%         amnt_up = -0.04;
%         f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); hold on
%         subplot('Position',[0.12 0.15 0.6 0.82]); hold on;
%         for ii = [1 2]
%             h = plot([band_x_plot{1}(numkpoints*ii) band_x_plot{1}(numkpoints*ii)],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
%         end
%         for bnd = 1:length(band_xaxis{1})
%             for kp = 2:length(band_xaxis)
%                 plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot{bnd}(kp-1)-cbm_vbm(2), band_plot{bnd}(kp)-cbm_vbm(2)],'color',[color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)],'linewidth',plotLineWidth)
%                 plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot_down{bnd}(kp-1)-cbm_vbm(2), band_plot_down{bnd}(kp)-cbm_vbm(2)],'color',[color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)],'linewidth',plotLineWidth)
%             end
%             for kp = 3:4:100
%                 t = text(([band_x_plot{bnd}(kp-2)]), band_plot{bnd}(kp-2)-cbm_vbm(2)-amnt_up,'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
%                 t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2)-amnt_down,'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
%             end
%             for kp = 103:8:150
%                 t = text(([band_x_plot{bnd}(kp-4)]), band_plot{bnd}(kp-4)-cbm_vbm(2)-amnt_up,'△'); t.Color = [color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
%                 t = text(([band_x_plot{bnd}(kp)  ]), band_plot_down{bnd}(kp)-cbm_vbm(2)-amnt_down,'▽'); t.Color = [color1_down{bnd}(kp), color2_down{bnd}(kp), color3_down{bnd}(kp)]; t.FontSize = textSize; t.FontWeight = 'b';
%             end
%         end
%         %         ▲▼△▽˅˄
%         %         plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
%         %         xlim([xlower xupper]);    ylim([actymin actymax]-cbm_vbm(2))
%         %         ylabel('E - E_F (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
%         %         xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
%         %         set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);
%         %         box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
%         %         set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
%         %         xtickz(1) = 0;
%         %         titl = title(strrep(filez(randvar).name(1:end-4),'_',' ')); set(titl,'fontsize',10,'fontweight','b');
%         %         %         xtickz(end+1) = max([band_xaxis{:}]);
%         %         set(gca,'xtick',xtickz,'xticklabel',xlabelz);
%         %
%         plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
%         xlim([xlower xupper]);    ylim([actymin actymax]-cbm_vbm(2))
%         ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
%         xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
%         set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);
%         set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
%         box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
%         set(gca,'xtick',xtickz,'xticklabel',xlabelz);
%         titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
%         
%         % Plotting DOS
%         colorz = [0 0 1; 1 0 0; 0 0.75 0];
%         subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
% 
%         p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
%         plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
% 
%         h3 = area(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{SiLoc}{1} DOS{SiLoc}{2} DOS{SiLoc}{3}],2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = 1;% Si
%         plot(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{SiLoc}{1} DOS{SiLoc}{2} DOS{SiLoc}{3}],2),'color',[0 0 1], 'linewidth',1); % Si
% 
%         h2 = area(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{HLoc}{1} DOS{HLoc}{2} DOS{HLoc}{3}],2)); h2.FaceColor = [0 1 0 ]; h2.FaceAlpha = 1;% H
%         plot(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{HLoc}{1} DOS{HLoc}{2} DOS{HLoc}{3}],2),'color',[0 1 0], 'linewidth',1); % H
% 
%         h1 = area(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{DLoc}{1} DOS{DLoc}{2} DOS{DLoc}{3}],2)); h1.FaceColor = [1 0 0]; h1.FaceAlpha = 1;% Dopant
%         plot(partial_up{1}(:,1)-cbm_vbm(2),sum([DOS{DLoc}{1} DOS{DLoc}{2} DOS{DLoc}{3}],2),'color',[1 0 0], 'linewidth',1); % Dopant
%         
%         
%         
% 
%         
%         p = area(total_down(:,1)-cbm_vbm(2),-total_down(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
%         plot(total_down(:,1)-cbm_vbm(2),-total_down(:,2),'color',[0 0 0], 'linewidth',1);
% 
%         h3 = area(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{SiLoc}{1} DOS_down{SiLoc}{2} DOS_down{SiLoc}{3}],2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = 1;% Si
%         plot(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{SiLoc}{1} DOS_down{SiLoc}{2} DOS_down{SiLoc}{3}],2),'color',[0 0 1], 'linewidth',1.5); % Si
% 
%         h2 = area(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{HLoc}{1} DOS_down{HLoc}{2} DOS_down{HLoc}{3}],2)); h2.FaceColor = [0 1 0 ]; h2.FaceAlpha = 1;% H
%         plot(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{HLoc}{1} DOS_down{HLoc}{2} DOS_down{HLoc}{3}],2),'color',[0 1 0], 'linewidth',1.5); % H
% 
%         h1 = area(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{DLoc}{1} DOS_down{DLoc}{2} DOS_down{DLoc}{3}],2)); h1.FaceColor = [1 0 0]; h1.FaceAlpha = 1;% Dopant
%         plot(partial_down{1}(:,1)-cbm_vbm(2),-sum([DOS_down{DLoc}{1} DOS_down{DLoc}{2} DOS_down{DLoc}{3}],2),'color',[1 0 0], 'linewidth',1.5); % Dopant
%       
% 
%         %         ylim([0 1.05.*max(total(:,2))]);    
%         xlim([actymin actymax]-cbm_vbm(2))
%         leg = legend([p h3 h2 h1],{'Total','Si','H',atomnames{DLoc}},'location','south','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
%         view([90 -90]);
%         ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2); set(gca,'xticklabel','')
%         set(gca,'linewidth',axLineWidth,'XMinorTick','on','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
%         box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
%         
%         saveas(f,[filez(randvar).name(1:end-4) '_AllAtoms' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_AllAtoms' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_AllAtoms' '.fig']);
%         close(f)
%         xtickz = [];
filez(randvar).name(1:end-4)
        clearvars -except filez randvar
        
    end
end


function  outStruct  = xml2struct(input)
%XML2STRUCT converts xml file into a MATLAB structure
% Originally written by W. Falkena, ASTI, TUDelft, 21-08-2010

errorMsg = ['%s is not in a supported format.\n\nInput has to be',...
    ' a java xml object, an xml file, or a string in xml format.'];

% check if input is a java xml object
if isa(input, 'org.apache.xerces.dom.DeferredDocumentImpl') ||...
        isa(input, 'org.apache.xerces.dom.DeferredElementImpl')
    xDoc = input;
else
    try
        if exist(input, 'file') == 2
            xDoc = xmlread(input);
        else
            try
                xDoc = xmlFromString(input);
            catch
                error(errorMsg, inputname(1));
            end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            error(errorMsg, inputname(1));
        else
            rethrow(ME)
        end
    end
end

% parse xDoc into a MATLAB structure
outStruct = parseChildNodes(xDoc);

end

% ----- Local function parseChildNodes -----
function [children, ptext, textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct;
textflag = 'Text';

if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);
    
    for count = 1:numChildNodes
        
        theChild = item(childNodes,count-1);
        [text, name, attr, childs, textflag] = getNodeData(theChild);
        
        if ~strcmp(name,'#text') && ~strcmp(name,'#comment') && ...
                ~strcmp(name,'#cdata_dash_section')
            % XML allows the same elements to be defined multiple times,
            % put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    % put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                % add new element
                children.(name){index} = childs;
                
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name){index}.(textfields{ii}) = ...
                            text.(textfields{ii});
                    end
                end
                if(~isempty(attr))
                    children.(name){index}.('Attributes') = attr;
                end
            else
                % add previously unknown (new) element to the structure
                
                children.(name) = childs;
                
                % add text data ( ptext returned by child node )
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name).(textfields{ii}) = text.(textfields{ii});
                    end
                end
                
                if(~isempty(attr))
                    children.(name).('Attributes') = attr;
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end
            
            % this is the text in an element (i.e., the parentNode)
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    % This is what happens when document is like this:
                    % <element>Text <!--Comment--> More text</element>
                    %
                    % text will be appended to existing ptext
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end
        
    end
end
end

% ----- Local function getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.end


%make sure name is allowed as structure name
name = char(getNodeName(theNode));
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');
name = strrep(name, '_', 'u_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
    attr = [];
end

%parse child nodes
[childs, text, textflag] = parseChildNodes(theNode);

% Get data from any childless nodes. This version is faster than below.
if isempty(fieldnames(childs)) && isempty(fieldnames(text))
    text.(textflag) = char(getTextContent(theNode));
end
end

% ----- Local function parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.
attributes = struct;
if hasAttributes(theNode)
    theAttributes = getAttributes(theNode);
    numAttributes = getLength(theAttributes);
    
    for count = 1:numAttributes
        % Suggestion of Adrian Wanner
        str = char(toString(item(theAttributes,count-1)));
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
    end
end
end

% ----- Local function xmlFromString -----
function xmlroot = xmlFromString(iString)
import org.xml.sax.InputSource
import java.io.*

iSource = InputSource();
iSource.setCharacterStream(StringReader(iString));
xmlroot = xmlread(iSource);
end