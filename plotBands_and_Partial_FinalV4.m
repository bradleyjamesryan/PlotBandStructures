close all; clear; clc;
cd '/home/anon/Documents/GraduateSchool/Data/Group IV/DFT/2021_03_08_1x1SiH/pDOS'
xlabelz = {'\Gamma' 'M' 'K' '\Gamma'};
filez = dir('*.xml');
ylimz = [-3 4]

for randvar = 1:length({filez(:).name})
if (exist([filez(randvar).name(1:end-4) '_SiP' '.png'])==2) & (exist([filez(randvar).name(1:end-4) '_SiH' '.png'])==2) & (exist([filez(randvar).name(1:end-4) '_Orbital' '.png'])==2)
continue
else
%% Reading
data = xml2struct(filez(randvar).name);
%% Misc values
if (exist(['natom_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['ntype_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['atom_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['atomnames_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['fermienergy_',filez(randvar).name(1:end-4),'.mat']) == 2)
load(['natom_',filez(randvar).name(1:end-4),'.mat'],'natom');
load(['ntype_',filez(randvar).name(1:end-4),'.mat'],'ntype');
load(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
load(['atomnames_',filez(randvar).name(1:end-4),'.mat'],'atomnames');
load(['fermienergy_',filez(randvar).name(1:end-4),'.mat'], 'fermienergy');
else
fermienergy = str2num(data.modeling.calculation.dos.i.Text);
natom = str2num(data.modeling.atominfo.atoms.Text); % Number of atoms in unit cell
count = 1;
for ii = 1:natom
atomz(ii) = str2num(data.modeling.atominfo.array{1}.set.rc{ii}.c{2}.Text); % Array of numeric identifiers of atoms
atomnames_dummy{ii} = [data.modeling.atominfo.array{1}.set.rc{ii}.c{1}.Text];
if ii > 1 && ~(atomz(ii) == atomz(ii-1))
atom{count} = ii-1;
count = count + 1;
end
end

atom{count} = ii-atom{count-1}; % Number of each atom type
ntype = sum(length(unique(atomz))); % Number of different typoes of atoms
atomnames = unique(atomnames_dummy,'stable');

save(['natom_',filez(randvar).name(1:end-4),'.mat'],'natom');
save(['ntype_',filez(randvar).name(1:end-4),'.mat'],'ntype');
save(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
save(['atomnames_',filez(randvar).name(1:end-4),'.mat'],'atomnames');
save(['fermienergy_',filez(randvar).name(1:end-4),'.mat'],'fermienergy');
end
%% Calculating VBM and CBM


%% x-axis for band structure
if (exist(['nband_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['data_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['fermienergy_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['band_xaxis_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['numkpoints_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['projected_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['partial_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['total_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['bands_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['atom_',filez(randvar).name(1:end-4),'.mat']) == 2)  && (exist(['recip_basis_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['symm_points_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['band_x_',filez(randvar).name(1:end-4),'.mat']) == 2) && (exist(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat']) == 2)
load(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat'], 'cbm_vbm');
load(['band_x_',filez(randvar).name(1:end-4),'.mat'], 'band_x');
load(['nband_',filez(randvar).name(1:end-4),'.mat'], 'nband');
load(['data_',filez(randvar).name(1:end-4),'.mat'], 'data');
load(['band_xaxis_',filez(randvar).name(1:end-4),'.mat'], 'band_xaxis');
load(['numkpoints_',filez(randvar).name(1:end-4),'.mat'], 'numkpoints');
load(['projected_',filez(randvar).name(1:end-4),'.mat'], 'projected');
load(['partial_',filez(randvar).name(1:end-4),'.mat'], 'partial');
load(['total_',filez(randvar).name(1:end-4),'.mat'], 'total');
load(['bands_',filez(randvar).name(1:end-4),'.mat'], 'bands');
load(['recip_basis_',filez(randvar).name(1:end-4),'.mat'],'recip_basis');
load(['symm_points_',filez(randvar).name(1:end-4),'.mat'],'symm_points');
fermienergy = fermienergy - cbm_vbm(2);
else
for ii = 1:length(data.modeling.kpoints.generation.v)
band_x{ii} = str2num(data.modeling.kpoints.generation.v{ii}.Text);
end
numkpoints = str2num(data.modeling.kpoints.generation.i.Text);

% Calculating the recriptocal lattice special points in the basis of the recriprocal lattice
for ii = 1:3
recip_basis{ii} =  str2num(data.modeling.structure{1}.crystal.varray{2}.v{ii}.Text);
end

for ii = 1:length(data.modeling.kpoints.generation.v)
symm_points{ii} = str2num(data.modeling.kpoints.generation.v{ii}.Text);
for jj = 1:3
testing(ii,jj) = symm_points{ii}(1)*recip_basis{1}(jj) + symm_points{ii}(2)*recip_basis{2}(jj) + symm_points{ii}(3)*recip_basis{3}(jj);
end
end

band_xaxis_old = zeros(numkpoints,1);
for ii = 1:length(data.modeling.kpoints.generation.v)-1
new_testing(ii) = sqrt((testing(ii+1,1) - testing(ii,1))^2 + (testing(ii+1,2) - testing(ii,2))^2 + (testing(ii+1,3) - testing(ii,3))^2)./(numkpoints-1);
newer_testing(:,ii) = linspace(0,new_testing(ii).*numkpoints,numkpoints+1);
band_xaxis_old((ii-1).*numkpoints+1:ii.*numkpoints,1) = newer_testing(1:end-1,ii) + band_xaxis_old(end);
end

%% Projected, DOS, and Bands
nband = length(data.modeling.calculation.projected.array.set.set.set{1}.set);
for kpoint = 1:length(data.modeling.calculation.projected.array.set.set.set)
for band = 1:nband
for atom = 1:length(data.modeling.calculation.projected.array.set.set.set{kpoint}.set{band}.r)
projected{kpoint}{band}(atom,:) = str2num(data.modeling.calculation.projected.array.set.set.set{kpoint}.set{band}.r{atom}.Text); % Projected: projected = [s py pz px dxy dyz dz2 dxz x2-y2]
end
end
for A = 1:length(data.modeling.calculation.eigenvalues.array.set.set.set{kpoint}.r)
bands_old{kpoint}{A} = str2num(data.modeling.calculation.eigenvalues.array.set.set.set{kpoint}.r{A}.Text);% Bands: bands_old = [energy, occypancy]
bands{kpoint}(A) = bands_old{kpoint}{A}(1);
end
band_xaxis{kpoint} = linspace(band_xaxis_old(kpoint),band_xaxis_old(kpoint),length(bands{kpoint}));
end

vbm = -1E10;
cbm = 1E10;
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
for atom = 1:length(data.modeling.calculation.dos.partial.array.set.set)
for A= 1:length(data.modeling.calculation.dos.partial.array.set.set{atom}.set.r)
partial{atom}(A,:) = str2num(data.modeling.calculation.dos.partial.array.set.set{atom}.set.r{A}.Text); % pDOS: partial = [f2]
total(A,:) = str2num(data.modeling.calculation.dos.total.array.set.set.r{A}.Text);  % Total DOS: total = [energy total integrated]
end
end
save(['cbm_vbm_',filez(randvar).name(1:end-4),'.mat'], 'cbm_vbm');
save(['band_x_',filez(randvar).name(1:end-4),'.mat'], 'band_x');
save(['nband_',filez(randvar).name(1:end-4),'.mat'], 'nband');
save(['data_',filez(randvar).name(1:end-4),'.mat'], 'data');
save(['band_xaxis_',filez(randvar).name(1:end-4),'.mat'], 'band_xaxis');
save(['numkpoints_',filez(randvar).name(1:end-4),'.mat'], 'numkpoints');
save(['projected_',filez(randvar).name(1:end-4),'.mat'], 'projected');
save(['partial_',filez(randvar).name(1:end-4),'.mat'], 'partial');
save(['total_',filez(randvar).name(1:end-4),'.mat'], 'total');
save(['bands_',filez(randvar).name(1:end-4),'.mat'], 'bands');
save(['recip_basis_',filez(randvar).name(1:end-4),'.mat'],'recip_basis');
save(['symm_points_',filez(randvar).name(1:end-4),'.mat'],'symm_points');
end
load(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');

%% Projected DOS for Si(px), Si(py), Si(pz)
s = 1; p = 2; d = 3;
spdf_vals = {2,3,4};
for kp = 1:length(projected)
for bnd = 1:length(projected{kp})
for atm = 1:length(atom)
if atm == 1
range{atm} = 1:atom{1};
else
range{atm} = (sum([atom{1:atm-1}])+1):(sum([atom{1:atm-1}])+atom{atm});
end
for spdf = 1:3
A = sum(sum(projected{kp}{bnd}(:,:)));
atomic_projected{kp}{bnd}{atm}(spdf) = sum(sum(projected{kp}{bnd}(range{atm},spdf_vals{spdf})));
percent_atomic_projected{kp}{bnd}{atm}(spdf) = atomic_projected{kp}{bnd}{atm}(spdf)./A;
end
end
end
end

set(0,'defaultAxesFontName', 'Times New Roman');
axLineWidth = 1.5; plotLineWidth = 2; gridTransparency = .15; textSize = 14;
LabelSize = 14; titleFontSize = 24; AxisNumberSize = 18; LegFontSize = 14;
f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); hold on
subplot('Position',[0.12 0.15 0.6 0.82]); hold on;

band_x_plot = cell(length(band_xaxis{1}),1);
band_plot = cell(length(band_xaxis{1}),1);
isymin(1) = min([bands{1}]); isymax(1) = max([bands{1}]);

color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
for bnd = 1:length(band_xaxis{1})
band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{1}(bnd)];
band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
color1{bnd} = [color1{bnd} percent_atomic_projected{1}{bnd}{1}(3)];
color2{bnd} = [color2{bnd} percent_atomic_projected{1}{bnd}{1}(1)];
color3{bnd} = [color3{bnd} percent_atomic_projected{1}{bnd}{1}(2)];
for kp = 2:length(band_xaxis)
band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{kp}(bnd)];
band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
color1{bnd} = [color1{bnd} percent_atomic_projected{kp}{bnd}{1}(3)];
color2{bnd} = [color2{bnd} percent_atomic_projected{kp}{bnd}{1}(1)];
color3{bnd} = [color3{bnd} percent_atomic_projected{kp}{bnd}{1}(2)];
plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot{bnd}(kp-1)-cbm_vbm(2), band_plot{bnd}(kp)-cbm_vbm(2)],'color',[color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)],'linewidth',plotLineWidth)
isymin(kp) = min([bands{kp}]); isymax(kp) = max([bands{kp}]);
end
end
actymin = min(isymin);    actymax = max(isymax);

A = [];
for ii = 1:ntype
for jj = 1:length(spdf_vals)
A = [];
for kk = 1:length(range{ii})
if kk == 1
A = partial{range{ii}(kk)}(:,spdf_vals{jj}+1);
else
A = A + partial{range{ii}(kk)}(:,spdf_vals{jj}+1);
end
end
y_plot{ii}{jj} = A;
end
end

xlim([min([band_xaxis{:}]) max([band_xaxis{:}])]);    ylim([ylimz(1) ylimz(2)])
ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);
set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
xtickz(1) = 0;
for ii = 1:length(band_x)-2
h = plot([band_xaxis{numkpoints*ii}(1) band_xaxis{numkpoints*ii}(1) ],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
xtickz(ii+1) = band_xaxis{numkpoints*ii}(1);
uistack(h,'bottom');
end
titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
xtickz(end+1) = max([band_xaxis{:}]);
set(gca,'xtick',xtickz,'xticklabel',xlabelz);

%         Plotting DOS
colorz = [1 0 0; 0 0 1; 0 1 0];
subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;% Total
d = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{3}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(px)
d = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(py)
d = area(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2)); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(pz)
h3 = area(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = 1;% Si(pz)
h2 = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = 0.6;% Si(py)
h1 = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{3}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = 0.6;% Si(px)
leg = legend([p h1 h2 h3],{'Total','Si(px)','Si(py)','Si(pz)'},'location','best','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
view([90 -90]);    ylim([0 1.05.*max(total(:,2))]);    xlim([ylimz(1) ylimz(2)])
ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',2); set(gca,'xticklabel','');        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
plot(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2),'color',[0 0.5 1], 'linewidth',1.5);;
plot(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1},'color',[0 1 0], 'linewidth',1.5);
plot(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{3},'color',[1 0.3 0.3], 'linewidth',1.5);
saveas(f,[filez(randvar).name(1:end-4) '_SiP' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_SiP' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_SiP' '.fig']);
%         close(f)

%% Projected DOS for H(s), Si(s), and Si(p)
s = 1; p = 2; d = 3;
spdf_vals = {1,2:4,5:9};
for kp = 1:length(projected)
for bnd = 1:length(projected{kp})
for atm = 1:length(atom)
if atm == 1
range{atm} = 1:atom{1};
else
range{atm} = (sum([atom{1:atm-1}])+1):(sum([atom{1:atm-1}])+atom{atm});
end
for spdf = 1:3
A = sum(sum(projected{kp}{bnd}(:,:)));
atomic_projected{kp}{bnd}{atm}(spdf) = sum(sum(projected{kp}{bnd}(range{atm},spdf_vals{spdf})));
percent_atomic_projected{kp}{bnd}{atm}(spdf) = atomic_projected{kp}{bnd}{atm}(spdf)./A;
end
end
end
end

f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); hold on
subplot('Position',[0.12 0.15 0.6 0.82]); hold on;

band_x_plot = cell(length(band_xaxis{1}),1);
band_plot = cell(length(band_xaxis{1}),1);
isymin(1) = min([bands{1}]); isymax(1) = max([bands{1}]);

color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
for bnd = 1:length(band_xaxis{1})
band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{1}(bnd)];
band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
color1{bnd} = [color1{bnd} percent_atomic_projected{1}{bnd}{2}(1)];
color2{bnd} = [color2{bnd} percent_atomic_projected{1}{bnd}{1}(1)];
color3{bnd} = [color3{bnd} percent_atomic_projected{1}{bnd}{1}(2)];
for kp = 2:length(band_xaxis)
band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{kp}(bnd)];
band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
color1{bnd} = [color1{bnd} percent_atomic_projected{kp}{bnd}{2}(1)];
color2{bnd} = [color2{bnd} percent_atomic_projected{kp}{bnd}{1}(1)];
color3{bnd} = [color3{bnd} percent_atomic_projected{kp}{bnd}{1}(2)];
plot([band_x_plot{bnd}(kp-1), band_x_plot{bnd}(kp)], [band_plot{bnd}(kp-1)-cbm_vbm(2), band_plot{bnd}(kp)-cbm_vbm(2)],'color',[color1{bnd}(kp), color2{bnd}(kp), color3{bnd}(kp)],'linewidth',plotLineWidth)
isymin(kp) = min([bands{kp}]); isymax(kp) = max([bands{kp}]);
end
end
%         actymin = min(isymin);    actymax = max(isymax);

A = [];
for ii = 1:ntype
for jj = 1:length(spdf_vals)
A = [];
for kk = 1:length(range{ii})
if kk == 1
A = partial{range{ii}(kk)}(:,spdf_vals{jj}+1);
else
A = A + partial{range{ii}(kk)}(:,spdf_vals{jj}+1);
end
end
y_plot{ii}{jj} = A;
end
end

xlim([min([band_xaxis{:}]) max([band_xaxis{:}])]);    ylim([ylimz(1) ylimz(2)])
ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
xtickz = [];
xtickz(1) = 0;
for ii = 1:length(band_x)-2
%             h = plot([band_xaxis{numkpoints*ii}(1) band_xaxis{numkpoints*ii}(1) ],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
xtickz(ii+1) = band_xaxis{numkpoints*ii}(1);
%             uistack(h,'bottom');
end
%         titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'k--','color',[0.5 0.5 0.5],'linewidth',1)
xtickz(end+1) = max([band_xaxis{:}]);
set(gca,'xtick',xtickz,'xticklabel',xlabelz);
%         ylim([-4 8])
set(gca,'yminortick','off')

% Plotting DOS
colorz = [1 0 0; 0 0 1; 0 1 0];
subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
d = area(partial{1}(:,1)-cbm_vbm(2),y_plot{2}{1}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% H(s)
d = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1}); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(s)
d = area(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2)); d.FaceColor = [1 1 1]; d.FaceAlpha = 1;% Si(p)
h3 = area(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2)); h3.FaceColor = [0 0 1]; h3.FaceAlpha = 1;% Si(p)
h2 = area(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1}); h2.FaceColor = [0 1 0]; h2.FaceAlpha = 1;% Si(s)
h1 = area(partial{1}(:,1)-cbm_vbm(2),y_plot{2}{1}); h1.FaceColor = [1 0 0]; h1.FaceAlpha = 1;% H(s)
leg = legend([p h1 h2 h3],{'Total','H(s)','Si(s)','Si(p)'},'location','best','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
view([90 -90]);    ylim([0 1.05.*max(total(:,2))]);    xlim([ylimz(1) ylimz(2)])
%         xlim([-4 8])
ylabel('DOS (eV^{-1})','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',2); set(gca,'xticklabel','');        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
%         plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
%         plot(partial{1}(:,1)-cbm_vbm(2),sum(y_plot{1}{2},2),'color',[0 0.5 1], 'linewidth',1.5);
%         plot(partial{1}(:,1)-cbm_vbm(2),y_plot{1}{1},'color',[0 1 0], 'linewidth',1.5);
%         plot(partial{1}(:,1)-cbm_vbm(2),y_plot{2}{1},'color',[1 0.3 0.3], 'linewidth',1.5);
%         saveas(f,[filez(randvar).name(1:end-4) '_Orbital' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_Orbital' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_Orbital' '.fig']);
%         close(f)

%% Plotting Partial DOS and Bands
% Storing the first atom data into 'partial_atom' so as to preallocate this space
partial_atom{1} = partial{1};
for ii = 2:length(atom)
partial_atom{ii} = partial{atom{ii-1}+1};
end

% Taking the sum of all atoms of a particular type for plotting (summing over all individual orbitals for each similar atom)
count = 2;
for ii = 1:ntype
for jj = 2:atom{ii}
partial_atom{ii}(:,2:end) = partial_atom{ii}(:,2:end) + partial{count}(:,2:end);
count = count + 1;
end
partial_atom{ii}(:,2) = sum(partial_atom{ii}(:,2:end),2);
end

total_sum = zeros(length(partial_atom{1}(:,1)),1);
for ii = 1:ntype
total_sum = total_sum + partial_atom{ii}(:,2);
end

norm_percent{2}= [];    norm_percent{3} = [];
for ii = 1:ntype
percent{ii}(:,1) = partial_atom{1}(:,1); % storing the energies
percent{ii}(:,2) = partial_atom{ii}(:,2)./total_sum;
norm_percent{ii}(:,1) = partial_atom{1}(:,1); % storing the energies
norm_percent{ii}(:,2) = (percent{ii}(:,2) - min(percent{ii}(:,2)))./(max(percent{ii}(:,2)) - min(percent{ii}(:,2)));
end

if isempty(norm_percent{2})
norm_percent{2} = zeros(size(norm_percent{1}));
norm_percent{3} = zeros(size(norm_percent{1}));
elseif isempty(norm_percent{3})
norm_percent{3} = zeros(size(norm_percent{1}));
end

band_plot2 = cell(nband,1);        x_plot2 = cell(nband,1);        x_plot2 = [];
for bnd = 1:nband
for kp = 1:length(projected)
band_plot2{bnd} = [band_plot2{bnd} bands{kp}(bnd)-cbm_vbm(2)];
if bnd ==1
x_plot2(kp) = band_xaxis{kp}(1);
end
end
isymin(bnd) = min([bands{bnd}]);        isymax(bnd) = max([bands{bnd}]);
end
%         actymin = min(isymin);    actymax = max(isymax);

f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); hold on
subplot('Position',[0.12 0.15 0.6 0.82]); hold on;
for bnd = 1:nband
for kp = 2:length(bands)
[loc, ~]=find((abs(band_plot2{bnd}(kp)  -  (partial_atom{1}(:,1)-cbm_vbm(2))))<0.1);
[~, loc2]=min(abs(partial_atom{1}(loc,1)-cbm_vbm(2)  -  band_plot2{bnd}(kp)));
plot([x_plot2(kp-1) x_plot2(kp)], [band_plot2{bnd}(kp-1) band_plot2{bnd}(kp)], 'color',[norm_percent{3}(loc(loc2),2) norm_percent{2}(loc(loc2),2) norm_percent{1}(loc(loc2),2)])
end
end
%         actymin = -10.558;
%         actymax = 5.5;

% Plotting grey vertical lines (still in bands plot)
xlim([min([band_xaxis{:}]) max([band_xaxis{:}])]);    ylim([ylimz(1) ylimz(2)])
ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
xtickz = []; xtickz(1) = 0;
for ii = 1:length(band_x)-2
h = plot([band_xaxis{numkpoints*ii}(1) band_xaxis{numkpoints*ii}(1) ],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
xtickz(ii+1) = band_xaxis{numkpoints*ii}(1);
uistack(h,'bottom');
end
titl = title(filez(randvar).name(1:end-4)); set(titl,'fontsize',10,'fontweight','b');
plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
xtickz(end+1) = max([band_xaxis{:}]);
set(gca,'xtick',xtickz,'xticklabel',xlabelz);

% Plotting DOS
colorz = [0 0 1; 1 0 0; 0 0.7 0];
subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
h = area(total(:,1)-cbm_vbm(2),total(:,2)); h.FaceColor = [0 0 0]; h.FaceAlpha = 0.75;
p = area(partial_atom{1}(:,1)-cbm_vbm(2),partial_atom{1}(:,2)); p.FaceColor = [1 1 1];
p = area(partial_atom{2}(:,1)-cbm_vbm(2),partial_atom{2}(:,2)); p.FaceColor = [1 1 1];
h1 = area(partial_atom{1}(:,1)-cbm_vbm(2),partial_atom{1}(:,2)); h1.FaceColor = colorz(1,:); h1.FaceAlpha = 0.75;
h2 = area(partial_atom{2}(:,1)-cbm_vbm(2),partial_atom{2}(:,2)); h2.FaceColor = colorz(3,:); h2.FaceAlpha = 0.75;
view([90 -90]);    ylim([0 1.05.*max(total(:,2))]);    xlim([ylimz(1) ylimz(2)])
ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
set(findall(gca, 'Type', 'Line'),'LineWidth',2); set(gca,'xticklabel','');        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');        leg = legend([h h1 h2],['Total',atomnames],'location','best','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
plot(partial_atom{1}(:,1)-cbm_vbm(2),partial_atom{1}(:,2),'color', [0 0.5 1],'linewidth',1.5);
plot(partial_atom{2}(:,1)-cbm_vbm(2),partial_atom{2}(:,2), 'color', [0 0.7 0],'linewidth',1.5);
saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.svg']);        saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.png']);        saveas(f,[filez(randvar).name(1:end-4) '_SiH' '.fig']);
close(f);

filez(randvar).name(1:end-4)
%         clearvars -except xlabelz filez randvar
end
end


function  outStruct  = xml2struct(input)
%XML2STRUCT converts xml file into a MATLAB structure
%
% outStruct = xml2struct2(input)
% 
% xml2struct2 takes either a java xml object, an xml file, or a string in
% xml format as input and returns a parsed xml tree in structure. 
% 
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Originally written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increase by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Chao-Yuan Yeh, August 2016

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
% Create structure of node info.

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

% This alterative to the above 'if' block will also work but very slowly.
% if any(strcmp(methods(theNode),'getData'))
%   text.(textflag) = char(getData(theNode));
% end

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

