clear; clc; close all;
tic
cd '/home/anon/Documents/GraduateSchool/Data/Group IV/DFT/2020_12_16_TwoInsertedOxygen4x4/Config2/pDOS_Organized'

xlabelz = {'K' '\Gamma' 'M' 'K'};
% xlabelz = {'\Gamma' 'X' 'Y' '\Gamma' 'Z' 'R' '\Gamma' 'T' 'U' '\Gamma' 'V'}
ylimz = [-3 4];
atoms_to_probe = {{'Si', 'H', 'O'}}; 
for ii = 1:length(atoms_to_probe)
    for jj = 1:3
        if length(atoms_to_probe{ii}{jj}) == 1
            atoms_to_probe{ii}{jj} = [atoms_to_probe{ii}{jj} ' '];
        end
    end
end
filez = dir('*.xml');
set(0,'defaultAxesFontName', 'Arial');
axLineWidth = 1.5; plotLineWidth = 2; gridTransparency = .15; textSize  = 14; LabelSize = 14; titleFontSize = 24; AxisNumberSize = 18; LegFontSize = 14; colorz = [1 0 0; 0 0 1; 0 1 0]; alpha = 1;

for randvar = 1
    axLineWidth = 1.5; plotLineWidth = 2; gridTransparency = .15; textSize = 14; LabelSize = 14; titleFontSize = 24; AxisNumberSize = 18; LegFontSize = 14; colorz = [1 0 0; 0 0 1; 0 1 0]; alpha = 1;
    direct = pwd;
    if (exist([strrep(direct(72:end-5),'/','_') '_' strrep([atoms_to_probe{1}{:}],' ','') '.svg'])==2)...
            & (exist([strrep(direct(72:end-5),'/','_') '_' strrep([atoms_to_probe{1}{:}],' ','') '.png'])==2)...
        continue
    else
        %% Reading
        %% Misc values
        
        if exist(['data_',filez(randvar).name(1:end-4),'.mat']) == 2
            load(['data_',filez(randvar).name(1:end-4),'.mat'])
        else
            data = xml2struct(filez(randvar).name);
            save(['data_',filez(randvar).name(1:end-4),'.mat'],'data');
        end
        
        natom = str2num(data.modeling.atominfo.atoms.Text); % Ntraumber of atoms in unit cell
        for ii = 1:natom
            atomz(ii) = str2num(data.modeling.atominfo.array{1}.set.rc{ii}.c{2}.Text); % Array of numeric identifiers of atoms
            atomnames_dummy{ii} = [data.modeling.atominfo.array{1}.set.rc{ii}.c{1}.Text];
        end
        
        ntype = length(unique(atomz)); % Number of different typoes of atom s
        
        for ii = 1:ntype
            atom{ii} =  [sum(atomz==ii)];
        end
        atomnames = unique(atomnames_dummy,'stable');
        save(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
        
        %% x-axis for band structure
        
        for ii = 1:length(data.modeling.kpoints.generation.v)
            band_x{ii} = str2num(data.modeling.kpoints.generation.v{ii}.Text);
        end
        numkpoints = str2num(data.modeling.kpoints.generation.i.Text);
        
        % Calculating the recriprocal-lattice-special-points in the basis of the recriprocal lattice
        for ii = 1:3
            recip_basis{ii} =  str2num(data.modeling.structure{2}.crystal.varray{2}.v{ii}.Text);
        end
        
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
        %% Calculating VBM and CBM
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
        cbm_vbm = [cbm vbm];            fermienergy = fermienergy - cbm_vbm(2);
        
        for atom = 1:length(data.modeling.calculation.dos.partial.array.set.set)
            for A = 1:length(data.modeling.calculation.dos.partial.array.set.set{atom}.set.r)
                partial{atom}(A,:) = str2num(data.modeling.calculation.dos.partial.array.set.set{atom}.set.r{A}.Text); % pDOS: partial = [f2]
                total(A,:) = str2num(data.modeling.calculation.dos.total.array.set.set.r{A}.Text);  % Total DOS: total = [energy total integrated]
            end
        end
        load(['atom_',filez(randvar).name(1:end-4),'.mat'],'atom');
        xlower = min([band_xaxis{:}]); xupper = max([band_xaxis{:}]);
        
        for ii = 1:length(atomnames)
            if atomnames{ii} == atoms_to_probe{1}{1}; SiLoc = ii; end
            if atomnames{ii} == atoms_to_probe{1}{2}; HLoc = ii;  end
            if atomnames{ii} == atoms_to_probe{1}{3}; DLoc = ii;  end
        end
        
        %% 3/3 Projected DOS for H, Si, and Dopant
        % Projected: orbital = [s py pz px dxy dyz dz2 dxz x2-y2]
        % Projected: projected{kpoint}{band}(atom,orbital)
        
        band_plot = cell(length(band_xaxis{1}),1);
        band_x_plot = cell(length(band_xaxis{1}),1);
        
        for bnd = 1:length(band_xaxis{1})
            band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{1}(bnd)];
            band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
            for kp = 2:length(band_xaxis)
                band_x_plot{bnd} = [band_x_plot{bnd} band_xaxis{kp}(bnd)];
                band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
            end
        end
        
        s = 1; p = 2; d = 3;
        spdf_vals = {1,2,3,4};
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
        
        set(0,'defaultAxesFontName', 'Times New Roman');
        axLineWidth = 1.5; plotLineWidth = 2; gridTransparency = .15; textSize = 14;
        LabelSize = 14; titleFontSize = 24; AxisNumberSize = 18; LegFontSize = 14;
        f = figure('rend','painters','pos',[10 10 799 795],'visible','on'); clf(f); hold on
        subplot('Position',[0.12 0.15 0.6 0.82]); hold on;
        
        isymin(1) = min([bands{1}]); isymax(1) = max([bands{1}]);
        color1 = cell(length(band_xaxis{1}),1); color2 = cell(length(band_xaxis{1}),1); color3 = cell(length(band_xaxis{1}),1);
        for bnd = 1:length(band_xaxis{1})
            band_plot{bnd} = [band_plot{bnd} bands{1}(bnd)];
            color1{bnd} = [color1{bnd} sum(percent_atomic_projected{1}{bnd}{DLoc}(:))];
            color2{bnd} = [color2{bnd} sum(percent_atomic_projected{2}{bnd}{HLoc}(:))];
            color3{bnd} = [color3{bnd} sum(percent_atomic_projected{3}{bnd}{SiLoc}(:))];
            for kp = 2:length(band_xaxis)
                band_plot{bnd} = [band_plot{bnd} bands{kp}(bnd)];
                color1{bnd} = [color1{bnd} sum(percent_atomic_projected{kp}{bnd}{DLoc}(:))];
                color2{bnd} = [color2{bnd} sum(percent_atomic_projected{kp}{bnd}{HLoc}(:))];
                color3{bnd} = [color3{bnd} sum(percent_atomic_projected{kp}{bnd}{SiLoc}(:))];
                
                if color1{bnd}(kp) > 1;                    color1{bnd}(kp) = 1;                 end
                if color2{bnd}(kp) > 1;                    color2{bnd}(kp) = 1;                 end
                if color3{bnd}(kp) > 1;                    color3{bnd}(kp) = 1;                 end
                
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
        
        xlim([xlower xupper]);    ylim(ylimz)
        ylabel('E - E_{VBM} (eV)','FontSize', LabelSize, 'fontweight','b','color','k');
        xlabel('Wave vector','FontSize', LabelSize, 'fontweight','b','color','k');
        set(findall(gca, 'Type', 'Line'),'LineWidth',plotLineWidth);
        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
        xtickz(1) = 0;
        
        xtickz(1) = 0;        count = 1;
        for ii = [1:length(xlabelz)-2]
            xtickz(count+1) = band_x_plot{1}(numkpoints*ii);
            count = count + 1;
        end
        xtickz(end+1) = max([band_xaxis{:}]);
        for ii = [1:length(xlabelz)-2]
            h = plot([band_x_plot{1}(numkpoints*ii) band_x_plot{1}(numkpoints*ii)],[actymin actymax]-cbm_vbm(2),'color',[0.5 0.5 0.5],'linewidth',1);
            uistack(h,'bottom');
        end
        titl = title(strrep(direct(72:end),'_','\_')); set(titl,'fontsize',10,'fontweight','b');
        plot([min([band_xaxis{:}]) max([band_xaxis{:}])],[fermienergy fermienergy],'--','color',[0.5 0.5 0.5],'linewidth',1)
        set(gca,'xtick',xtickz,'xticklabel',xlabelz);
        
        % Plotting DOS
        colorz = [0 0 1; 1 0 0; 0 0.75 0];
        subplot('Position',[0.75 0.15 0.2 0.82]); hold on;
        p = area(total(:,1)-cbm_vbm(2),total(:,2)); p.FaceColor = [0 0 0]; p.FaceAlpha = 0.75;
        h3 = plot(partial{1}(:,1)-cbm_vbm(2),sum([y_plot{SiLoc}{1} y_plot{SiLoc}{2} y_plot{SiLoc}{3}],2),'color',[0 0 1], 'linewidth',1.5); % Si
        h2 = plot(partial{1}(:,1)-cbm_vbm(2),sum([y_plot{HLoc}{1} y_plot{HLoc}{2} y_plot{HLoc}{3}],2),'color',[0 1 0], 'linewidth',1.5); % H
        h1 = plot(partial{1}(:,1)-cbm_vbm(2),sum([y_plot{DLoc}{1} y_plot{DLoc}{2} y_plot{DLoc}{3}],2),'color',[1 0 0], 'linewidth',1.5); % Dopant
        
        leg = legend([p h3 h2 h1],{'Total',atomnames{SiLoc},atomnames{HLoc},atomnames{DLoc}},'location','best','AutoUpdate','off');  set(leg,'fontsize',LegFontSize,'fontweight','b');
        view([90 -90]);
        ylim([0 1.05.*max(total(:,2))]);
        xlim(ylimz);
        ylabel('DOS (States/eV)','FontSize', LabelSize, 'fontweight','b','color','k');
        set(findall(gca, 'Type', 'Line'),'LineWidth',2); set(gca,'xticklabel','')
        set(gca,'linewidth',axLineWidth,'XMinorTick','off','YMinorTick','on','fontsize',AxisNumberSize,'tickdir','out','MinorGridLineStyle','-');
        box on; grid off; ax = gca;ax.XColor = 'k';ax.YColor = 'k';
        xlim(ylimz)
        
        plot(total(:,1)-cbm_vbm(2),total(:,2),'color',[0 0 0], 'linewidth',1);
        colorz  = [0 0.5 1; 1 0.3 0.3; 0 1 0];
        
        saveas(f,[strrep(direct(72:end-5),'/','_') '_' strrep([atoms_to_probe{1}{:}],' ','') '.svg']);
        saveas(f,[strrep(direct(72:end-5),'/','_') '_' strrep([atoms_to_probe{1}{:}],' ','') '.png']);
%         saveas(f,[strrep(direct(72:end-5),'/','_') '_' strrep([atoms_to_probe{1}{:}],' ','') '.fig']);
        
        close(f)
        xtickz = [];
        
        filez(randvar).name(1:end-4)
        clearvars -except xlabelz filez randvar atoms_to_probe ylimz
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