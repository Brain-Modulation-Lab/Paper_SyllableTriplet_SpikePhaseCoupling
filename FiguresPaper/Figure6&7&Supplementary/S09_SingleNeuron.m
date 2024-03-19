clc
clear all
close all

% set paths
PATH_GROUPLEVEL = 'W:\Users\MV1019\PhaseLocking\groupanalyses\Output\Results_new\default\group-level';
PATH_FRDATA = fullfile(PATH_GROUPLEVEL,'FRDATA_group.mat');
PATH_FRSTATS = fullfile(PATH_GROUPLEVEL,'FRSTATS_group.mat');
PATH_ANALYSIS = fullfile('W:\Users\MV1019\PhaseLocking\Supplementary_Analysis');
PATH_FIGURES = fullfile(PATH_ANALYSIS,'Figures');
PATH_RESULTS = fullfile(PATH_ANALYSIS,'Results');
% load firing rate data & stats
load(PATH_FRDATA);
load(PATH_FRSTATS);

%%

unit_list = [188 31 21 34]; % -, + , -+, x neurons 
subjid_list = {DATA(unit_list).SubjectID};
track_list = {DATA(unit_list).electrode};
session_list = [DATA(unit_list).session];

plot_example_neuron(DATA,unit_list);
figname = fullfile(PATH_FIGURES,'example_singleneurons_firing_rate');
saveas(gcf,[figname,'.fig'])
saveas(gcf,[figname,'.png'])
saveas(gcf,[figname,'.pdf'])
saveas(gcf,[figname,'.svg'])

% get MNI coordinates of these neurons
load(fullfile(PATH_GROUPLEVEL,'MNI_group.mat'));
MNIcoords = get_MNIcoords(MNI_all,subjid_list,session_list,track_list);

%% plot MNI coords on STN atlas using LeadDBS;

atlas = 'DISTAL Minimal (Ewert 2017)';
colors = [0 0 1; 1 0 0; 0 .5 0];
load([ea_space([],'atlases'),atlas,filesep,'atlas_index.mat']); % manually load definition of DISTAL atlas.

plot_MERLocations_LeadDBS(MNIcoords,atlases,colors);


function plot_MERLocations_LeadDBS(MNIcoords,atlases,colors)

ea_mnifigure(); % open up Elvis viewer
lSTN=atlases.roi{1, 2}.fv; % do the same for the left STN.
lSTN=reducepatch(lSTN,0.5);
patch('Faces',lSTN.faces,'Vertices',lSTN.vertices,'facecolor','#EDB120','facealpha',.25,'edgecolor','w');
hold on
s = plot_colored_spheres(MNIcoords,[1 2 3],0.3, colors);
% for loci = 1 : numel(colors)
%     scatter3(MNIcoords(loci,1),MNIcoords(loci,2),MNIcoords(loci,3),175,colors(loci),'filled');
% end
ea_setplanes(nan,nan,nan); % set planes to a nice view.

% now define a specific camera view:
v.az= 126.0430;
v.el= -34.9893;
v.camva= 0.4648;
v.camup= [0 0 1];
v.camproj= 'orthographic';
v.camtarget= [-1.2419 -23.0911 3.9242];
v.campos= [-1.0973e+03 1.0671e+03 -1.0542e+03];
ea_view(v); % apply view.

end

% list functions

function MNIcoords = get_MNIcoords(MNI_all,subjid_list,session_list,track_list)
nNeurons = numel(session_list);
MNIcoords = nan(nNeurons,3);
for unit_i = 1 :nNeurons
    subjcoords = MNI_all(contains({MNI_all.subj_id},subjid_list{unit_i})).S;
    MNIcoords(unit_i,:) = subjcoords{subjcoords.session_id == session_list(unit_i)  & contains(subjcoords.electrode, track_list{unit_i}) ,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}};
end

end


function plot_example_neuron(DATA,unit_list)
nNeurons = numel(unit_list);
figure('renderer','painters','position',[300 100 600*(nNeurons) 500])
tiledlayout(3,nNeurons)

gcp = 1;
for unit_i = unit_list
    UnitType =DATA(unit_i).RecType;
    UnitGrade =DATA(unit_i).grade;
    
    ST = DATA(unit_i).SpeechOnset.DD; % spike train
    IFRmean = DATA(unit_i).SpeechOnset.meanifr;
    [nTrials,nSamples] = size(ST);
    IFRsem = DATA(unit_i).SpeechOnset.stdifr;
    SignFlag = [DATA(unit_i).SpeechOnset.sig_excit; DATA(unit_i).SpeechOnset.sig_inhib];



    IFRbase = DATA(unit_i).SpeechOnset.basemeanIFR;

    events = [DATA(unit_i).coding_events.syl1_onset - DATA(unit_i).coding_events.syl1_onset DATA(unit_i).coding_events.syl1_offset - DATA(unit_i).coding_events.syl1_onset , ...
        DATA(unit_i).coding_events.syl2_onset - DATA(unit_i).coding_events.syl1_onset DATA(unit_i).coding_events.syl2_offset - DATA(unit_i).coding_events.syl1_onset, ...
        DATA(unit_i).coding_events.syl3_onset - DATA(unit_i).coding_events.syl1_onset DATA(unit_i).coding_events.syl3_offset - DATA(unit_i).coding_events.syl1_onset];
    SpeechDur = DATA(unit_i).SpDur;
    [~,idx_ordtrials] = sort(SpeechDur);
    
    FRpctchange = 100*((IFRmean - IFRsem) /IFRbase - 1);
    FRzscoretime = (IFRmean - IFRbase)/std(mean(DATA(ui).SpeechOnset.IFRbase,1),[],2);
    
    Time = linspace(DATA(unit_i).SpeechOnset.respInterval(1),DATA(unit_i).SpeechOnset.respInterval(2),nSamples);
    % raster plot


    % raster plot - order trials by duration of speech prodcution epoch
    nexttile(gcp,[1 1])
    title(sprintf('%s-%s <IFR_{baseline}> = %1.2f spk/s',UnitType, UnitGrade, IFRbase))
    hold on
    ST_ord = ST(idx_ordtrials,:);
    events_ord = events(idx_ordtrials,:);
    for trial_i = 1 : nTrials
        scatter(Time(ST_ord(trial_i,:) == 1), trial_i*ones(1,sum(ST_ord(trial_i,:)==1)),1,'k','filled')
        plot([events_ord(trial_i,1) events_ord(trial_i,2)],[trial_i,trial_i],'color',[1 0 0 .3])
        plot([events_ord(trial_i,3) events_ord(trial_i,4)],[trial_i,trial_i],'color',[1 0 0 .3])
        plot([events_ord(trial_i,5) events_ord(trial_i,6)],[trial_i,trial_i],'color',[1 0 0 .3])
    end
    ylabel('ordered trial [#]')
    axis("tight")

    % firing rate change plot w.r.t. baseline
    nexttile(gcp + nNeurons,[2 1])
    toPlot = FRpctchange;
    xShaded = [Time fliplr(Time)];
    yShaded = [toPlot  fliplr(toPlot )];
    fill(xShaded, yShaded,'k','FaceAlpha',.3);
    hold on
    plot(Time, toPlot ...
        ,'k','linewidth',1.3)
    yline(0,'--')
    % create filling stats
    for sign_type = 1 : 2 % + and -
        sign_bin = SignFlag(sign_type,:) == 1;
        if any(sign_bin)
            segments = bwconncomp(sign_bin);
            for seg = 1 : segments.NumObjects
                segment = segments.PixelIdxList{seg};
                xFill = [Time(segment) fliplr(Time(segment))];
                switch sign_type
                    case 1
                        yFill = [zeros(1,numel(segment)) fliplr(100*((IFRmean(segment) - IFRsem(segment)) /IFRbase - 1) ) ];
                        fill(xFill, yFill,'r','FaceAlpha',.5)
                    case 2
                        yFill = [100*((IFRmean(segment) + IFRsem(segment)) /IFRbase - 1)  zeros(1,numel(segment))];
                        fill(xFill, yFill,'b','FaceAlpha',.5)
                end
            end

        end
    end
    box off
    %
    xlabel('Time [s] ')
    ylabel('IFR change [%]')
    axis("tight")
    gcp = gcp + 1;
end
end