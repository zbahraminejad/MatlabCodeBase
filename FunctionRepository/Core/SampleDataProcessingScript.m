
path = '/Users/Zahrabahrami/Documents/20141010-Zahra  Fix after live imaging_Plate_459/';
load('/Users/Zahrabahrami/Documents/20141010-Zahra  Fix after live imaging_Plate_459/20150127_Zahra.mat');

site_num = 1;

eyfp = handles.Measurements.Nuclei.Intensity_MeanIntensity_CorrEYFP;
texas = handles.Measurements.Nuclei.Intensity_MeanIntensity_CorrTexasRed;

wellnames = handles.Measurements.Image.FileName_EYFP;
for i=1:length(wellnames); Wellnames{i}=wellnames{i}(39:43); end
%%
[meanTexas,distTexas] = cleanData(texas,site_num);
[meaneYfp,disteYfp] = cleanData(eyfp,site_num); %first input is "data" for clean fxn; second input is site number set above (# images/well)

row = 4;col = 10;
reShapePlate = @(data)reshape(data',col,row); %defines how to reshape data below (mini-function);

distmcitrine = reShapePlate(disteYfp)';
meanmcitrine = reShapePlate(meaneYfp)';
distPPARG = reShapePlate(distTexas)';
meanPPARG = reShapePlate(meanTexas)';
WellNames = reShapePlate(Wellnames)'
%%
%Clone 1 
conditions={
    %'title',row,col;
   'No stimulus',1,2;
   'Rosi',2,2;
   'Rosi+DMI',3,2;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,2;
   'Rosi',2,2;
   'Rosi+DMI',3,2;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 mCitrine');
filename=[path '20141010_clone1mCitrine.eps'];
print('-depsc',filename)
%%
%Clone 2
conditions={
    %'title',row,col;
   'No stimulus',1,3;
   'Rosi',2,3;
   'Rosi+DMI',3,3;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone2 Antibody Staining');
filename=[path '20141010_clone2Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,3;
   'Rosi',2,3;
   'Rosi+DMI',3,3;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone2 mCitrine');
filename=[path '20141010_clone2mCitrine.eps'];
print('-depsc',filename)
%%
%clone3
conditions={
    %'title',row,col;
   'No stimulus',1,4;
   'Rosi',2,4;
   'Rosi+DMI',3,4;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone3 Antibody Staining');
filename=[path '20141010_clone3Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,4;
   'Rosi',2,4;
   'Rosi+DMI',3,4;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone3 mCitrine');
filename=[path '20141010_clone3mCitrine.eps'];
print('-depsc',filename)
%%
%clone 4
conditions={
    %'title',row,col;
   'No stimulus',1,5;
   'Rosi',2,5;
   'Rosi+DMI',3,5;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone4 Antibody Staining');
filename=[path '20141010_clone4Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,5;
   'Rosi',2,5;
   'Rosi+DMI',3,5;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone4 mCitrine');
filename=[path '20141010_clone4mCitrine.eps'];
print('-depsc',filename)
%%
%clone5
conditions={
    %'title',row,col;
   'No stimulus',1,6;
   'Rosi',2,6;
   'Rosi+DMI',3,6;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,6;
   'Rosi',2,6;
   'Rosi+DMI',3,6;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone5 mCitrine');
filename=[path '20141010_clone5mCitrine.eps'];
print('-depsc',filename)

%%
%clone6
conditions={
    %'title',row,col;
   'No stimulus',1,7;
   'Rosi',2,7;
   'Rosi+DMI',3,7;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,7;
   'Rosi',2,7;
   'Rosi+DMI',3,7;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone6 mCitrine');
filename=[path '20141010_clone6mCitrine.eps'];
print('-depsc',filename)

%%
%clone7
conditions={
    %'title',row,col;
   'No stimulus',1,8;
   'Rosi',2,8;
   'Rosi+DMI',3,8;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,6;
   'Rosi',2,6;
   'Rosi+DMI',3,6;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone8 mCitrine');
filename=[path '20141010_clone8mCitrine.eps'];
print('-depsc',filename)

%%
%clone8
conditions={
    %'title',row,col;
   'No stimulus',1,9;
   'Rosi',2,9;
   'Rosi+DMI',3,9;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,9;
   'Rosi',2,9;
   'Rosi+DMI',3,9;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone8 mCitrine');
filename=[path '20141010_clone8mCitrine.eps'];
print('-depsc',filename)
%%
%clone9
conditions={
    %'title',row,col;
   'No stimulus',1,10;
   'Rosi',2,10;
   'Rosi+DMI',3,10;
   };

   drawHistTaka(conditions,distPPARG,1e-3,1.5e-2)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone1 Antibody Staining');
filename=[path '20141010_clone1Texas.eps'];
print('-depsc',filename)

figure()
conditions={
    %'title',row,col;
   'No stimulus',1,10;
   'Rosi',2,10;
   'Rosi+DMI',3,10;
   };

   drawHistTaka(conditions,distmcitrine,5e-4,4e-3)%conditions,distData,binMin,binMax,binNum,'log'
    set(gca, 'xscale', 'lin', 'fontsize', 18)
    ylabel('Percent'); xlabel('PPARG Fluoresence (Mean)'); box on; title('Clone10 mCitrine');
filename=[path '20141010_clone10mCitrine.eps'];
print('-depsc',filename)
%%
thresholdPPARG=4.7e-3; %Fluorescence value cutoff for PPARG+
   for r=1:row; for c=1:col;      
        a=distPPARG{r,c}>thresholdPPARG;
        percentPPARG(r,c)=sum(a)/length(distPPARG{r,c})*100
       end; end;
%%



%%
pulses=[4 3 2 1];

figure(); set(gcf, 'color', 'w');
hold on
errorbar(pulses, mean(percentPPARG(1:4,7:9),2), std(percentPPARG(1:4,7:9),1,2)/sqrt(3), '-or', 'linewidth', 3)
errorbar(pulses, mean(percentPPARG(1:4,1:3),2), std(percentPPARG(1:4,1:3),1,2)/sqrt(3), '-om', 'linewidth', 3)
errorbar(pulses, mean(percentPPARG(1:4,4:6),2), std(percentPPARG(1:4,4:6),1,2)/sqrt(3), '-og', 'linewidth', 3)
errorbar(pulses, mean(percentPPARG(1:4,10:12),2), std(percentPPARG(1:4,10:12),1,2)/sqrt(3), '-ob', 'linewidth', 3)

set(gca, 'yscale', 'lin', 'fontsize', 22, 'ytick', [20 40 60 80 100])
axis([0 5 0 100]); box on; xlabel('Number of ~12h Pulses'); ylabel('% PPARG+ at Day 4');
%legend('DMI Pulses', 'Dex Pulses', 'IBMX Pulses', 'Rosi Pulses')
filename=[path '20141021_PPARG_PercentPulses.eps'];
print('-depsc',filename)
