function BackgroundSurvey_Sample
filepath='H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
filepath2='H:\Images\2013-12-12_BackgroundCharacterization\';
filetag='_Cy3_';
shot1='1_3_1';
shot2='4_2_1';

filename1=[filepath,'Raw\',shot1,filetag,'stain.tif'];
a=single(imread(filename1));
%highthresh=callthresh(a,3);
highthresh=10000;
a=removesmears(a,highthresh);
filename2=[filepath2,'Raw\',shot2,filetag,'stain.tif'];
b=single(imread(filename2));
b=removesmears(b,highthresh);

%%% alternative calculation for self shift %%%%%%%%%%
% [height,width]=size(a);
% heightshift=round(height/10); widthshift=round(width/10);
% aorg=a;
% a=aorg(1+heightshift:end,1+widthshift:end);
% b=aorg(1:end-heightshift,1:end-widthshift);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aorg=a; borg=b;
a=a(:); b=b(:);
nanidx=isnan(a) | isnan(b);
a(nanidx)=[]; b(nanidx)=[];
%as=a(:); meana=mean(as); stda=std(as); norma=(as-meana)/stda; bs=b(:); meanb=mean(bs); stdb=std(bs); normb=(bs-meanb)/stdb; axb=norma.*normb; c=sum(axb(:))/(numel(as));
as=a-mean(a); bs=b-mean(b); c=dot(as/norm(as),bs/norm(bs));
fprintf('c = %8.4f\n',c);

figure; imagesc(aorg); colorbar; %default [600,750]
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end