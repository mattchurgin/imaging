% Matt Churgin, August 2018
clear all
close all

homeDir='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR';
cd(homeDir)

% lookup table of gloms identified in gh146-gal4
% from V. Grabe et al. antenna lobe atlas
gh146glomsR{1}='Or69a';
gh146glomsR{2}='0r67d';
gh146glomsR{3}='Or56a';
gh146glomsR{4}='Or43a';
gh146glomsR{5}='Or19a';
gh146glomsR{6}='Or13a';
gh146glomsR{7}='Or83c';
gh146glomsR{8}='Or10a';
gh146glomsR{9}='Or65a';
gh146glomsR{10}='Or49a';
gh146glomsR{11}='Or7a';
gh146glomsR{12}='Or42b';
gh146glomsR{13}='Or22a';
gh146glomsR{14}='Or47a';
gh146glomsR{15}='Or59b';
gh146glomsR{16}='Or33b';
gh146glomsR{17}='Or67a';
gh146glomsR{18}='Ir64a';
gh146glomsR{19}='Or88a';
gh146glomsR{20}='Or47b';
gh146glomsR{21}='Or92a';
gh146glomsR{22}='Or67b';
gh146glomsR{23}='Or85d';
gh146glomsR{24}='Or49b';
gh146glomsR{25}='Or82a';
gh146glomsR{26}='Or46a';
gh146glomsR{27}='Or33c';
gh146glomsR{28}='Or71a';
gh146glomsR{29}='Ir84a';
gh146glomsR{30}='Ir31a';
gh146glomsR{31}='Ir92a';
gh146glomsR{32}='Or43b';
gh146glomsR{33}='Or9a';
gh146glomsR{34}='Ir76a';
gh146glomsR{35}='Or42a';
gh146glomsR{36}='Or59c';
% adding three additional gloms listed as hit by GH146 in Jeanne et al.
% 2018
gh146glomsR{37}='Ir75d';
gh146glomsR{38}='Unknown';
gh146glomsR{39}='Or23a';

gh146gloms{1}='D';
gh146gloms{2}='DA1';
gh146gloms{3}='DA2';
gh146gloms{4}='DA4l';
gh146gloms{5}='DC1';
gh146gloms{6}='DC2';
gh146gloms{7}='DC3';
gh146gloms{8}='DL1';
gh146gloms{9}='DL3';
gh146gloms{10}='DL4';
gh146gloms{11}='DL5';
gh146gloms{12}='DM1';
gh146gloms{13}='DM2';
gh146gloms{14}='DM3';
gh146gloms{15}='DM4';
gh146gloms{16}='DM5';
gh146gloms{17}='DM6';
gh146gloms{18}='DP1m';
gh146gloms{19}='VA1d';
gh146gloms{20}='VA1v';
gh146gloms{21}='VA2';
gh146gloms{22}='VA3';
gh146gloms{23}='VA4';
gh146gloms{24}='VA5';
gh146gloms{25}='VA6';
gh146gloms{26}='VA7l';
gh146gloms{27}='VC1';
gh146gloms{28}='VC2';
gh146gloms{29}='VL2a';
gh146gloms{30}='VL2p';
gh146gloms{31}='VM1';
gh146gloms{32}='VM2';
gh146gloms{33}='VM3';
gh146gloms{34}='VM4';
gh146gloms{35}='VM7d';
gh146gloms{36}='VM7v';
% adding three additional gloms listed as hit by GH146 in Jeanne et al.
% 2018
gh146gloms{37}='VL1';
gh146gloms{38}='VA7m';
gh146gloms{39}='DA3';

% lookup table of gloms identified in orco-gal4
% from V. Grabe et al. antenna lobe atlas
orcoglomsR{1}='Or69a';
orcogloms{1}='D';
orcoglomsR{2}='Or67d';
orcogloms{2}='DA1';
orcoglomsR{3}='Or56a';
orcogloms{3}='DA2';
orcoglomsR{4}='Or23a';
orcogloms{4}='DA3';
orcoglomsR{5}='Or43a';
orcogloms{5}='DA4l';
orcoglomsR{6}='Or2a';
orcogloms{6}='DA4m';
orcoglomsR{7}='Or19a';
orcogloms{7}='DC1';
orcoglomsR{8}='Or13a';
orcogloms{8}='DC2';
orcoglomsR{9}='Or83c';
orcogloms{9}='DC3';
orcoglomsR{10}='Or10a';
orcogloms{10}='DL1';
orcoglomsR{11}='Or65a';
orcogloms{11}='DL3';
orcoglomsR{12}='Or49a';
orcogloms{12}='DL4';
orcoglomsR{13}='Or7a';
orcogloms{13}='DL5';
orcoglomsR{14}='Or42b';
orcogloms{14}='DM1';
orcoglomsR{15}='Or22a';
orcogloms{15}='DM2';
orcoglomsR{16}='Or47a';
orcogloms{16}='DM3';
orcoglomsR{17}='Or59b';
orcogloms{17}='DM4';
orcoglomsR{18}='Or33b';
orcogloms{18}='DM5';
orcoglomsR{19}='Or67a';
orcogloms{19}='DM6';
orcoglomsR{20}='Gr21a';
orcogloms{20}='V';
orcoglomsR{21}='Or88a';
orcogloms{21}='VA1d';
orcoglomsR{22}='Or47b';
orcogloms{22}='VA1v';
orcoglomsR{23}='Or92a';
orcogloms{23}='VA2';
orcoglomsR{24}='Or67b';
orcogloms{24}='VA3';
orcoglomsR{25}='Or85d';
orcogloms{25}='VA4';
orcoglomsR{26}='Or49b';
orcogloms{26}='VA5';
orcoglomsR{27}='Or82a';
orcogloms{27}='VA6';
orcoglomsR{28}='Or46aA';
orcogloms{28}='VA7l';
orcoglomsR{29}='Unknown';
orcogloms{29}='VA7m';
orcoglomsR{30}='Or33c';
orcogloms{30}='VC1';
orcoglomsR{31}='Or71a';
orcogloms{31}='VC2';
orcoglomsR{32}='Or35a';
orcogloms{32}='VC3';
orcoglomsR{33}='Or67c';
orcogloms{33}='VC4';
orcoglomsR{34}='Or43b';
orcogloms{34}='VM2';
orcoglomsR{35}='Or9a';
orcogloms{35}='VM3';
orcoglomsR{36}='Or85b';
orcogloms{36}='VM5d';
orcoglomsR{37}='Or98a';
orcogloms{37}='VM5v';
orcoglomsR{38}='Or42a';
orcogloms{38}='VM7d';
orcoglomsR{39}='Or59c';
orcogloms{39}='VM7v';

% load extracted glomeruli centroids

filesToLoad{1}='extractedALcentroids.csv';
filesToLoad{2}='fullExtractedALdata.csv';

% load AL centroids
fid = fopen(filesToLoad{1});
currdata = textscan(fid, '%s','Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
fclose(fid);

for j=7:5:length(currdata{1})
    glomLocationNames(str2num(str2num(currdata{1}{j-1})))=str2num(currdata{1}{j});
    xc{str2num(str2num(currdata{1}{j-1}))}=str2num(currdata{1}{j+1});
    yc{str2num(str2num(currdata{1}{j-1}))}=str2num(currdata{1}{j+2});
    zc{str2num(str2num(currdata{1}{j-1}))}=str2num(currdata{1}{j+3});
end

% load all AL data
fid = fopen(filesToLoad{2});
currdata = textscan(fid, '%s','Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
fclose(fid);

tempi=1;
for j=6:5:length(currdata{1})
    allALdataGlomname{tempi}=str2num(currdata{1}{j+1});
    allALdataX(tempi)=str2num(currdata{1}{j+2});
    allALdataY(tempi)=str2num(currdata{1}{j+3});
    allALdataZ(tempi)=str2num(currdata{1}{j+4});
    tempi=tempi+1;
end

% save each AL's data into a unique cell array if it is a gh146 glom
xAllDatagh146=cell(1,length(gh146gloms));
yAllDatagh146=cell(1,length(gh146gloms));
zAllDatagh146=cell(1,length(gh146gloms));
% remove glomeruli not tagged by gh146
for i=1:length(allALdataGlomname)
    for j=1:length(gh146gloms)
        if strcmp(allALdataGlomname{i},gh146gloms{j})
            xAllDatagh146{j}=[xAllDatagh146{j} allALdataX(i)];
            yAllDatagh146{j}=[yAllDatagh146{j} allALdataY(i)];
            zAllDatagh146{j}=[zAllDatagh146{j} allALdataZ(i)];
            break
        end
    end
    if mod(i,100)==0
       display(['processed element ' num2str(i) ' of ' num2str(length(allALdataGlomname))]) 
    end
end

% save odor panel data in unique folder (in case other odor panels exist)
dname=uigetdir;
cd(dname)
try
    splits=strfind(dname,'/');
    savename=[dname((splits(end)+1):end) '_DoORData'];
catch
    splits=strfind(dname,'\');
    savename=[dname((splits(end)+1):end) '_DoORData'];
end

currFiles=dir(pwd);

dummyi=1;
for j=1:length(currFiles)
    if strfind(currFiles(j).name,'csv')
        filesToLoad{dummyi}=currFiles(j).name;
        dummyi=dummyi+1;
    end
end

numGlom=78;
% initialize odor matrix
odorResponses=NaN*zeros(numGlom,dummyi-1);
glomNames=cell(1,numGlom);

for i=1:dummyi-1
    fid = fopen(filesToLoad{i});
    currdata = textscan(fid, '%s','Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
    fclose(fid);
    for j=5:4:length(currdata{1})
        try
            odorResponses(str2num(str2num(currdata{1}{j})),i)=str2num(currdata{1}{j+3});
            glomNames{str2num(str2num(currdata{1}{j}))}=str2num(currdata{1}{j+1});
        catch
            odorResponses(str2num(currdata{1}{j}),i)=str2num(currdata{1}{j+3});
            glomNames{str2num(currdata{1}{j})}=currdata{1}{j+1};
        end
    end
end

publishedOR.rawOdorResponses=odorResponses;

% set inhibitory responses to zero since we are comparing to calcium data
% which doesn't exhibit much inhibitory response
odorResponses(odorResponses<0)=0;

glomLocationNamesgh146=cell(1,length(gh146gloms));
xcgh146=NaN*zeros(1,length(gh146gloms));
ycgh146=NaN*zeros(1,length(gh146gloms));
zcgh146=NaN*zeros(1,length(gh146gloms));
% remove glomeruli not tagged by gh146
for i=1:length(glomLocationNames)
    for j=1:length(gh146gloms)
        if strcmp(glomLocationNames{i},gh146gloms{j})
            glomLocationNamesgh146{j}=glomLocationNames{i};
            xcgh146(j)=xc{i};
            ycgh146(j)=yc{i};
            zcgh146(j)=zc{i};
            break
        end
    end
end

glomNamesgh146=cell(1,length(gh146gloms));
odorResponsesgh146=NaN*zeros(length(gh146gloms),size(odorResponses,2));
% remove glomeruli not tagged by gh146 
for i=1:length(glomNames)
    for j=1:length(gh146glomsR)
        if strfind(glomNames{i},gh146glomsR{j})
            glomNamesgh146{j}=glomNames{i};
            odorResponsesgh146(j,:)=odorResponses(i,:);
            break
        end
        %         if j==length(gh146glomsR)
        %             odorResponses(i,:)=NaN;
        %             glomNames{i}=[];
        %         end
    end
end


receptorNamesOrco=cell(1,length(orcogloms));
glomNamesOrco=cell(1,length(orcogloms));
odorResponsesOrco=NaN*zeros(length(orcogloms),size(odorResponses,2));
% remove glomeruli not tagged by gh146 
for i=1:length(glomNames)
    for j=1:length(orcoglomsR)
        if strfind(glomNames{i},orcoglomsR{j})
            receptorNamesOrco{j}=glomNames{i};
            glomNamesOrco{j}=orcogloms{j};
            odorResponsesOrco(j,:)=odorResponses(i,:);
            break
        end
        %         if j==length(gh146glomsR)
        %             odorResponses(i,:)=NaN;
        %             glomNames{i}=[];
        %         end
    end
end

publishedOR.gh146response=odorResponsesgh146;
publishedOR.gh146receptorNames=glomNamesgh146;
publishedOR.orcoresponse=odorResponsesOrco;
publishedOR.orcoreceptorNames=receptorNamesOrco;
publishedOR.orcoGlomerulusNames=glomNamesOrco;
publishedOR.gh146glomerulusNames=glomLocationNamesgh146;
publishedOR.gh146xCentroid=xcgh146;
publishedOR.gh146yCentroid=ycgh146;
publishedOR.gh146zCentroid=zcgh146;
publishedOR.gh146glomBorderX=xAllDatagh146;
publishedOR.gh146glomBorderY=yAllDatagh146;
publishedOR.gh146glomBorderZ=zAllDatagh146;
save(savename,'publishedOR')