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

% load extracted glomeruli centroids
currFiles=dir(pwd);

dummyi=1;
for j=1:length(currFiles)
    if strfind(currFiles(j).name,'csv')
        filesToLoad{dummyi}=currFiles(j).name;
        dummyi=dummyi+1;
    end
end

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
        odorResponses(str2num(str2num(currdata{1}{j})),i)=str2num(currdata{1}{j+3});
        glomNames{str2num(str2num(currdata{1}{j}))}=str2num(currdata{1}{j+1});
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

publishedOR.gh146response=odorResponsesgh146;
publishedOR.gh146receptorNames=glomNamesgh146;
publishedOR.gh146glomerulusNames=glomLocationNamesgh146;
publishedOR.gh146xCentroid=xcgh146;
publishedOR.gh146yCentroid=ycgh146;
publishedOR.gh146zCentroid=zcgh146;
publishedOR.gh146glomBorderX=xAllDatagh146;
publishedOR.gh146glomBorderY=yAllDatagh146;
publishedOR.gh146glomBorderZ=zAllDatagh146;
save(savename,'publishedOR')