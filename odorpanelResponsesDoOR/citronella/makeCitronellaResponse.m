% Matt Churgin, September 2018
clear all
close all


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

numOdors=2;
filesToLoad{1}='geraniol.csv';
filesToLoad{2}='limonene.csv';
% weight of each odor to citronella CEYLON TYPE (from wikipedia) https://en.wikipedia.org/wiki/Citronella_oil
% omitted other others because they have very little response according to
% door
weights=[.2 .1];

numGlom=78;
% initialize odor matrix
odorResponses=NaN*zeros(numGlom,numOdors);
glomNames=cell(1,numGlom);

for i=1:numOdors
    fid = fopen(filesToLoad{i});
    currdata = textscan(fid, '%s','Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
    fclose(fid);
    for j=5:4:length(currdata{1})
        odorResponses(str2num(str2num(currdata{1}{j})),i)=str2num(currdata{1}{j+3});
        glomNames{str2num(str2num(currdata{1}{j}))}=str2num(currdata{1}{j+1});
    end
end

citronellaResponses=nansum(odorResponses.*weights,2)/sum(weights);

% if a glom has no response to all odors, set to nan
inds=find((any(odorResponses,2)));
toremove=~any(odorResponses,2);
citronellaResponses(toremove)=[];
glomNames(toremove)=[];
glomNamesT=glomNames';

indsPlusResponse(:,1)=inds;
indsPlusResponse(:,2)=citronellaResponses;

fname='citronella_Output.csv';
fid = fopen(fname, 'w') ;
fprintf(fid, '%s\n', glomNamesT{1:end-1}) ;
fprintf(fid, '%s\n', glomNamesT{end}) ;
fclose(fid) ;
dlmwrite(fname, [indsPlusResponse], '-append','roffset',0,'coffset',1);
