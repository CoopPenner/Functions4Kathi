GCaseTable=readtable('/Users/pennerc/Downloads/NY_data.txt');
dx=GCaseTable.Disease;
age=GCaseTable.Age_at_sample;
ageAtOnset=GCaseTable.Age_at_onset_PD;
gcaseActivity=GCaseTable.GCase_activity;
riskSNP=GCaseTable.rs199347;
N370S=GCaseTable.N370S;
E326K=GCaseTable.E326K;
T369M=GCaseTable.T369M;


Parkinson= (contains(dx, 'PD'));
HC= (contains(dx, 'Control'));


overGPNMB=riskSNP==0; %the major allele
het=riskSNP==1;
underGPNMB= riskSNP==2; %the minor allele



%%
pt2Plot=[ Parkinson  ];
val2Plot=gcaseActivity;
SNP=riskSNP;
remVals= isnan(val2Plot) | isempty(SNP) | val2Plot<1  ;





figure
b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot )   )  )  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('Gcase Activity')
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & het & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot))+.5, val2Plot(~remVals & overGPNMB & pt2Plot), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(~remVals & het & pt2Plot))+2.5, val2Plot(~remVals & het & pt2Plot), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(~remVals & underGPNMB & pt2Plot))+4.5, val2Plot(~remVals & underGPNMB & pt2Plot), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GA','GG (under Production)'})



remVals= isnan(gcaseActivity) | isempty(SNP)| ~Parkinson;



ageVar=age(~remVals);
ageAtOnsetVar=ageAtOnset(~remVals);
gcaseActivityVar=gcaseActivity(~remVals);
snpVar=categorical(SNP(~remVals));
dxVar=categorical(dx(~remVals));

%% 

varNames=["Age","ageAtOnset", "gcase","rs199347",'dx'];

glmeTable = table(ageVar,ageAtOnsetVar,gcaseActivityVar,snpVar,dxVar, 'Variablenames',varNames);




glme = fitglme(glmeTable,...
		'gcase ~ 1  + rs199347 + Age + ageAtOnset ',...
		'Distribution','normal','Link','identity','FitMethod','Laplace',...
		'DummyVarCoding','reference')


glmeHeldOut=fitglme(glmeTable,...
		'cogSlope ~ 1 + Sex + startScore + (1|ageAtTest)   ',...
		'Distribution','normal','Link','identity','FitMethod','Laplace',...
		'DummyVarCoding','reference');
