close all
% clear all

loadModules();

ACRaw=readRawData('rawAGISCIGTS.xlsx');
[ACInterp, raw]=interpolateData(ACRaw);
ACProg=labelProgression(ACInterp);

[ OHTSProgTrain,OHTSProgTest,OHTSRawTrain,OHTSRawTest ] = splitTrainTest( ACProg,ACRaw, .99);

% %Fix dimensions of files so col 8 would be the same as other training and
% %testing files
% OHTSProgTrain(:,[8])=[];
% OHTSProgTest(:,[8])=[];

table1=Table1(OHTSRawTrain,OHTSRawTest,OHTSProgTrain,OHTSProgTest)

[ A0 C0 Q0 R0 INITX0 INITV0 ]=initializeEM(OHTSProgTrain);
[A_OHTS, C_OHTS, Q_OHTS, R_OHTS, INITX_OHTS, INITV_OHTS, LL_OHTS] = learn_kalman(OHTSProgTrain(2:end,3), A0, C0, Q0, R0, INITX0, INITV0,100);

getRegModel( A_OHTS, C_OHTS, Q_OHTS, R_OHTS, INITX_OHTS, INITV_OHTS, OHTSProgTrain);
o=readRegCoeff();

% [acc,dd]= paretoAnalysis( A_AC,C_AC,Q_AC,R_AC,INITX_AC,INITV_AC,o, ACProgTrain,ACProgTest,0.1,'AC Training','AC Testing');
% 
% 
% beta = 0.98;
% rho=0.37;

save 'OHTSdata.mat'
%  [~,fd]=forwardTNT(A, C, Q, R, INITV, INITX,o,beta,rho,ACProgTest,6);
% plotTrjectory([],fd,ACRawTest,'');

