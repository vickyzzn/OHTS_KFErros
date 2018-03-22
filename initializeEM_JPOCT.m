function [ A0 C0 Q0 R0 INITX0 INITV0 ] = initializeEM_JPOCT( data, type )
%INITIALIZEEM Summary of this function goes here
%   Detailed explanation goes here
%   type: 0 corresponds to MD, IOP, PSD
%         1 corresponds to MD, IOP, PSD, OCT
%         2 corresponds to MD, IOP, PSD, OCT, DH
numPat=size(data,1)-1;


if type == 0
    dcol = 3; %data column
    numvars = 9; %number of Kalamn Filter variables
elseif type == 1
    dcol = 13;
    numvars = 10;
else
    dcol = 14;
    numvars = 11;
end



INITX0=zeros(numvars,1);
for i =1:numPat
    INITX0=INITX0+data{i+1,dcol}(:,1); %include OCT
end
INITX0=INITX0/numPat;
[TM, ERR ] = kalman_est_cov_JPOCT( data(2:end,:), type );
A0 = TM;
C0 = eye(numvars,numvars);
INITV0 = ERR;
Q0 = ERR./2;
R0 = ERR./2;

Q0=triu(tril(Q0));

R0=triu(tril(R0));


end

