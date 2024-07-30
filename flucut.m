function [Xnew,EmNew,ExNew] = flucut(X,EmAx,ExAx,Cut1,Cut2)
%[Xnew,EmNew,ExNew] = flucut(X,EmAx,ExAx,Cut1,Cut2)
% 
% Inserts a mixture of NaN and 0 values outside the data area of 
% an EEM (Excitation-Emission Matrix).
%
% INPUT:
%    X      The X-array from the fluoressence spectra
%    EmAx   The Emitation Axis
%    ExAx   The Excitation Axis
%    Cut1+  Two numbers indicating how far away from the Rayleigh scatter
%     Cut2   line NaN and 0 values should be insterted. Negative values
%            indicate that the NaN's or 0's should be inserted above the
%            Rayleigh scatter line. First number is for NaN values the 
%            the second for 0's. NaN values indicate that nothing should 
%            be insterted. E.g. [NaN 10] means that no NaN-values should be
%            inserted, but 0's should be inserted 10nm from the Rayleigh
%            scatter line. Cut1 is for 1st order Rayleigh, Cut2 is for 2nd
%            order Rayleigh.
%
% OUTPUT:
%    Xnew   New X-array with inserted NaN and 0-values. Columns with only missing
%            values are removed
%    EmNew  Same reduced dimension as Xnew
%    ExNew  Same reduced dimension as Xnew
%
% See also: cutflu

%Based on: Lisbeth Garbrecht Thygesen, Åsmund Rinnan, Søren Barsberg, Jens
%           K.S. Møller (2004): Stabilizing the PARAFAC decomposition of
%           fluorescence spectra by insertion of zeros outside the data
%           area, Chemometrics and Intelligent Laboratory Systems, 71,
%           97-106. DOI: 10.1016/j.chemolab.2003.12.012

% OMAL 210504

%Checks if the vectors are correctly given.
if length(EmAx)~=size(X,2)
    if length(EmAx)==size(X,2)
        error('You have mixed Em and Ex')
    else
        error('The correct Em and Ex are not given')
    end
end

%Lines 39-50 inserts values into Cut1 and Cut2 if these are not given.
if nargin<4
   Cut1=ones(1,2)*NaN;
end
if length(Cut1)<2
    Cut1=[Cut1 NaN];
end
if nargin<5
   Cut2=ones(1,2)*NaN;
end
if length(Cut2)<2
    Cut2=[Cut2 NaN];
end

%Runs an insertion twice. First with NaN's and then with 0's.
in=[NaN 0];
for i=1:2
    %If the Cut parameter is set to NaN noting is done.
    if ~isnan( Cut1(i))
        for j = 1:length(ExAx)
            k = find(EmAx<(ExAx(j)-Cut1(i)));
            X(:,k,j)=in(i);
        end
    end
    
    if ~isnan( Cut2(i))
        for j = 1:length(ExAx)
            k = find(EmAx>(ExAx(j)*2+Cut2(i)));
            X(:,k,j)=in(i);
        end
    end
end

%The matrix is reduced in size in the cases when NaN's have been introduced
%into whole slabs in the 2nd or 3rd mode of the original data. Nothing
%happens here in the cases where 0's have been inserted.
n=squeeze(isnan(X(1,:,:)));
i=find(sum(n)<size(X,2));
j=find(sum(n')<size(X,3));
Xnew=X(:,j,i);
EmNew=EmAx(j);
ExNew=ExAx(i);