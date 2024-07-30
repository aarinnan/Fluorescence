function wid=fluwid(x,em,ex);
% wid=fluwid(x,em,ex) - utility function for flumod
%
% Estimates the width of the 1st order Rayleigh peak. Can also be used as
% an estimate for 2nd order Rayleigh peak. (However, it is only calculated
% on the 1st order Rayleigh peak.)
%
% INPUT:
%  x     EEM (Dimension: Samples x Emission x Excitation)
%  em    Emission wavelengths
%  ex    Excitation wavelengths
%
% OUTPUT:
%  wid   Estimated width (can be used in either fluweight or flufac)

% Copyright, 2005 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Åsmund Rinnan
% E-mail asmundrinnan@gmail.com

dim=size(x);
if length(em)~=dim(2) | length(ex)~=dim(3)
    id=[find(dim~=length(em) & dim~=length(ex)) find(dim==length(em)) find(dim==length(ex))];
    disp('The EEM should be in the format Samples x Emission x Excitation')
    disp(['Type ''X=reshape(nshape(X,[' num2str(id) ']),[' num2str(dim(id)) ']);'''])
    disp('PS! Remember to change ''X'' to the variable name you are using!')
    disp(' ')
    error('Wrong dimension of X')
end

%Finds the part of the matrix which is below the 1st order Rayleigh line
X=x;
if sum(isnan(vec(x)))>0
    X(isnan(X))=max(vec(X));
end
x=squeeze(sum(X));
for i=1:length(ex)
    j=find(em>ex(i));
    x(j,i)=NaN;
end
%Finds the average value in this part of the data matrix
ref=nanmean(x(:));

%Finds the maximum distance from the 1st order Rayleigh line which is
%recorded. Limits this to 50nm
mx=floor((ex(end)-em(1))/5)*5;
if mx>50
    mx=50;
end

%Stepwise inserting missing values into the data matrix starting closest
%to the 1st order Rayleigh line
for i=5:5:mx
    for j=1:length(ex)
        k=find(em>ex(j)-i);
        x(k,j)=NaN;
    end
    %Calculating the average value on the data matrix not set to missing
    ray(i/5)=nanmean(x(:));
end
ray=[ref ray];

%Estimate the 1st order Rayleigh peak as the place where the average of the
%data matrix decreases with less than 15%
ch=(ray(1:end-1)-ray(2:end))./ray(1:end-1);
wid=0:5:mx;
if isempty(find(ch<.15))
    [i,j]=min(ch);
    wid=wid(j);
else
    wid=wid(1,min(find(abs(ch)<.15)));
end

%-------------------
function xnew=vec(x)
%xnew=vec(x)
%
% Vectorices the matrix x

xnew=x(:);