function [expand,extend]=ExtLogic(onOff,dimI,dimJ,step)
%
%   Expand the logical '1' in the matrix OnOff
%   
%   Lingbo Zhang: 2017-2018
%   Email          : lingboz2015@gmail.com
%   Last updated   : 22/04/2017 with MATLAB R2014a
%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
   error('Insufficient inputs');
end
if(~islogical(onOff))
   error('Input is not a logical array');
end 
%% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   onOff=reshape(onOff,[dimI,dimJ]);
   expand=onOff;
%% Expand the logical array
   for i=1:step
       mLeft=[expand(:,2:dimJ),expand(:,1)];
       mRight=[expand(:,dimJ),expand(:,1:dimJ-1)];
       expand=expand|mLeft|mRight;
       mUp=[expand(2:dimI,:);expand(1,:)];
       mDown=[expand(dimI,:);expand(1:dimI-1,:)];
       expand=expand|mUp|mDown;
   end
%% Obtain the extend logical array
   extend=~onOff&expand;
   extend=reshape(extend,[dimI*dimJ,1]);
   expand=reshape(expand,[dimI*dimJ,1]);
end