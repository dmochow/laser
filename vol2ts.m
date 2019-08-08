function bold2D=vol2ts(bold,mask)
% bold must be nx x ny x nz x nTR
% TIME DIMENSION MUST BE FOURTH
%[nx,ny,nz,nTR]=size(bold);
bold2D=permute(bold,[4 1 2 3]);  % move time to first dimension
bold2D=bold2D(:,:);
if nargin>1
    try
        bold2D=bold2D(:,mask);
    catch
        bold2D=bold2D(:,logical(mask));
    end
end