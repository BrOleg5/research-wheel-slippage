function colvec = anyvec2colvec(anyvec)
% ANYVEC2COLVEC Convert any vector to column vector
%   colvec = ANYVEC2COLVEC(anyvec) Convert any vector to column vector.
%   If anyvec is a column vector, then ANYVEC2COLVEC leaves it unchanged.
%
% See also: ANYVEC2ROWVEC

arguments
    anyvec {mustBeVector}
end

if(iscolumn(anyvec))
    colvec = anyvec;
elseif(isrow(anyvec))
    colvec = anyvec.';
end
end