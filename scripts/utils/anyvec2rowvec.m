function rowvec = anyvec2rowvec(anyvec)
% ANYVEC2ROWVEC Convert any vector to row vector
%   rowvec = ANYVEC2ROWVEC(anyvec) Convert any vector to row vector.
%   If anyvec is a row vector, then ANYVEC2ROWVEC leaves it unchanged.
%
% See also: ANYVEC2COLUMNVEC

arguments
    anyvec {mustBeVector}
end

if(isrow(anyvec))
    rowvec = anyvec;
elseif(iscolumn(anyvec))
    rowvec = anyvec.';
end
end