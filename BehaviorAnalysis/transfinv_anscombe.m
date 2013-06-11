function origData = transfinv_anscombe(inputTransformedData)
%TRANSFINV_ANSCOMBE (sacvector2): invert the var-stabilizing Anscombe transf
%   origData = TRANSFINV_ANSCOMBE(inputTransformedData)
%
%   Can be used by PREF_DIR_BOTH_SEQUENTIAL, DECODE_RATE_TO_VECTOR
%
%  MH - http://github.com/histed/tools-mh

origData = sign(inputTransformedData) ...
           .* ((inputTransformedData.^2) - (3/8));

