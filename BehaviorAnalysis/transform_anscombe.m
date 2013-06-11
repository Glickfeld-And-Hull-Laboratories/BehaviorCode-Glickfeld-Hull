function transformedData = transform_anscombe(inputData)
%TRANSFORM_ANSCOMBE (sacvector2): do a var-stabilizing Anscombe transformation
%   transformedData = TRANSFORM_ANSCOMBE(inputData)
%
%   Can be used by PREF_DIR_BOTH_SEQUENTIAL, DECODE_RATE_TO_VECTOR
%
%  MH - http://github.com/histed/tools-mh

transformedData = sqrt(inputData + (3/8));
