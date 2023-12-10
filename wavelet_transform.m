

%   Usage:  c = fwt(f,w,J);
%           c = fwt(f,w,J,dim);
%           [c,info] = fwt(...);
% 
%   Input parameters:
%         f     : Input data.
%         w     : Wavelet definition.
%         J     : Number of filterbank iterations.
%         dim   : Dimension to along which to apply the transform.
% 
%   Output parameters:
%         c      : Coefficient vector.
%         info   : Transform parameters struct.


% Wavelet w options:
% Several formats of the basic filterbank definition w are recognized.
% One of them is a text string formed by a concatenation of a function 
% name with the wfilt_ prefix followed by a list of numerical arguments
% delimited by :. For example 'db10' will result in a call to 
% wfilt_db(10) or 'spline4:4' in call to wfilt_spline(4,4) etc.
% All filter defining functions can be listed by running
% dir([ltfatbasepath,filesep,'wavelets',filesep,'wfilt_*']);
% Please see help of the respective functions and follow references
% therein.
% 
% For other recognized formats of w please see FWTINIT.

pkg load signal
pkg load ltfat

[f,fs] = greasy;
J = 10;
[c,info] = fwt(f,'db8',J);
plotwavelets(c,info,fs,'dynrange',90);


wdef='db10';
fwtinit(wdef);

fwt()

% FWTINIT  Wavelet Filterbank Structure Initialization
%   Usage:  w = fwtinit(wdef);
%           w = fwtinit(wdef,prefix);
%           [w,info]=fwtinit(...)
% 
%   Input parameters:
%         wdef   : Wavelet filters specification.
%         prefix : Function name prefix
% 
%   Output parameters:
%         w    : Structure defining the filterbank.