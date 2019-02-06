% from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;f15c35e.0707

function h = fast_fslgamma(t,meanlag,stddev)
% h = fast_fslgamma(t,<meanlag>,<stddev>)
%
% FSL's Gamma Function
%
% meanlag = 6 sec by default
% stddev  = 3 sec by default
% No normalization is done (yet)
%
% $Id: fast_fslgamma.m,v 1.1 2007/06/19 18:49:38 greve Exp $
%

% fast_fslgamma
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense%
% General inquiries: [log in to unmask]
% Bug reports: [log in to unmask]
%

h = [];
if(nargin < 1 | nargin > 3)
    fprintf('h = fast_fslgamma(t,meanlag,stddev)\n');
    return;
end
if(~exist('meanlag')) meanlag = 6; end
if(~exist('stddev'))  stddev = 3; end
a = (meanlag/stddev)^2;
b = meanlag/(stddev^2);
h = pdf_gamma(t,a,b);
return;

function pdfx = pdf_gamma(x,a,b)
% pdfx = pdf_gamma(x,a,b)
% Gamma distribution evaluated at x for parameters a and b
%
% pdfx = (b.^2) .* (x.^(a-1)) .* exp(-b.*x) ./ gamma(a);
% mean = a/b
% var  = a/(b^2)
% mode = (a-1)/b --> this the value of x at the peak
%
% When creating a hemodynamic response (eg, as in FSL), then
%   a = (tm/sigma)^2
%   b = tm/(sigma^2)
% where sigma = "Stddev" and tm = "Mean lag"
%
%
%
% pdf_gamma.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:'


pdfx = [];

if(nargin ~= 3)
    fprintf('pdfx = pdf_gamma(x,a,b)\n');
    return;
end

nx = length(x);
pdfx = zeros(nx,1);
ind = find(x>0);
pdfx(ind) = (b.^2) .* (x(ind).^(a-1)) .* exp(-b.*x(ind)) ./ gamma(a);

return;