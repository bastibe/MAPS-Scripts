function r=yin2(p,fileinfo) 
% YIN2 - fundamental frequency estimator
% new version (feb 2003)
%
%


% process signal a chunk at a time
idx=0;
totalhops=round(fileinfo.nsamples / p.hop);
r1=nan*zeros(1,totalhops);r2=nan*zeros(1,totalhops);
r3=nan*zeros(1,totalhops);r4=nan*zeros(1,totalhops);
idx2=0+round(p.wsize/2/p.hop);
while (1)
	xx=sf_wave(fileinfo, [idx+1, idx+p.bufsize], []);
	xx=xx(:,1); 			% first channel if multichannel
	[prd,ap0,ap,pwr]=yin_helper(xx,p);
	n=size(prd ,2);
	if (~n) break; end;
	idx=idx+n*p.hop;

	r1(idx2+1:idx2+n)= prd;	
	r2(idx2+1:idx2+n)= ap0;
	r3(idx2+1:idx2+n)= ap;
	r4(idx2+1:idx2+n)= pwr;
	idx2=idx2+n;	
end
size(r1)
r.r1=r1; % period estimate
r.r2=r2; % gross aperiodicity measure
r.r3=r3; % fine aperiodicity measure
r.r4=r4; % power
sf_cleanup(fileinfo);
% end of program



% Estimate F0 of a chunk of signal
function [prd,ap0,ap,pwr]=yin_helper(x,p,dd)
[m,n]=size(x);
smooth=ceil(p.sr/p.lpf);
x=rsmooth(x,smooth); 	% light low-pass smoothing
x=x(smooth+1:end);
maxlag = ceil(p.sr/p.minf0); 	
minlag = floor(p.sr/p.maxf0);

hops=floor((m-maxlag-p.wsize)/p.hop);
prd=zeros(1,hops);
ap0=zeros(1,hops);
ap=zeros(1,hops);
pwr=zeros(1,hops);
if hops<1; return; end

% difference function matrix
dd=zeros(ceil((m-maxlag-2)/p.hop),maxlag+2); % +2 to improve interp near maxlag
lags=[zeros(1,maxlag+2); 1:maxlag+2]; 
rdiff_inplace(x,x,dd,lags,p.hop);
rsum_inplace(dd,round(p.wsize/p.hop));
dd=dd';

% parabolic interpolation near min, then cumulative mean normalization
[dd,ddx]=minparabolic(dd);
cumnorm_inplace(dd);

for j=0:hops-1
	idx=j*p.hop;
	d=dd(:,j+1);
	dx=ddx(:,j+1);
	
	% estimate period
	pd=dftoperiod(d,[minlag,maxlag],p.thresh*2); 

	% gross aperiodicity is value of cumnormed df at minimum:
	ap0(j+1)=d(pd)/2;
	
	% fine tune period based on parabolic interpolation
	pd=pd+dx(pd+1)+1;

	% power estimates
	k=(1:ceil(pd))';
	x1=x(k+idx);
	x2=k+idx+pd-1;
	interp_inplace(x,x2);
	x3=x2-x1;
	x4=x2+x1;
	x1=x1.^2; rsum_inplace(x1,pd);
	x3=x3.^2; rsum_inplace(x3,pd);
	x4=x4.^2; rsum_inplace(x4,pd);
	
	x1=x1(1)/pd;
	x2=x2(1)/pd;
	x3=x3(1)/pd;
	x4=x4(1)/pd;
	% total power
	pwr(j+1)=x1;

	% fine aperiodicity
	ap(j+1)=(eps+x3)/(eps+(x3+x4)); % accurate, only for valid min

	prd(j+1)=pd;
end


