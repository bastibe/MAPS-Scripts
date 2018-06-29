function r=yin2(x,p);
% YIN2: a simple implementation of the yin period-estimation algorithm
%
%  yin2(x) : plot the period, power, and aperiodicity as a function of time
%
%  r=yin2(x,p): use parameters in p, return result in r:
%
%    r.prd: period 
%    r.ap: aperiodicity measure
%  
%    p.maxprd: samples, maximum of search range [default 100]
%    p.minprd: samples, minimum of search range [default 2]
%    p.wsize: samples, window size [default maxprd]
%    p.hop: samples, frame period [default wsize]
%    p.thresh: threshold for period minimum [default: 0.1]
%    p.smooth: samples, size of low-pass smoothing window [default: minprd/2]

if nargin>2; error('!'); end
if nargin<2; p=[]; end
if nargin<1; error('!'); end

% defaults
if ~isfield(p, 'maxprd'); p.maxprd=100; end
if ~isfield(p, 'minprd'); p.minprd=2; end
if ~isfield(p, 'wsize'); p.wsize=p.maxprd; end
if ~isfield(p, 'hop'); p.hop=p.wsize; end
if ~isfield(p, 'thresh'); p.thresh=0.1; end
if ~isfield(p, 'smooth'); p.smooth=p.minprd/2; end

if min(size(x)) ~= 1; error('data should be 1D'); end
x=x(:);
nsamples=numel(x);

nframes=floor((nsamples-p.maxprd-p.wsize)/p.hop);
pwr=zeros(1,nframes);
prd=zeros(1,nframes);
ap=zeros(1,nframes);

% shifted data
x=convmtx(x,p.maxprd+1);
x=x(p.maxprd:end-p.maxprd,:);

for k=1:nframes
   
    start=(k-1)*p.hop; % offset of frame
    xx=x(start+1:start+p.wsize,:);
    d=mean( (xx - repmat(xx(:,1),1,p.maxprd+1)).^2 )/2;     % squared difference function
    dd= d(2:end) ./ (cumsum(d(2:end)) ./ (1:(p.maxprd)));   % cumulative mean - normalized
    
    % parabolic interpolation of all triplets to refine local minima
    min_pos=1:numel(dd);    % nominal position of each sample
    x1=dd(1:end-2);
    x2=dd(2:end-1);
    x3=dd(3:end);
    a=(x1+x3-2*x2)/2;
    b=(x3-x1)/2;
    shift=-b./(2*a);        % offset of interpolated minimum re current sample
    val=x2-b.^2./(4*a);     % value of interpolated minimum
    
    % replace all local minima by their interpolated value, 
    idx= 1 + find(x2<x1 & x2<x3);
    dd(idx)=val(idx-1);
    min_pos(idx)=min_pos(idx-1)+shift(idx-1);
    
    % find index of first min below threshold
    a=dd<p.thresh;
    if isempty(find(a))
        [~,prd0]=min(dd); % none below threshold, take global min instead
    else
        b=min(find(a)); % left edge
        c=min(b*2,numel(a));
        [~,prd0]=min(dd(b:(c-1))); 
        prd0=b+prd0-1;
    end
    
    prd=min_pos(prd0)+1;
        
    if prd>2 & prd<numel(dd) & d(prd0)<d(prd0-1) & d(prd0)<d(prd0+1)
        
        % refine by parabolic interpolation of raw difference function 
        x1=d(prd-1);
        x2=d(prd);
        x3=d(prd+1);
        a=(x1+x3-2*x2)/2;
        b=(x3-x1)/2;
        shift=-b./(2*a);        % offset of interpolated minimum re current sample
        val=x2-b.^2./(4*a);     % value of interpolated minimum
        prd=prd+shift-1;
    end

    % aperiodicity
    frac=prd-floor(prd);
    if frac==0
        yy=xx(:,prd);
    else
        yy=(1-frac)*xx(:,floor(prd+1))+frac*xx(:,floor(prd+1)+1); % linear interpolation
    end
    pwr=(mean(xx(:,1).^2) + mean(yy.^2))/2; % average power over fixed and shifted windows
    res=mean(((xx(:,1) - yy)).^2) / 2;
    ap=res/pwr;
    
        
    r.prd(k)=prd;
    r.ap(k)=ap;

end

if nargout==0; 
    subplot 211; plot(r.prd); title('period'); xlabel('frame'); ylabel('samples');
    subplot 212; plot(r.ap); title('periodicity'); xlabel('frame');
    r=[]; 
end

