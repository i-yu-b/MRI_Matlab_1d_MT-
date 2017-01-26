function sf = apodize2d(s,alpha)

if ~exist('alpha','var'), alpha = 0.5; end;
sf = zeros(size(s));
%nr=np,nc=nv,ne=ne,nd=ns
[nr,nc,ne,nd,nx] = size(s); %3D=[np/2,nv,nv2,ne]
%What is nx?

for kx = 1:nx %none
    for kd = 1:nd %none
        for ke = 1:ne
            w = zeros(nr,nc); %2D image size
            ss = s(:,:,ke,kd,kx); %raw signal (2D kspace)
            [mx,rax] = max(ss); %mx=max value, rax=index of max value in 1st dimension
            [mx2,cmx] = max(mx); %max value of maximums array for 2nd dim
            rmx = rax(cmx); %index in vector array 1st dimension for max of maxes
            lr = 2*(nr/2-abs(nr/2-rmx)); %2*(midpt of read-midpt-max index)=2*where to center window
            dr = (nr-lr)/2; %(read length-(2*peak))/2=how many units off center
            wr = tukeywin(lr,alpha); %define window length, portion of taper
            %2nd dimension
            lc = 2*(nc/2-abs(nc/2-cmx)); 
            dc = (nc-lc)/2;
            wc = tukeywin(lc,alpha);
            w(dr+1:end-dr,dc+1:end-dc) = wr*wc'; %defines window length of lr and lc
            sf(:,:,ke,kd,kx) = ss.*w; %multiply data by 2D window

        end
    end
end
