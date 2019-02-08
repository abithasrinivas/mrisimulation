%MRI Project, 2017
%Author : Sai Abitha Srinivas
%all times in ms, all freqs in kHz, all distances in cm, all B's in T

clc;
clear all;
close all; 

Fovx=20;
Fovy=10;
Nx=64;
Ny=32;

complexobj = 1; % 1 if using complex data
 if complexobj
    load object17;
 else
    xobj = 2;		% cm
    yobj = 2;		% cm 
    zobj = 0;		% cm
    T1obj = 1000;	% ms
    T2obj = 100;    % ms
 end
 
nobj = length(xobj);
gambar = 42570;               % gamma/2pi in kHz/T
gam = gambar*2*pi;            % gamma in kiloradians/T

% simulation values
dt = 32e-3;                   % ms
te = 10.0;                    % ms
endtime =14;                  % ms
time = [0:dt:endtime]';       % ms
npts = length(time);          % number of time points for simulation
bx = zeros([npts nobj]);      % place holders
by = zeros([npts nobj]);      % place holders
bz = zeros([npts nobj]);      % place holders

% define 90 RF
pwrf90 = 2.56;                % in ms
sincper = pwrf90/4;           % in ms
rfsteps = pwrf90/dt;
rftime = [-(rfsteps-1)/2:(rfsteps-1)/2]'.*dt;
rfshape = hanning(rfsteps).*sinc(rftime./sincper);
a_rf90 = (pi/2)./(gam.*sincper);                  % in T 
b1_90 = a_rf90.*[rfshape; zeros([npts-rfsteps 1])]; % in T

kxpos=[];% place holders
kypos=[];% place holders
M=[];% place holders

% define gz
slthick =1;  % in cm
rfbw=1./(sincper); % kHz
a_gz1 = rfbw./(gambar.*slthick); % in T/cm
pwgz1 = pwrf90;
a_gz2 = -a_gz1;
pwgz2 = pwrf90/2;
gz =  (time < pwgz1) .* a_gz1 ...
          + (time >= pwgz1).*(time < (pwgz1+pwgz2)) .* a_gz2;

%define gx 
pwgx2=4.096;
deltat=pwgx2./Nx;
RcvBw=1./(deltat);
a_gx2=(2*pi*RcvBw)./(gam.*20);
a_gx1=-2*a_gx2;
pwgx1=1.024;
gx=(time>(te-(pwgx2/2)-pwgx1)).*(time<te-(pwgx2./2)).*a_gx1...
        +(time >= te-(pwgx2/2)).*(time<=te+(pwgx2/2)).*a_gx2;    

dely=Fovy./Ny;
kymax=1./(2.* dely);
delky=(2.* kymax)./Ny;
pwgy=1.024;
         
% set the initial magnitization
m0 = [0 0 1]'*ones([1 nobj]);
         
for slice = 1:2  % slice loop
    npe=32;
     for pe = 1:npe % phase encode loop
        disp(sprintf('PE %d of %d, Slice %d', pe, npe,slice))       
        Gymaxarea=((kymax)-pe.*delky)./(gambar);
        a_gy=-(Gymaxarea./pwgy);
        gy=(time>(te-(pwgx2/2)-pwgx1)).*(time<te-(pwgx2./2)).*a_gy;
        bx =repmat(b1_90,[1 nobj]);
        by =zeros(size(bx));
        bz= gz .*zobj + gx .*xobj + gy .* yobj;

        [mx my mz] = blochsim2(m0,bx,by,bz,T1obj,T2obj,dt);

        gz=repmat(gz,[1 nobj]);
        gy=repmat(gy,[1 nobj]);
        gx=repmat(gx,[1 nobj]);

        for m=1:438
        bz(m,:)=gz(m,:).*zobj+gy(m,:).*yobj+gx(m,:).*xobj;
        end

        gx=gx(:,1);
        gy=gy(:,1);
        gz=gz(:,1);

        % set the initial magnitization
        m0 = [0 0 1]'*ones([1 nobj]);
        
        % do Bloch equation simulation 
        [mx my mz] = blochsim2(m0,bx,by,bz,T1obj,T2obj,dt);

        nread=64;
        temp=(time >= te-(pwgx2/2)).*(time<=te+(pwgx2/2));
        ind=find(temp==1);
        ind=ind(1:2:end);
        Gxt=gx(gx ~= 0);
        
        if isempty(Gxt)
            Gxt=zeros(Nx,1);
        end

        Gyt=gy(gy~=0);
        if isempty(Gyt)
            Gyt=zeros(Ny,1);
        end
        kxt(1)=gambar.*Gxt(1).*dt;

        for itr=2:length(Gxt)
            kxt(itr)=kxt(itr-1)+gambar.*Gxt(itr).*dt;
        end

        kyt(1)=gambar.*Gyt(1).*dt;

        for itr=2:length(Gyt)
            kyt(itr)=kyt(itr-1)+gambar.*Gyt(itr).*dt;
        end

        kxread = kxt(160-127:end);
        kxread=kxread(1:1:end);
        kyread=kyt(end).*ones(size(kxread));
        kxpos= kxread;
        kypos=[kypos kyread];
        sig=(mx(ind,:)+1i.*my(ind,:));
        M(pe,:,:)=sig;
        
    end  % phase encode loop
end 

subplot(3,1,1)
plot(time,mx);
xlabel('time (ms)');
ylabel('Mx');
%             axis([0 15 -1 1]);
title('magnetization v/s time (T1 = 1000ms , T2 = 100ms), with gz and gx');

subplot(3,1,2)
plot(time,my);
xlabel('time (ms)');
ylabel('My');
%             axis([0 15 -1 1]);

subplot(3,1,3)
plot(time,mz);
xlabel('time (ms)');
ylabel('Mz');
%             axis([0 15 -1 1]);

NewM = sum (M, 3);
M = NewM;
xpos = [-nread/2:nread/2-1]/nread*20;
ypos = [-Ny/2:Ny/2-1]/Ny*10;
im=fftshift(ifft2(M));

figure; imagesc(kxpos,kypos,abs(M)); colormap gray; axis('image'); axis('xy')
xlabel('kx (cm-1)');
ylabel('ky (cm-1)');
title('abs(M(Kx,Ky))');

figure;imagesc(kxpos,kypos,imag(M)); colormap gray; axis('image'); axis('xy')
xlabel('kx (cm-1)');
ylabel('ky (cm-1)');
title('imag(M(Kx,Ky))')

figure;imagesc(xpos,ypos, abs(im)); colormap gray; axis('image'); axis('xy')
xlabel('x (cm)');
ylabel('y (cm)');
title('abs(image)')
