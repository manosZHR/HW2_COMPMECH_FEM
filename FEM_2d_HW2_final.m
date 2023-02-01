function [xpt ypt u sk r1]=FEM_2d_HW2_final(nex,ney,theta)
% 2-dimensional linear problem. 
% Quadratic finite element basis functions.
%
% MATLAB version of the original FORTRAN code 1d_FEM_quad.f90
%
% Syntax
%  [x y u]=FEM_2d_quad(nex,ney)
%
% Parent function: FEM_2d_quad
% Nested functions: xydiscr, nodnumb, xycoord, abfind, tsfun
%
% This code uses nested functions.
% The primary difference between nested functions and other types
% of functions is that they can access and modify variables
% that are defined in their parent functions.
% http://www.mathworks.com/help/matlab/matlab_prog/nested-functions.html

%Global variables
nnx = 2*nex+1;
nny = 2*ney+1;
ne  = nex*ney;
np  = nnx*nny;
h=1;
To=20.;
sk = []; r1 = []; u = [];
xorigin = []; yorigin = []; thetaorigin=[]; xlast = []; xc=[]; ylast = []; thetalast=[]; 
deltax = []; deltay = []; deltatheta=[]; xpt = []; ypt = []; theta_i=[];
nop = []; ntop = []; nright = []; nleft = [];
ncod = []; bc = [];
w = []; gp =[];
phi = []; phic = []; phie = [];


fprintf('\t 2-D problem. Biquadratic basis functions \n')

xydiscr();
nodnumb();
xycoord();

fprintf('nex=%d, ney=%d, ne=%d, np=%d\n',nex,ney,ne,np)

% prepare for essential boundary conditions
ncod = zeros(np,1);
bc   = zeros(np,1);


ncod(1:nny:np-nny+1) = 1;
bc(1:nny:np-nny+1)   = 200.;

% prepare for natural boundary conditions
ntop = zeros(ne,1);


ntop(ney:ney:ne)  = 1;
nright(ne-ney+1:ne) = 1;
nleft(1:ney) = 1;

% initialization
r1 = zeros(np,1);
sk = zeros(np,np);

% matrix assembly
for nell=1:ne
    abfind(nell)
end

% impose essential boundary conditions
for i=1:np
    if(ncod(i)==1)
        r1(i)=bc(i);
        sk(i,1:np)=0.;
        sk(i,i)=1.;
    end
end

% solve the system of equations by Gauss elimination
u = sk\r1;

    function xydiscr()
        xorigin = 0.;
        yorigin = 3.;
        thetaorigin = 90.;
        
        xc = 6-3/tand(theta);
        xlast   = xc;
        ylast   = 0.;
        thetalast = theta;
        
        deltax = (xlast-xorigin)/nex;
        deltay = abs(ylast-yorigin)/ney;
        deltatheta = (thetaorigin-thetalast)/nex;

    end

    function nodnumb()
    % *** nodal numbering
        nel=0;
        for i=1:nex
            for j=1:ney
                nel=nel+1;
                for k=1:3
                    l=3*k-2;
                    nop(nel,l)=nny*(2*i+k-3)+2*j-1;
                    nop(nel,l+1)=nop(nel,l)+1;
                    nop(nel,l+2)=nop(nel,l)+2;
                end
            end
        end
    end

    function xycoord()
    % *** (x,y) coordinates of each node
        xpt=zeros(np,1);
        ypt=zeros(np,1);
        theta_i=zeros(np,1);
        
        xpt(1)=xorigin;
        ypt(nny)=yorigin;
        theta_i(1)=thetaorigin;
        
        for i=1:nnx
            nnode=(i-1)*nny+nny;
            xpt(nnode)=xpt(1)+(i-1)*deltax/2.;
            ypt(nnode)=ypt(nny);
            theta_i(nnode)=theta_i(1)-(i-1)*deltatheta/2;
            
            for j=2:nny
                theta_i(nnode+j-1)=theta_i(nnode); 
                
                if i~=1
                    dx1=deltay/(2*tand(theta_i(nnode+j-1)));

                    xpt(nnode-j+1) = xpt(nnode)+(j-1)*dx1;
                    ypt(nnode-j+1) = ypt(nnode)-(j-1)*deltay/2.;

                elseif i==1
                    xpt(nnode-j+1)=xpt(nnode);
                    ypt(nnode-j+1)=ypt(nnode)-(j-1)*deltay/2.;

                end               

            end
        end

        
    end


    
    function tsfun(x,y)
        l1  =@(c)  2.*c^2-3.*c+1.;
        l2  =@(c) -4.*c^2+4.*c;
        l3  =@(c)  2.*c^2-c;
        dl1 =@(c)  4.*c-3.;
        dl2 =@(c) -8.*c+4.;
        dl3 =@(c)  4.*c-1.;
        
        phi(1)=l1(x)*l1(y);
        phi(2)=l1(x)*l2(y);
        phi(3)=l1(x)*l3(y);
        phi(4)=l2(x)*l1(y);
        phi(5)=l2(x)*l2(y);
        phi(6)=l2(x)*l3(y);
        phi(7)=l3(x)*l1(y);
        phi(8)=l3(x)*l2(y);
        phi(9)=l3(x)*l3(y);
        phic(1)=l1(y)*dl1(x);
        phic(2)=l2(y)*dl1(x);
        phic(3)=l3(y)*dl1(x);
        phic(4)=l1(y)*dl2(x);
        phic(5)=l2(y)*dl2(x);
        phic(6)=l3(y)*dl2(x);
        phic(7)=l1(y)*dl3(x);
        phic(8)=l2(y)*dl3(x);
        phic(9)=l3(y)*dl3(x);
        phie(1)=l1(x)*dl1(y);
        phie(2)=l1(x)*dl2(y);
        phie(3)=l1(x)*dl3(y);
        phie(4)=l2(x)*dl1(y);
        phie(5)=l2(x)*dl2(y);
        phie(6)=l2(x)*dl3(y);
        phie(7)=l3(x)*dl1(y);
        phie(8)=l3(x)*dl2(y);
        phie(9)=l3(x)*dl3(y);
    end
    
    function abfind(nell)
        w  = [0.27777777777778, 0.444444444444, 0.27777777777778];
        gp = [0.1127016654    , 0.5           , 0.8872983346    ];
        
        ngl = nop(nell,1:9);
        
        for j=1:3 %LOOP j
            for k=1:3 %LOOPk
                tsfun(gp(j),gp(k))
                % *** isoparametric transformation
                x1=0;x2=0;y1=0;y2=0;
                for n=1:9
                    x1=x1+xpt(ngl(n))*phic(n);
                    x2=x2+xpt(ngl(n))*phie(n);
                    y1=y1+ypt(ngl(n))*phic(n);
                    y2=y2+ypt(ngl(n))*phie(n);
                end
                dett=x1*y2-x2*y1;
                for i=1:9
                    tphx(i)=(y2*phic(i)-y1*phie(i))/dett;
                    tphy(i)=(x1*phie(i)-x2*phic(i))/dett;
                end
                % *** residuals
                for l=1:9
                    for m=1:9
                        sk(ngl(l),ngl(m)) = sk(ngl(l),ngl(m)) ...
                                            -w(j)*w(k)*dett*(tphx(l)*tphx(m)+tphy(l)*tphy(m));
                    end
                end
            end %LOOP k
        end %LOOP j
        
        if(ntop(nell)==1)
            for k1=1:3
                tsfun(gp(k1),1.)
                x=0.; x1=0;
                for n=1:9
                    x=x+xpt(ngl(n))*phi(n);
                    x1=x1+xpt(ngl(n))*phic(n); %theta(x)/theta(ksi)
                end
                for k11=3:3:9
                    r1(ngl(k11))=r1(ngl(k11))+w(k1)*x1*phi(k11)*h*To;
                    for k12=3:3:9
                        sk(ngl(k11),ngl(k12))=sk(ngl(k11),ngl(k12))... 
                            - h*w(k1)*x1*phi(k11)*phi(k12);
                    end
                end
            end
        end
        

        
    end
end
