%
% Mesh generator for OFE_COUPLING 2
%
clear all;
% input
%
ndiv=100;
h=2.0/ndiv;
%
nump=(ndiv+1)*(ndiv+1);
nquad=ndiv*ndiv;
%
POINTS=zeros(nump,3);
POINTS_ONOFF1=false(nump,1);
POINTS_ONOFF2=false(nump,1);
QUADS=zeros(nquad,4);
QUADS_ONOFF1=false(nquad,1);
QUADS_ONOFF2=false(nquad,1);
%
n=0;
for i=1:ndiv+1
    for j=1:ndiv+1
        n=n+1;
        POINTS(n,1)=n;
        POINTS(n,2)=(i-1)*h;
        POINTS(n,3)=(j-1)*h;
        if(i==1)
           POINTS_ONOFF1(n)=true;
        end
    end
end
%
n=0;
for i=1:ndiv
    for j=1:ndiv
        n1=(i-1)*(ndiv+1)+j;
        n2=(i)*(ndiv+1)+j;
        n3=n2+1;
        n4=n1+1;
%
        n=n+1;
        QUADS(n,1)=n1;
        QUADS(n,2)=n2;
        QUADS(n,3)=n3;
        QUADS(n,4)=n4;
%
        if (i==1)|(i==ndiv)|(j==1)|(j==ndiv)
            QUADS_ONOFF1(n)=true;
        end
    end
end
n=0;
for i=1:ndiv
    for j=1:ndiv
        n=n+1;
        if (QUADS_ONOFF1(n))
            n1=QUADS(n,1);
            n2=QUADS(n,2);
            n3=QUADS(n,3);
            n4=QUADS(n,4);
%
            POINTS_ONOFF2(n1)=true;
            POINTS_ONOFF2(n2)=true;
            POINTS_ONOFF2(n3)=true;
            POINTS_ONOFF2(n4)=true;
        end
    end
end
[QUADS_ONOFF1,QUADS_ONOFF2]=ExtLogic(QUADS_ONOFF1,ndiv,ndiv,1);
%---------------------------------------------------------------
% Generate input data
%---------------------------------------------------------------
% NODE
NODE=zeros(nump,11);
n=0;
for i=1:ndiv+1
    for j=1:ndiv+1
        n=n+1;
%
        NODE(n,1)=n;
        NODE(n,4)=1;
        NODE(n,5)=POINTS(n,2);
        NODE(n,6)=POINTS(n,3);
%
        if(~POINTS_ONOFF2(n))
            NODE(n,8)=1;
            NODE(n,9)=1;
        else
            if(POINTS_ONOFF1(n))
               NODE(n,8)=3;
               NODE(n,9)=3;
            else
               NODE(n,8)=6;
               NODE(n,9)=6;
            end
        end
%             
    end
end
% ELE
n=0;
for i=1:ndiv
    for j=1:ndiv
        n=n+1;
        if (QUADS_ONOFF1(n))
            n1=QUADS(n,1);
            n2=QUADS(n,2);
            n3=QUADS(n,3);
            n4=QUADS(n,4);
%
            POINTS_ONOFF2(n1)=true;
            POINTS_ONOFF2(n2)=true;
            POINTS_ONOFF2(n3)=true;
            POINTS_ONOFF2(n4)=true;
        end
    end
end
points=POINTS(POINTS_ONOFF2,:);
numele=size(points,1);
nd_ele=zeros(nump,1);
ELE=zeros(numele,7);
for i=1:numele
    n=points(i,1);
    ELE(i,1)=i;
    ELE(i,2)=n;
    nd_ele(n)=i;
    if(POINTS_ONOFF1(n))
        ELE(i,3)=1;
        ELE(i,4)=1;     
    end
%
    ELE(i,5)=h;
    ELE(i,6)=h;
end
% Int-Tri
quads=QUADS(QUADS_ONOFF1&~QUADS_ONOFF2,:);
numquads=size(quads,1);
INTTRI=zeros(2*numquads,5);
n=0;
for i=1:numquads
    n1=quads(i,1);
    n2=quads(i,2);
    n3=quads(i,3);
    n4=quads(i,4);
%
    x1=POINTS(n1,2);
    x2=POINTS(n2,2);
    x3=POINTS(n3,2);
    x4=POINTS(n4,2);
%
    y1=POINTS(n1,3);
    y2=POINTS(n2,3);
    y3=POINTS(n3,3);
    y4=POINTS(n4,3);  
%
%    patch([x1 x2 x3],[y1 y2 y3],'y');
%    patch([x1 x3 x4],[y1 y3 y4],'y');    
%
    ele1=nd_ele(n1);
    ele2=nd_ele(n2);
    ele3=nd_ele(n3);
    ele4=nd_ele(n4);
%
    n=n+1;
    INTTRI(n,1)=n;
    INTTRI(n,2)=ele1;
    INTTRI(n,3)=ele2;
    INTTRI(n,4)=ele3;
    
%
    n=n+1;
    INTTRI(n,1)=n;
    INTTRI(n,2)=ele1;
    INTTRI(n,3)=ele3;
    INTTRI(n,4)=ele4;    
end
nn=n;
% Int-Rectangle
quads=QUADS(QUADS_ONOFF2,:);
numquads=size(quads,1);
INTRECT=zeros(numquads,6);
FEM2=zeros(numquads,6);
INTTRI2=zeros(2*numquads,5);
n=0;
for i=1:numquads
    n1=quads(i,1);
    n2=quads(i,2);
    n3=quads(i,3);
    n4=quads(i,4);
%
    x1=POINTS(n1,2);
    x2=POINTS(n2,2);
    x3=POINTS(n3,2);
    x4=POINTS(n4,2);
%
    y1=POINTS(n1,3);
    y2=POINTS(n2,3);
    y3=POINTS(n3,3);
    y4=POINTS(n4,3);
%
    patch([x1 x2 x3 x4],[y1 y2 y3 y4],'g');
%
    ele1=nd_ele(n1);
    ele2=nd_ele(n2);
    ele3=nd_ele(n3);
    ele4=nd_ele(n4);
%
    n=n+1;
    INTRECT(n,1)=n;
    INTRECT(n,2)=ele1;
    INTRECT(n,3)=ele2;
    INTRECT(n,4)=ele3;  
    INTRECT(n,5)=ele4;
%   
    FEM2(n,1)=n;
    FEM2(n,2)=n1;
    FEM2(n,3)=n2;
    FEM2(n,4)=n3;
    FEM2(n,5)=n4;
    FEM2(n,6)=1;
%  
    nn=nn+1;
    INTTRI2(nn,1)=nn;
    INTTRI2(nn,2)=ele1;
    INTTRI2(nn,3)=ele2;
    INTTRI2(nn,4)=ele3;    
%  
    nn=nn+1;
    INTTRI2(nn,1)=nn;
    INTTRI2(nn,2)=ele1;
    INTTRI2(nn,3)=ele3;
    INTTRI2(nn,4)=ele4;     
end
% BOUND
numbound=ndiv*2;
BOUND=zeros(numbound,10);
n=0;
for i=1:ndiv
    n=n+1;
    n1=i+1;
    n2=i;
%
    ele1=nd_ele(n1);
    ele2=nd_ele(n2);
%
    BOUND(n,1)=n;
    BOUND(n,2)=ele1;
    BOUND(n,3)=ele2;
    BOUND(n,4)=1;
    BOUND(n,5)=1;
    BOUND(n,8)=1;
    BOUND(n,9)=1;
%   
end
%
for i=1:ndiv
    n=n+1;
    n1=i*(ndiv+1)+ndiv+1;
    n2=(i-1)*(ndiv+1)+ndiv+1;
%
    ele1=nd_ele(n1);
    ele2=nd_ele(n2);
%
    BOUND(n,1)=n;
    BOUND(n,2)=ele1;
    BOUND(n,3)=ele2;
    BOUND(n,4)=2;
    BOUND(n,5)=2;
    BOUND(n,7)=1;
    BOUND(n,8)=1;
    BOUND(n,9)=1;
%   
end
% Generate quads
quads=QUADS(~QUADS_ONOFF1,:);
numfem=size(quads,1);
FEM=zeros(numfem,6);
for i=1:numfem
    n1=quads(i,1);
    n2=quads(i,2);
    n3=quads(i,3);
    n4=quads(i,4);
%
    x1=POINTS(n1,2);
    x2=POINTS(n2,2);
    x3=POINTS(n3,2);
    x4=POINTS(n4,2);
%
    y1=POINTS(n1,3);
    y2=POINTS(n2,3);
    y3=POINTS(n3,3);
    y4=POINTS(n4,3);   
%
   patch([x1 x2 x3 x4],[y1 y2 y3 y4],'b');
%
    FEM(i,1)=i;
    FEM(i,2)=n1;
    FEM(i,3)=n2;
    FEM(i,4)=n3;
    FEM(i,5)=n4;
    FEM(i,6)=1;
end

