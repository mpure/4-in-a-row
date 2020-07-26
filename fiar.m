function fiar(varargin)
m = 7;
n = 9;
X = zeros(m,n);
t = 1;
drawbd(X);
while 1
    p = varargin{t};
    flag = -1;
    while flag == -1
        if t==2
            X = 3-X;
            X(X==3) = 0;
        end
        j = feval(p,X);
        [i,flag] = loc(X,j);
    end
    if t==2
        X = 3-X;
        X(X==3) = 0;
    end
    X(i,j) = t;
    drawbd(X);
    if checkbd(X,i,j)
        title([p ' win!'],'fontsize',18);
        return;
    end
    t = 3 - t;
end

function j = person(~)
[x,y] = ginput(1);
j = ceil(x);

function j=ai4iar(X)
n1=0;
n2=0;
Z=X;
m=7;n=9;
for i=1:m*n
    if X(i)==1
        n1=n1+1;
    end
    if X(i)==2
        n2=n2+1;
    end
end
if n2>n1
    for i=1:m*n
        if X(i)==1
            Z(i)=2;
        end
        if X(i)==2
            Z(i)=1;
        end
    end
end
X=Z;
A=[];
B=nan(1,10000);
C=nan(1,10000);
D=nan(1,10000);
E=nan(1,10000);
b=0;
c=0;
for a=1:m*n
    if X(a)==1
        b=b+1;
    end
    if X(a)==2
        c=c+1;
    end
end
if b==c
    player=1;
else
    player=2;
end
for a=1:n
    for b=1:n
        for c=1:n
            for d=1:n
                Y=X;
                [Y,flag]=operatefour(Y,a);
                if flag==0
                    [Y,flag]=operatefour(Y,b);
                    if flag==0
                        [Y,flag]=operatefour(Y,c);
                        if flag==0
                            [Y,flag]=operatefour(Y,d);
                            p=a*1000+b*100+c*10+d;
                            if flag==0
                                A(p)=estimatefour(Y);
                            else A(p)=NaN;
                            end
                        else p=a*1000+b*100+c*10;
                            A(p)=NaN;
                        end
                    else p=a*1000+b*100;
                        A(p)=NaN;
                    end
                else p=a*1000;
                    A(p)=NaN;
                end
            end
        end
    end
end
if player==1
    for a=1:n
        if ~isnan(A(1000*a))
            for b=1:n
                if ~isnan(A(1000*a+100*b))
                    for c=1:n
                        if ~isnan(A(1000*a+100*b+10*c))
                            for d=1:n
                                if ~isnan(A(1000*a+100*b+10*c+d))
                                    B(1000*a+100*b+10*c)=min(A(1000*a+100*b+10*c+d),B(1000*a+100*b+10*c));
                                end
                            end
                            C(1000*a+100*b)=max(B(1000*a+100*b+10*c),C(1000*a+100*b));
                        end
                    end
                    D(1000*a)=min(C(1000*a+100*b),D(1000*a));
                end
            end
            E(1)=max(D(1000*a),E(1));
        end
    end
    for f=1:n
        if E(1)==D(1000*f)
            j=f;
            return;
        end
    end
end
if player==2
    for a=1:n
        if ~isnan(A(1000*a))
            for b=1:n
                if ~isnan(A(1000*a+100*b))
                    for c=1:n
                        if ~isnan(A(1000*a+100*b+10*c))
                            for d=1:n
                                if ~isnan(A(1000*a+100*b+10*c+d))
                                    B(1000*a+100*b+10*c)=max(A(1000*a+100*b+10*c+d),B(1000*a+100*b+10*c));
                                end
                            end
                            C(1000*a+100*b)=min(B(1000*a+100*b+10*c),C(1000*a+100*b));
                        end
                    end
                    D(1000*a)=max(C(1000*a+100*b),D(1000*a));
                end
            end
            E(1)=min(D(1000*a),E(1));
        end
    end
    for f=1:n
        if E(1)==D(1000*f)
            j=f;
            return;
        end
    end
end

function [Y,j]=operatefour(X,i)
b=0;
c=0;
m=7;n=9;
for a=1:m*n
    if X(a)==1
        b=b+1;
    end
    if X(a)==2
        c=c+1;
    end
end
if b==c
    player=1;
else player=2;
end

if X(1,i)==0
    j=0;
    Y=X;
    for a=m:-1:1
        if X(a,i)==0
           Y(a,i)=player;
           break
        end
    end
else j=1;
    Y=X;
end

a1=[];
a2=[];
A=X;
B=X;
for i=1:m*n
    if X(i)==1
        a1=[a1,i];
    end
    if X(i)==2
        a2=[a2,i];
    end
    if X(i)==0
        A(i)=1;
        B(i)=2;
    end
end
length1=length(a1);
length2=length(a2);

if length1>3
    b1=nchoosek(a1,4);
    [lb1,~]=size(b1);
    for i=1:lb1
        if b1(i,2)-b1(i,1)==m&&b1(i,3)-b1(i,2)==m&&b1(i,4)-b1(i,3)==m
            j=0;
            Y=X;
            return;
        end
        if b1(i,2)-b1(i,1)==(m-1)&&b1(i,3)-b1(i,2)==(m-1)&&b1(i,4)-b1(i,3)==(m-1)
            if mod(b1(i,1),7)>3||mod(b1(i,1),7)==0
                j=0;
            Y=X;
                return;
            end
        end
        if b1(i,2)-b1(i,1)==(m+1)&&b1(i,3)-b1(i,2)==(m+1)&&b1(i,4)-b1(i,3)==(m+1)&&mod(b1(i,1),m)>0&&mod(b1(i,1),m)<(m-2)
            j=0;
            Y=X;
            return;
        end
        if b1(i,2)-b1(i,1)==1&&b1(i,3)-b1(i,2)==1&&b1(i,4)-b1(i,3)==1&&mod(b1(i,1),m)>0&&mod(b1(i,1),m)<(m-2)
            j=0;
            Y=X;
            return;
        end
    end
end
if length2>3
    b2=nchoosek(a2,4);
    [lb2,~]=size(b2);
    for i=1:lb2
        if b2(i,2)-b2(i,1)==m&&b2(i,3)-b2(i,2)==m&&b2(i,4)-b2(i,3)==m
            j=0;
            Y=X;
            return;
        end
        if b2(i,2)-b2(i,1)==(m-1)&&b2(i,3)-b2(i,2)==(m-1)&&b2(i,4)-b2(i,3)==(m-1)
            if mod(b2(i,1),m)>3||mod(b2(i,1),m)==0
                j=0;
            Y=X;
            end
            return;
        end
        if b2(i,2)-b2(i,1)==(m+1)&&b2(i,3)-b2(i,2)==(m+1)&&b2(i,4)-b2(i,3)==(m+1)&&mod(b2(i,1),m)>0&&mod(b2(i,1),m)<(m-2)
            j=0;
            Y=X;
            return;
        end
        if b2(i,2)-b2(i,1)==1&&b2(i,3)-b2(i,2)==1&&b2(i,4)-b2(i,3)==1&&mod(b2(i,1),m)>0&&mod(b2(i,1),m)<(m-2)
            j=0;
            Y=X;
            return;
        end
    end
end

function x=estimatefour(X)
x=0;
a1=[];
a2=[];
A=X;
B=X;
m=7;n=9;
for i=1:m*n
    if X(i)==1
        a1=[a1,i];
    end
    if X(i)==2
        a2=[a2,i];
    end
    if X(i)==0
        A(i)=1;
        B(i)=2;
    end
end
T=zeros(m+1,n+2);
T(1:m,2:(n+1))=X;
T1=T;
T2=T;
for i=(m+2):(m*n+m+n)
    if T(i)==0
        if T(i+1)~=0||T(i+(m+1))~=0||T(i-(m+1))~=0
            T1(i)=1;
            T2(i)=2;
        end
    end
end
length1=length(a1);
length2=length(a2);
for i=1:m
    for j=1:(n-3)
        if X(i,j)==1&&X(i,j+1)==1&&X(i,j+2)==1&&X(i,j+3)==1
            x=1000000-length1;
            return
        end
        if X(i,j)==2&&X(i,j+1)==2&&X(i,j+2)==2&&X(i,j+3)==2
            x=-1000000+length2;
            return
        end
        if sum(A(i,j:j+3))==4&&sum(X(i,j:j+3))==3
            for k=0:3
                if X(i,j+k)==0
                    if i==m
                        x=x+10;
                    else if X(i+1,j+k)~=0
                            x=x+10;
                        else if mod(i,2)==0
                                x=x+20;
                            else x=x+50;
                            end
                        end
                    end
                    break;
                end
            end
        end
        if sum(B(i,j:j+3))==8&&sum(X(i,j:j+3))==6
            for k=0:3
                if X(i,j+k)==0
                    if i==m
                        x=x-10;
                    else if X(i+1,j+k)~=0
                            x=x-10;
                        else if mod(i,2)==0
                                x=x-50;
                            else x=x-20;
                            end
                        end
                    end
                    break;
                end
            end
        end
        if sum(A(i,j:j+3))==4&&sum(X(i,j:j+3))==2
            x=x+3;
        end
        if sum(B(i,j:j+3))==8&&sum(X(i,j:j+3))==4
            x=x-3;
        end
        if sum(A(i,j:j+3))==4&&sum(X(i,j:j+3))==1
            x=x+1;
        end
        if sum(B(i,j:j+3))==8&&sum(X(i,j:j+3))==2
            x=x-1;
        end
    end
end
for i=1:(m-3)
    for j=1:n
        if X(i,j)==1&&X(i+1,j)==1&&X(i+2,j)==1&&X(i+3,j)==1
            x=1000000-length1;
            return
        end
        if X(i,j)==2&&X(i+1,j)==2&&X(i+2,j)==2&&X(i+3,j)==2
            x=-1000000+length2;
            return
        end
        if sum(A(i:i+3,j))==4&&sum(X(i:i+3,j))==3
            x=x+10;
        end
        if sum(B(i:i+3,j))==8&&sum(X(i:i+3,j))==6
            x=x-10;
        end
        if sum(A(i:i+3,j))==4&&sum(X(i:i+3,j))==2
            x=x+3;
        end
        if sum(B(i:i+3,j))==8&&sum(X(i:i+3,j))==4
            x=x-3;
        end
        if sum(A(i:i+3,j))==4&&sum(X(i:i+3,j))==1
            x=x+1;
        end
        if sum(B(i:i+3,j))==8&&sum(X(i:i+3,j))==2
            x=x-1;
        end
    end
end
for i=1:(m-3)
    for j=1:(n-3)
        if X(i,j)==1&&X(i+1,j+1)==1&&X(i+2,j+2)==1&&X(i+3,j+3)==1
            x=1000000-length1;
            return
        end
        if X(i,j)==2&&X(i+1,j+1)==2&&X(i+2,j+2)==2&&X(i+3,j+3)==2
            x=-1000000+length2;
            return
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==3&&A(i,j)+A(i+1,j+1)+A(i+2,j+2)+A(i+3,j+3)==4
            for k=0:3
                if X(i+k,j+k)==0
                    if i+k==m
                        x=x+10;
                    else if X(i+k+1,j+k)~=0
                            x=x+10;
                        else if mod(i+k,2)==0
                                x=x+20;
                            else x=x+50;
                            end
                        end
                    end
                    break
                end
            end
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==6&&B(i,j)+B(i+1,j+1)+B(i+2,j+2)+B(i+3,j+3)==8
            for k=0:3
                if X(i+k,j+k)==0
                    if i+k==m
                        x=x-10;
                    else if X(i+k+1,j+k)~=0
                            x=x-10;
                        else if mod(i+k,2)==0
                                x=x-50;
                            else x=x-20;
                            end
                        end
                    end
                    break;
                end
            end
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==2&&A(i,j)+A(i+1,j+1)+A(i+2,j+2)+A(i+3,j+3)==4
            x=x+3;
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==4&&B(i,j)+B(i+1,j+1)+B(i+2,j+2)+B(i+3,j+3)==8
            x=x-3;
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==1&&A(i,j)+A(i+1,j+1)+A(i+2,j+2)+A(i+3,j+3)==4
            x=x+1;
        end
        if X(i,j)+X(i+1,j+1)+X(i+2,j+2)+X(i+3,j+3)==2&&B(i,j)+B(i+1,j+1)+B(i+2,j+2)+B(i+3,j+3)==8
            x=x-1;
        end
    end
end
for i=4:m
    for j=1:(n-3)
        if X(i,j)==1&&X(i-1,j+1)==1&&X(i-2,j+2)==1&&X(i-3,j+3)==1
            x=1000000-length1;
            return
        end
        if X(i,j)==2&&X(i-1,j+1)==2&&X(i-2,j+2)==2&&X(i-3,j+3)==2
            x=-1000000+length2;
            return
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==3&&A(i,j)+A(i-1,j+1)+A(i-2,j+2)+A(i-3,j+3)==4
            for k=0:3
                if X(i-k,j+k)==0
                    if i-k==m
                        x=x+10;
                    else
                        if X(i-k+1,j+k)~=0
                            x=x+10;
                        else
                            if mod(i-k,2)==0
                                x=x+20;
                            else x=x+50;
                            end
                        end
                    end
                    break
                end
            end
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==6&&B(i,j)+B(i-1,j+1)+B(i-2,j+2)+B(i-3,j+3)==8
            for k=0:3
                if X(i-k,j+k)==0
                    if i-k==m
                        x=x-10;
                    else
                        if X(i-k+1,j+k)~=0
                            x=x-10;
                        else
                            if mod(i-k,2)==0
                                x=x-50;
                            else x=x-20;
                            end
                        end
                    end
                    break
                end
            end
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==2&&A(i,j)+A(i-1,j+1)+A(i-2,j+2)+A(i-3,j+3)==4
            x=x+3;
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==4&&B(i,j)+B(i-1,j+1)+B(i-2,j+2)+B(i-3,j+3)==8
            x=x-3;
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==1&&A(i,j)+A(i-1,j+1)+A(i-2,j+2)+A(i-3,j+3)==4
            x=x+1;
        end
        if X(i,j)+X(i-1,j+1)+X(i-2,j+2)+X(i-3,j+3)==2&&B(i,j)+B(i-1,j+1)+B(i-2,j+2)+B(i-3,j+3)==8
            x=x-1;
        end
    end
end

function flg = checkbd(X,i,j)
jp = size(X,2)+1-j;
l1 = checkiar(X(i,:)',j);
l2 = checkiar(X(:,j),i);
l3 = checkiar(diag(X,j-i),min(i,j));
l4 = checkiar(diag(fliplr(X),jp-i), min(i,jp));
flg = any([l1 l2 l3 l4]>=4);

function l = checkiar(x,i)
x = [-1; x; -1];
i = i + 1;
ip = i;
while x(ip)==x(i)
    ip = ip-1;
end
im = i;
while x(im)==x(i)
    im = im + 1;
end
l = im-ip-1;

function [i,flg] = loc(X,j)
if all(X(:,j))
    flg = -1;
    i=0;
else
    flg = 1;
    i = find(X(:,j)==0,1,'last');
end

function drawbd(X)
[m,n] = size(X);
clf;
hold on;
t = linspace(0,2*pi);
r = 0.4;
s = sin(t)*r;
c = cos(t)*r;
col = 'wrb';
patch([0 1 1 0 0]*n,[0 0 1 1 0]*m,'k');
for i = 0:m
    plot([0 n],[i i],'r-','linewidth',2);
end
for j = 0:n
    plot([j j],[0 m],'r-','linewidth',2);
end
for i = 1:m
    for j = 1:n
        patch(j-.5+c,i-.5+s,col(X(m+1-i,j)+1));
    end
end
title('FOUR in a ROW','color','r','fontsize',18);
axis('off');