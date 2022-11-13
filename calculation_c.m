function A=calculation_c(x,y,rhop)
% Function that finds a self-cross point on the phase diagram
% INPUT [x(n*1),y(n*1)]
flag1=0;
Num1=zeros(1,2);
for k1=1:length(x)-1
    d1 = x(k1+1)-x(k1);
    switch flag1
        case 0
            if(d1<0)
                flag1 =1 ;
                Num1(1)=k1;
            end
        case 1
            if(d1>0)
                Num1(2)=k1;
                flag1 = 2;
            end
        otherwise
            break;
    end
end

if (flag1==0)
    A=-1;
    return;
end
xx1 = x(1:Num1(1));
mid = x(Num1(2));
p1=abs(abs(xx1)-abs(mid));
[~,b1]=min(p1);
X1 = [b1,Num1(1)];
yy2 = y(Num1(2):end);
mid = y(Num1(1));
p1=abs(abs(yy2)-abs(mid));
[~,b1]=min(p1);
X2= [Num1(2),b1+Num1(2)];

% Test for the correctness of the points
if y(X2(1))>y(X1(1))
    for Nt=X2(1):-1:1
        if y(Nt)<y(Nt-1)
            Ntt=Nt;
            X2(1)=Ntt;
            if y(Ntt)>y(X1(1)+10)
                A=-2;
                return
            else 
                break;
            
            end
            
        end
    end
end

A=b2_s_fp(rhop(X1(1):X1(2)),rhop(X2(1):X2(2)));





