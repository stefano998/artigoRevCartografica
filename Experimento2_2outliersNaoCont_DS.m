warning ("off");
W=3.29;
if 1
mem=[1	1	4
2	1	5
3	1	7
4	1	8
5	4	5
6	4	6
7	4	8
8	5	6
9	5	7
10	6	7
11	6	8
12	7	8
13	12	16
14	12	17
15	13	14
16	13	20
17	14	15
18	14	18
19	14	20
20	15	16
21	15	18
22	15	19
23	16	17
24	16	19];
  end
for s=1:4
    if s==1
    ini_int_erro=3; fim_int_erro=6;
  elseif s==2
    ini_int_erro=6; fim_int_erro=12;
  elseif s==3
    ini_int_erro=12; fim_int_erro=25;
  elseif s==4
    ini_int_erro=25; fim_int_erro=100;end
#rand("state",[3]);randn("state",[3]);
qtd_itr=2000; cont_acerto=0;
m1=20; n=10;m=m1;
A=[1	0	0	0	0	0	0	0	0	0;
-1	0	0	0	0	0	0	0	0	1;
0	0	0	0	0	1	0	0	0	-1;
0	0	0	0	0	1	-1	0	0	0;
0	0	0	0	0	0	1	0	0	0;
0	0	0	0	0	0	-1	1	0	0;
0	0	0	0	0	0	0	1	0	0;
0	0	0	0	0	0	0	1	-1	0;
-1	0	0	0	0	0	0	0	1	0;
0	0	0	0	0	0	0	0	1	-1;
0	0	0	0	0	-1	0	0	1	0;
0	0	0	0	1	0	0	0	0	-1;
-1	1	0	0	0	0	0	0	0	0;
0	-1	1	0	0	0	0	0	0	0;
0	0	-1	1	0	0	0	0	0	0;
0	0	0	-1	1	0	0	0	0	0;
0	0	0	0	1	-1	0	0	0	0;
-1	0	1	0	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	-1;
0 1 0 0 0 0 0 0 0 -1];

d=[49 41 38 34 22 13 23 48 15 24 62 49 35 43 20 28 19 39 27 21];
dp=1.*sqrt(d);p=dp.^(-2);
P=diag(p);

for q=1:qtd_itr
    #Entrada L (rede perfeita) e simulaçoes em L
L=[163854.9;6446.2;57037.0;126209.5;101128.6;296885.8;398014.4;60449.1;173710.4;167264.2;
110227.2;155928.2;52875.0;62904.2;3889.5;42705.7;98891.2;115779.2;113222.5;46428.8];
    m=m1;   
    for z=1:m
         do a=randn(1);
         until (a<=3)
       Laleat(z,1)=L(z)+dp(z)*a;
    end

    
    for k=1:100   
    Lgross=Laleat;m=m1;   
    
    f=randi([1 24]);
    j=mem(f,2); p=mem(f,3);
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(j)=L(j)+choice*dp(j)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(p)=L(p)+choice*dp(p)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    x=inv(A'*P*A)*A'*P*Lgross;
    v=A*x-Lgross;
    Ev=inv(P)-A*inv(A'*P*A)*A';
    ww=zeros(m, 1);
    for i=1:m
      ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;
    if max(ww)>W
      if max(ww)==ww(j) | max(ww)==ww(p)
          if max(ww)==ww(j)
            out1=j;out2=p; 
            if p>j
              out2=out2-1;end
          else out1=p; out2=j;
            if j>p
              out2=out2-1; end
          end
                
        m=m-1;
        A2=A;Lgross2=Lgross;P2=P;
        A2(out1,:)=[];Lgross2(out1,:)=[];P2(out1,:)=[];P2(:,out1)=[];
        x=inv(A2'*P2*A2)*A2'*P2*Lgross2;
        v=A2*x-Lgross2;
        Ev=inv(P2)-A2*inv(A2'*P2*A2)*A2';
        ww=zeros(m, 1);
        for i=1:m
          ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;
        if max(ww)>W
         if max(ww)==ww(out2)
       
        m=m-1;
        A3=A2;Lgross3=Lgross2;P3=P2;
        A3(out2,:)=[];Lgross3(out2,:)=[];P3(out2,:)=[];P3(:,out2)=[];
        x=inv(A3'*P3*A3)*A3'*P3*Lgross3;
        v=A3*x-Lgross3;
        Ev=inv(P3)-A3*inv(A3'*P3*A3)*A3';
        ww=zeros(m, 1);
        for i=1:m
          ww(i,1)=abs(v(i,1)/sqrt(abs(Ev(i,i)))); 
          if ww(i)>10**5     
            ww(i)=3; end     
          end;   
          
        if max(ww)<W
          cont_acerto=cont_acerto+1;end;
end
end
end
end
    
end
end
ini_int_erro
cont_acerto
TS=(cont_acerto/200000)*100
end
