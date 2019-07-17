t0 = clock ();
W=3.29;
for s=1:4
#rand("state",[0]);randn("state",[0]);
  if s==1
    ini_int_erro=3; fim_int_erro=6;
  elseif s==2
    ini_int_erro=6; fim_int_erro=12;
  elseif s==3
    ini_int_erro=12; fim_int_erro=25;
  elseif s==4
    ini_int_erro=25; fim_int_erro=100;end
qtd_itr=2000; cont_idt_correto=0;cont_idt_erro2=0;
%qtd_itr se refere aos cenarios sem outliers. para cada um deles, serão 100 cenarios com outliers.
%assim, a qtd de iteracoes serah qtd_itr*100

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
#p=ones(1,m);          ##desmarcar aqui para testar DS com pesos unitarios#
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
    j=randi([1 m]);
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(j)=L(j)+choice*dp(j)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    x=inv(A'*P*A)*A'*P*Lgross;
    v=A*x-Lgross;
    Ev=inv(P)-A*inv(A'*P*A)*A';
    ww=zeros(m, 1);r=zeros(m, 1);
    for i=1:m
      ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));
      end;
    if max(ww)>W
    if max(ww)==ww(j)
        cont_idt_correto=cont_idt_correto+1; 
        
        m=m-1;
        A2=A;Lgross2=Lgross;P2=P;
        A2(j,:)=[];Lgross2(j,:)=[];P2(j,:)=[];P2(:,j)=[];
        x=inv(A2'*P2*A2)*A2'*P2*Lgross2;
        v=A2*x-Lgross2;
        Ev=inv(P2)-A2*inv(A2'*P2*A2)*A2';
        ww=zeros(m, 1);
        for i=1:m
          ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;
        if max(ww)>W
          cont_idt_erro2=cont_idt_erro2+1;end;
        end; end;         
        
   end;
end
ini_int_erro
acerto_real=cont_idt_correto-cont_idt_erro2 
TS=(acerto_real/200000)*100
end
tempo_medio = etime (clock (), t0)/800000