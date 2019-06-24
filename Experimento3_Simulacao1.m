warning ("off");
randn("state",[5]); qtd_redes=2000;
tolerance=10**(-3) # 2 casas decimais além da última casa decimal das observações da rede # 
m=20; n=10;com_teste_qui=1;valor_qui=20.48; valor_qui_min=3.25;  #teste bicaudal significancia 5%#
cont_max_v1=0;cont_mean_v1=0;cont_std_v1=0;cont_min_v1=0;cont_range_v1=0;
cont_max_DPx=0;cont_mean_DPx=0;cont_std_DPx=0;cont_min_DPx=0;cont_range_DPx=0;
cont_max_DPv=0;cont_mean_DPv=0;cont_std_DPv=0;cont_min_DPv=0;cont_range_DPv=0;
LInfTotal=0;IterTotal=0;
qui_sim_usu=0;qui_sim_prop=0;

L_perf=[163854.9;6446.2;57037.0;126209.5;101128.6;296885.8;398014.4;60449.1;173710.4;167264.2;
110227.2;155928.2;52875.0;62904.2;3889.5;42705.7;98891.2;115779.2;113222.5;46428.8];

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

for sss = 1:qtd_redes

L=L_perf;  
d=[49 41 38 34 22 13 23 48 15 24 62 49 35 43 20 28 19 39 27 21]; 
dp=1.0.*sqrt(d);
for q=1:m
    do a=randn(1);
    until (a<=3)
    aleat(q)=dp(q)*a;
    L(q)=L(q)+aleat(q);
end

##################################################### calculo do MRA
p=ones(m,1);
P=diag(p); 

I = eye(m);
O = zeros(m,2*n);
A1 = [A -A -I I zeros(m,1)];
A2 = [O P -P (-1)*ones(m,1)];
A3 = [O -P P (-1)*ones(m,1)];
A4 = cat(1,A1, A2, A3);
L1 = cat(1,L,zeros(2*m,1));
c=zeros(1,2*(n+m)+1);c(1,2*(n+m)+1)=1;

ctype=[];
for i=1:(m)
  ctype = [ctype, "S"];
 end
 for i=1:(2*m)
  ctype = [ctype, "U"];
 end
vartype=[];
for i=1:(2*(n+m)+1)
  vartype = [vartype, "C"];
 end
param.dual=1; 
param.lpsolver=1; 

[xopt, fopt, erro, extra] = glpk (c, A4, L1, lb=[], ub=[], ctype, vartype, s=1, param);
xopt;fopt;erro;   extra.status;
for i=1:n
    x(i)=xopt(i)-xopt(n+i);
end
for i=1:m
    v(i)=xopt(2*n+i)-xopt(2*n+i+m);
end
v;
LInf=max(abs(v));
LInfTotal=LInfTotal+LInf;
##################################################### iteracoes do MMQ
p=(LInf^(-2)).*ones(m,1);
P=diag(p);
x=inv(A'*P*A)*A'*P*L;
v=A*x-L;

  for k = 1:10000
      pant=p;
      for i = 1:m
          if abs(v(i))-LInf>tolerance
            p(i)=p(i)*(abs(v(i))/LInf);
          end;
      end;
      P=diag(p);
      x=(A'*P*A)\A'*P*L;
      v=A*x-L;
      if (pant==p),
        if k>130 ####################isso eh soh pra saber se passou de 130 iteracoes do MMQ em algum cenario
          k end
        IterTotal=IterTotal+k; break; end; %stop if necessary
      
  end;
p_iter=p;
######################################################## 1) MMQ com mod. estocastico usual
#"MMQ com p=1/d"
var=dp.^2;
p=1./var;
ELb=diag(var);
P=diag(p);
x=inv(A'*P*A)*A'*P*L;
v=A*x-L;

if com_teste_qui ==1
    sigma2_pos=(v'*P*v)/(m-n);
    QuiQuadCalc=sigma2_pos*(m-n);
    LimSupQuiQuad=valor_qui;
    LimInfQuiQuad=valor_qui_min;
    if QuiQuadCalc<LimSupQuiQuad
      if QuiQuadCalc>LimInfQuiQuad
        qui_sim_usu=qui_sim_usu+1;
        end end
end

sigma2_pos=(v'*P*v)/(m-n);
Qx=inv(A'*P*A);
Ex=Qx.*sigma2_pos;
ELa=A*Ex*A';
Ev=(ELb.*sigma2_pos)-ELa;

###output
v1=abs(v);
diagDPx=sqrt(diag(Ex)); 
diagDPv=sqrt(abs(diag(Ev)));

max_v1=max(v1);
mean_v1=mean(v1);
std_v1=std(v1);
max_DPx=max(diagDPx);
mean_DPx=mean(diagDPx);
std_DPx=std(diagDPx);
max_DPv=max(diagDPv);
mean_DPv=mean(diagDPv);
std_DPv=std(diagDPv);
min_v1=min(v1);
range_v1=max_v1-min_v1;
min_DPx=min(diagDPx);
range_DPx=max_DPx-min_DPx;
min_DPv=min(diagDPv);
range_DPv=max_DPv-min_DPv;
######################################################## 2) MMQ com mod. estocastico proposto
#"MMQ com mod. estocastico proposto"
p=p_iter;
daux=1./p;
ELb=diag(daux);
P=diag(p);
x=inv(A'*P*A)*A'*P*L;
v=A*x-L;

if com_teste_qui ==1
    sigma2_pos=(v'*P*v)/(m-n);
    QuiQuadCalc=sigma2_pos*(m-n);
    LimSupQuiQuad=valor_qui;
    LimInfQuiQuad=valor_qui_min;
    if QuiQuadCalc<LimSupQuiQuad
      if QuiQuadCalc>LimInfQuiQuad
        qui_sim_prop=qui_sim_prop+1;
        end end
end

sigma2_pos=(v'*P*v)/(m-n);
Qx=inv(A'*P*A);
Ex=Qx.*sigma2_pos;
ELa=A*Ex*A';
Ev=(ELb.*sigma2_pos)-ELa;

###output
v1=abs(v);
diagDPx=sqrt(diag(Ex)); 
diagDPv=sqrt(abs(diag(Ev)));

max_v1_prop=max(v1);
mean_v1_prop=mean(v1);
std_v1_prop=std(v1);
max_DPx_prop=max(diagDPx);
mean_DPx_prop=mean(diagDPx);
std_DPx_prop=std(diagDPx);
max_DPv_prop=max(diagDPv);
mean_DPv_prop=mean(diagDPv);
std_DPv_prop=std(diagDPv);
min_v1_prop=min(v1);
range_v1_prop=max_v1_prop-min_v1_prop;
min_DPx_prop=min(diagDPx);
range_DPx_prop=max_DPx_prop-min_DPx_prop;
min_DPv_prop=min(diagDPv);
range_DPv_prop=max_DPv_prop-min_DPv_prop;
###comparacoes
  if max_v1_prop<max_v1
  cont_max_v1=cont_max_v1+1;
  end
  if mean_v1_prop<mean_v1
  cont_mean_v1=cont_mean_v1+1;
  end
  if std_v1_prop<std_v1
  cont_std_v1=cont_std_v1+1;
  end
  if max_DPx_prop<max_DPx
  cont_max_DPx=cont_max_DPx+1;
  end
  if mean_DPx_prop<mean_DPx
  cont_mean_DPx=cont_mean_DPx+1;
  end
  if std_DPx_prop<std_DPx
  cont_std_DPx=cont_std_DPx+1;
  end
  if max_DPv_prop<max_DPv
  cont_max_DPv=cont_max_DPv+1;
  end
  if mean_DPv_prop<mean_DPv
  cont_mean_DPv=cont_mean_DPv+1;
  end
  if std_DPv_prop<std_DPv
  cont_std_DPv=cont_std_DPv+1;
  end
  if min_v1_prop<min_v1
  cont_min_v1=cont_min_v1+1;
  end
  if range_v1_prop<range_v1
  cont_range_v1=cont_range_v1+1;
  end
  if min_DPx_prop<min_DPx
  cont_min_DPx=cont_min_DPx+1;
  end
  if range_DPx_prop<range_DPx
  cont_range_DPx=cont_range_DPx+1;
  end
  if min_DPv_prop<min_DPv
  cont_min_DPv=cont_min_DPv+1;
  end
  if range_DPv_prop<range_DPv
  cont_range_DPv=cont_range_DPv+1;
  end
 end
mediaLInf=LInfTotal/qtd_redes
mediaIter=IterTotal/qtd_redes
qui_sim_usu
qui_sim_prop
"Contagem ME proposto foi menor (logo MELHOR):"
cont_max_v1
cont_mean_v1
cont_std_v1
cont_max_DPx
cont_mean_DPx
cont_std_DPx
cont_max_DPv
cont_mean_DPv
cont_std_DPv

#cont_min_v1
#cont_range_v1
#cont_min_DPx
#cont_range_DPx
#cont_min_DPv
#cont_range_DPv