clear all;global Mdx Mdt Mdx2 pmit1 t  sgma fr ro alfa 
StartTime=cputime;
ro=1;
n=50;
alfa=0.75;
a=0.75;
sgma=1;
%---------------------------
na=n+1;nt=na
    hx=0.2;
 for i=1:nt
    x(i)=-5+(i-1)*hx;
    %x(i)=1;
 end
    x'
ht=0.1;
for j=1:nt
t(j)=(j-1)*ht;
%t(j)=ht;
end
t'

for j=1:nt
for i=1:na
 xx=x(i);
 tt=t(j);
L0=1;L1=(tt^alfa)/(gamma(alfa+1));
L2=(tt^(2*alfa))/(gamma(2*alfa+1));
             L3=(tt^(3*alfa))/(gamma(3*alfa+1));
              f0=1/(1+exp(-xx/sqrt(2)));
             f1=0.5*(1-2*a)*exp(-xx/sqrt(2))*(1+exp(-xx/sqrt(2)))^(-2);
 f2=0.25*(1-2*a)*(1-2*a)*exp(-xx/sqrt(2))*(1-exp(-xx/sqrt(2)))*(1+exp(-xx/sqrt(2)))^(-3);
         f3=((3*exp(-2^(1/2)*xx)*exp(-(2^(1/2)*xx)/2)*(2*a - 1)*(a/2 - 1/4))/(exp(-(2^(1/2)*xx)/2) + 1)^4 - (3*exp(-2^(1/2)*xx)*(2*a - 1)*(a/2 - 1/4))/(2*(exp(-(2^(1/2)*xx)/2) + 1)^3) + (9*exp(-2^(1/2)*xx)*(2*a - 1)*(a/2 - 1/4)*(exp(-(2^(1/2)*xx)/2) - 1))/(2*(exp(-(2^(1/2)*xx)/2) + 1)^4) - (exp(-(2^(1/2)*xx)/2)*(2*a - 1)*(a/2 - 1/4)*(exp(-(2^(1/2)*xx)/2) - 1))/(2*(exp(-(2^(1/2)*xx)/2) + 1)^3) - (6*exp(-2^(1/2)*xx)*exp(-(2^(1/2)*xx)/2)*(2*a - 1)*(a/2 - 1/4)*(exp(-(2^(1/2)*xx)/2) - 1))/(exp(-(2^(1/2)*xx)/2) + 1)^5)+(a+1)*(f1*f1)*gamma(2*alfa+1)/(gamma(alfa+1)*gamma(alfa+1))+2*(a+1)*(f2*f0)-(a*f2)-(3*f0)*(f1*f1)*gamma(2*alfa+1)/(gamma(alfa+1)*gamma(alfa+1))-(3*(f0*f0)*f2);
           ua(i,j)=L0*f0+L1*f1+L2*f2+L3*f3;

end
end

for i=1:nt 
    for j=1:nt
ue(i,j)=Exact(x(i),t(j),nt,na,a,alfa);
    end
end
    uer=abs(ua-ue);
    fprintf('    x   t       ua         ue       er \n ')
for i=1:nt
    for j=1:nt
fprintf('%9.6f   %9.6f %9.6f %9.6f %10.5e \n',x(i),t(j),ua(i,j),ue(i,j),uer(i,j))
 end
end 
    uerr=norm2(uer,nt)
    fprintf('norm err  u   = %10.5e \n',uerr);
    %z=three_dim_grpah(x,t,ue,ua)
    z=three_dim_grpah_error(x,t,ua)
   % plot(t,ua)
    EndTime=cputime;
      tof=EndTime-StartTime;
      fprintf('excution time of  lfopc is %5.2f (second)\n',tof)
    fprintf('--------------------------------------------------------------\n',tof)
    Minute1=tof/60;
      fprintf('excution time  in minutes  is %15.8f (Minutes)\n',Minute1)
    hold on
