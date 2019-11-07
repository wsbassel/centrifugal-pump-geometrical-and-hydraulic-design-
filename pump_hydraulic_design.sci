function [fi]=nasa_criterion_pump_design(r1, r2, wedth, beta2d ,z2,Q,rpm) 


nmax=2001 // number of points of pump blade curve
x1=0
y1=r1
u1=2*%pi*r1*rpm/60 //periphery speed at inlet
cm1=Q/(%pi*r1*width)// inlet axial absolute velocity
beta1=atan(cm1/u1) // inlet flow angle radians
beta2=beta2d*%pi/180  // exit flow angle radians
fi=-.1  // initial guess , fi angle between r1 and r2
delta=.1    // inital guess , delta angle between L12 and r2
L12=r1   // initail guess 
f=zeros(3,1);
f(1)=-beta1-fi-2*delta+beta2;
f(2)-r1*sin(fi)+L12*sin(delta);
f(3)=-r2^2+r1^2+L12^2-2*r1*L12*cos(%pi-delta-fi);
//START NEWTON RAPHSON FOR BLADE GEOMETRY
iter =1
while max(abs(f))> 1e-7
j=zeros(3,3);
j(1,1)=-1;
j(1,2)=-2;
j(2,1)=-r1*cos(fi);
j(2,2)=L12*cos(delta);
j(2,3)=sin(delta)
j(3,1)=-2*r1*L12*sin(%pi-delta-fi);
j(3,2)=-2*r1*L12*sin(%pi-delta-fi);
j(3,3)=-2*r1*cos(%pi-delta-fi)+2*L12;
del=j\f;
nf=10
fi=fi-del(1)/nf;
delta=delta-del(2)/nf;
L12=L12-del(3)/nf;
f=zeros(3,1);
f(1)=-beta1-fi-2*delta+beta2;
f(2)-r1*sin(fi)+L12*sin(delta);
f(3)=-r2^2+r1^2+L12^2-2*r1*L12*cos(%pi-delta-fi);
iter=iter+1
end
disp f_geomet
disp (f)


x2=r2*sin(fi)
y2=r2*cos(fi)
L12_r=sqrt((x1-x2)^2+(y1-y2)^2)
ang1=%pi/2-beta1-fi-delta
ang2=%pi/2-beta2+delta
dif_ang12=ang2+ang1
for n=1:nmax
bta(n)=ang1-(n-1)*dif_ang12/(nmax-1);
end
dyL=L12/(nmax-1)
yL(1)=0
for n= 2 :nmax
yL(n)=(n-1)*dyL;
end
dxL(1)=0
xL(1)=0
for n=2:nmax    
dxL(n)=dyL*tan(bta(n));
xL(n)=xL(n-1)+dxL(n);
end
for n=1:nmax
d_lenth(n)=sqrt(dxL(n)^2+dyL^2);
end
lenth=sum(d_lenth)
disp lenth
disp (lenth)
dx_i=(x2-x1)/(nmax-1)
dy_i=(y2-y1)/(nmax-1)
for n=1:nmax
x_i(n)=x1+(n-1)*dx_i;
y_i(n)=y1+(n-1)*dy_i;
end
m12=(y2-y1)/(x2-x1);
ang_i= atan(m12)
for n=1 :nmax
xr(n)=x_i(n)+xL(n)*sin(ang_i);
yr(n)=y_i(n)-xL(n)*cos(ang_i);
end
for n=1:73
teta_c(n)=(n-1)*5/180*%pi
xr1(n)=r1*cos(teta_c(n))
yr1(n)=r1*sin(teta_c(n))
xr2(n)=r2*cos(teta_c(n))
yr2(n)=r2*sin(teta_c(n))
end
for n=1 : 21
xr_r(n)=xr(100*(n-1)+1)
yr_r(n)=yr(100*(n-1)+1)
end
coord(:,1)=xr_r
coord(:,2)=yr_r
disp 'coordentes x and y'
disp (coord)
thick=.02
for n=2:21
R(n)=sqrt(xr_r(n)^2+yr_r(n)^2)
atang(n)=atan(yr_r(n)/xr_r(n))
teta(n)=atang(n)+2*%pi/z2     
xr_r2(n)=R(n)*cos(teta(n))
yr_r2(n)=R(n)*sin(teta(n))
end
xr_r2(1)=r1*cos(2*%pi/z2)
yr_r2(1)=r1*sin(2*%pi/z2)
figure
plot(xr1,yr1, xr2, yr2, xr_r ,yr_r)
DELTA=r1/r2
X=lenth/(2*r2)
um_sigma=1+(1+.6*sin(beta2))/(z2*(1+DELTA)*X^2+.25*(1-DELTA)^2)^.5
sigma=1/um_sigma

area_out=%pi*2*r2*width
cm2=Q/area_out
beta2=beta2d/180*%pi
u2=2*%pi*r2*rpm/60
slip=u2*(1-sigma)
cu2th=u2+cm2/tan(beta2) 
cu2=-slip+cu2th     //u2-wu2th-slip
disp slip
disp (slip)
disp cu2th
disp (cu2th)
disp cu2
disp (cu2)
disp u2
disp (u2)
Hth=u2*cu2/9.81*3.2803
Qgm=Q*1.58503e4
con=(32.12)^(.75)
conns=%pi*rpm/30*Qgm^.5   //H^3.75
ns=.1
eta=.8
H=.8*Hth
// START NEWTON RAPHSON FOR PUMP HYDRAULIC DESIGN
f2=zeros(3,1)
f2(1)=-ns+%pi*rpm/30/(32.12)^.75*H^(-.75)*Qgm^.5
//1Jhyd,design = 0.41989 + 2.1524 Ns - 3.1434 Ns 2 + 1.5673 Ns 3 (24)
if ns < .8 then
f2(2)=-eta+(.41989+2.1524*ns-3.1434*ns^2+1.5673*ns^3)
elseif ns >.8 & ns< 8
    f2(2)=-eta+1.02-.12*ns
else 
    eta=.04
    end
f2(3)=-H+Hth*eta
H=100
iter2=1
while max(abs(f2))>1e-7
j2=zeros(3,3)
j2(1,1)=-1
j2(1,2)=(%pi*rpm/30/(32.12)^.75*(-.75)*H^(-1.75)*Qgm^.5)
if ns <.8 then
j2(2,1)=(2.1524-3.1434*2*ns+1.5673*3*ns^2)

elseif ns<.8 & ns>8
    j2(2,1)=-.12
    

   
    end
j2(2,3)=-1
j2(3,2)=-1
j2(3,3)=Hth
delta2=j2\f2
ns=ns-delta2(1)/100
H=H-delta2(2)/100
eta=eta-delta2(3)/100
f2=zeros(3,1)
f2(1)=-ns+%pi*rpm/30/(32.12)^.75*H^(-.75)*Qgm^.5
//1Jhyd,design = 0.41989 + 2.1524 Ns - 3.1434 Ns 2 + 1.5673 Ns 3 (24)
if ns < .8 then
f2(2)=-eta+(.41989+2.1524*ns-3.1434*ns^2+1.5673*ns^3)
elseif ns >.8 & ns< 8
    f2(2)=-eta+1.02-.12*ns
else 
    eta=.04
    end
f2(3)=-H+Hth*eta
iter2=iter2+1
end
disp f2
disp (f2)
disp ns
disp (ns)
disp 'Head feet'
disp (H)
disp eta
disp (eta)

if ns > 3.5 then
    disp ' pumps - (centrifugal or axial)- do not operate with ns greater than 3.5 - the result is mathematical calculation has no physical meaning '
    disp ' do not use if ns greater tham 3.5 '
end

endfunction 
r1=.05
r2=.25
z2=20
width=.02;
beta2d=60
Q=.01
rpm=5000
[fi]=nasa_criterion_pump_design(r1, r2, width, beta2d,z2,Q,rpm) 
