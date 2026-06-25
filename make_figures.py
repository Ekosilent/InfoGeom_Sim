import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.fft import dstn, idstn
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
os.makedirs("figs", exist_ok=True)

def make_solver(NX,NY,XMIN=-10,XMAX=10,YMIN=-5,YMAX=5):
    Lx,Ly=XMAX-XMIN,YMAX-YMIN
    dx=Lx/(NX+1); dy=Ly/(NY+1)
    x=XMIN+dx*np.arange(1,NX+1); y=YMIN+dy*np.arange(1,NY+1)
    X,Y=np.meshgrid(x,y)
    ii=np.arange(1,NX+1); jj=np.arange(1,NY+1)
    lx=(2*np.cos(np.pi*ii/(NX+1))-2)/dx**2
    ly=(2*np.cos(np.pi*jj/(NY+1))-2)/dy**2
    LAM=ly[:,None]+lx[None,:]
    return x,y,dx,dy,X,Y,(lambda f: idstn(dstn(f,type=1)/LAM,type=1))

# FIG 1
NX,NY=300,150
x,y,dx,dy,X,Y,solve=make_solver(NX,NY)
A,B,ALPHA=1.0,0.5,1.0
src=[(3.0,0.5,1.0,0.6),(-2.0,-1.0,0.6,1.2),(0.5,2.0,0.4,0.8)]
rho=np.zeros_like(X)
for sx,sy,amp,sig in src: rho+=amp*np.exp(-((X-sx)**2+(Y-sy)**2)/sig**2)
gy,gx=np.gradient(rho,dy,dx)
I=A*rho+B*(gx**2+gy**2)
Phi=solve(ALPHA*I); h00=-2*Phi
gy2,gx2=np.gradient(Phi,dy,dx)
Fx=RegularGridInterpolator((y,x),-gx2,bounds_error=False,fill_value=0.0)
Fy=RegularGridInterpolator((y,x),-gy2,bounds_error=False,fill_value=0.0)
def rhs(t,s):
    xp,yp,vx,vy=s; return [vx,vy,float(Fx((yp,xp))),float(Fy((yp,xp)))]
trajs=[solve_ivp(rhs,(0,45),[-9.0,y0,0.6,0.0],t_eval=np.linspace(0,45,900),
        method="RK45",rtol=1e-6,atol=1e-6).y for y0 in (-2.5,-1.0,0.5,2.0)]
fig,ax=plt.subplots(figsize=(8,4))
im=ax.imshow(h00,extent=[x.min(),x.max(),y.min(),y.max()],origin="lower",cmap="inferno",aspect="auto")
ax.contour(X,Y,h00,levels=12,colors="white",alpha=0.25,linewidths=0.4)
for t in trajs: ax.plot(t[0],t[1],color="#ffd27a",lw=1.2)
for sx,sy,_,_ in src: ax.plot(sx,sy,"o",mfc="none",mec="#ffd27a",ms=8,mew=1.3)
ax.set_xlabel("x"); ax.set_ylabel("y")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.02); cb.set_label(r"$h_{00}=-2\Phi$")
plt.tight_layout(); plt.savefig("figs/field.pdf"); plt.close()

# FIG 2 MMS
def mms(NX,NY):
    x,y,dx,dy,X,Y,solve=make_solver(NX,NY)
    Lx,Ly=20.0,10.0
    u=(X+10)/Lx; v=(Y+5)/Ly
    Pex=np.sin(np.pi*u)*np.sin(np.pi*v)
    f=-(np.pi**2/Lx**2+np.pi**2/Ly**2)*Pex
    Pn=solve(f)
    return dx,np.sqrt(np.mean((Pn-Pex)**2))/np.sqrt(np.mean(Pex**2)),X,Y,Pn-Pex
hs=[];es=[]
for N in [(50,25),(100,50),(200,100),(400,200),(800,400)]:
    h,e,_,_,_=mms(*N); hs.append(h); es.append(e)
hs=np.array(hs);es=np.array(es)
slope=np.polyfit(np.log(hs),np.log(es),1)[0]
_,_,Xe,Ye,emap=mms(200,100)
fig,(a1,a2)=plt.subplots(1,2,figsize=(9,3.5))
pc=a1.pcolormesh(Xe,Ye,emap,cmap="RdBu_r",shading="auto")
a1.set_title(r"Errore $\Phi_{\rm num}-\Phi_{\rm ex}$  (200$\times$100)",fontsize=10)
a1.set_xlabel("x"); a1.set_ylabel("y"); fig.colorbar(pc,ax=a1,fraction=0.046,pad=0.02)
a2.loglog(hs,es,"o-",color="#1f4e8c",label="errore L2 relativo")
a2.loglog(hs,es[0]*(hs/hs[0])**2,"--",color="gray",label=r"riferimento $O(h^2)$")
a2.set_xlabel("passo griglia $h$"); a2.set_ylabel("errore L2 rel.")
a2.set_title(f"Convergenza: pendenza $\\approx${slope:.2f}",fontsize=10)
a2.legend(fontsize=8); a2.grid(True,which="both",alpha=0.3)
plt.tight_layout(); plt.savefig("figs/mms.pdf"); plt.close()
print("OK slope",round(slope,3),"h00max",round(h00.max(),3))
