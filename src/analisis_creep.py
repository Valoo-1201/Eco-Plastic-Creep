"""Análisis de creep de filamento de PET reciclado: deformación, ajuste de modelos (Burgers, Findley, Prony),
validación predictiva y gráficas. Uso:  python src/analisis_creep.py
"""
import numpy as np, json, csv
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
# Raíz del repositorio (este script vive en src/)
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S=os.path.join(ROOT,'resultados')
DATOS=os.path.join(ROOT,'datos')
BLUE,ORANGE,GRAY,INK,MUTED='#2a78d6','#eb6834','#b9b8b2','#0b0b0b','#52514e'
plt.rcParams.update({'font.size':10,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
  'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e5e0','grid.linewidth':0.6,'figure.dpi':150})
g=9.81; d_mm=0.75; A=np.pi*(d_mm*1e-3/2)**2
def load(f):
    d=np.genfromtxt(f,delimiter=',',names=True); return d
def clean(d, hu_ref, hl_ref, tol=0.8):
    L=d['yl']-d['yu']
    ok=~np.isnan(L)&(np.abs(d['hu']-hu_ref)<tol)&(np.abs(d['hl']-hl_ref)<tol)
    return L, ok
def logbin(t,y,t0,t1,n):
    e=np.logspace(np.log10(t0),np.log10(t1),n+1); tb=[];yb=[];sb=[]
    for a,b in zip(e[:-1],e[1:]):
        m=(t>=a)&(t<b)
        if m.sum()>=2: tb.append(np.median(t[m])); yb.append(np.median(y[m])); sb.append(np.std(y[m])/np.sqrt(m.sum())+1e-6)
    return np.array(tb),np.array(yb),np.array(sb)
def burgers(t,e0,e2,tau,r): return e0+e2*(1-np.exp(-t/tau))+r*t
def findley(t,e0,a,n): return e0+a*t**n
def prony2(t,e0,a1,t1,a2,t2,r): return e0+a1*(1-np.exp(-t/t1))+a2*(1-np.exp(-t/t2))+r*t
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
out={}

# ---------------- PRUEBA 2 ----------------
d=load(f'{DATOS}/prueba2/seguimiento_crudo_px.csv'); t=d['t']
L,ok=clean(d,11.57,12.75)
ref=(t>=2.0)&(t<=4.2); L0=L[ref].mean(); sL0=L[ref].std()
t_load=4.4                                   # inicio de descarga de la mano (s de video)
PXMM=2.856                                   # px/mm vertical (regla de aluminio, FFT de marcas mm)
sig=1.0*g/A
tr=t-t_load; eps=(L-L0)/L0
# serie limpia para exportar (1 dato/s, mediana)
rows=[]
for s in np.arange(0,int(tr.max())):
    m=ok&(tr>=s)&(tr<s+1)
    if m.sum(): rows.append([s+0.5, np.median(L[m]), np.median(L[m])/PXMM, np.median(eps[m])*100])
with open(f'{S}/prueba2_deformacion_1s.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_s','L_px','L_mm','deformacion_%']); w.writerows(np.round(rows,5))
# ajuste: desde 8 s (peso ya liberado) — bins log
fit=ok&(tr>=8)
tb,eb,sb=logbin(tr[fit],eps[fit],8,tr.max(),40)
pB,cB=curve_fit(burgers,tb,eb,p0=[0.004,0.002,100,1e-6],sigma=sb,bounds=([0,0,1,0],[0.05,0.05,1e4,1e-3]),maxfev=20000)
pF,cF=curve_fit(findley,tb,eb,p0=[0.004,1e-4,0.3],sigma=sb,bounds=([0,0,0.01],[0.05,0.05,1]),maxfev=20000)
pP,cP=curve_fit(prony2,tb,eb,p0=[0.004,0.001,10,0.002,300,1e-7],sigma=sb,bounds=([0,0,1,0,20,0],[0.05,0.05,200,0.05,1e4,1e-3]),maxfev=50000)
e0,e2,tau,rr=pB
E1=sig/e0; E2=sig/e2; eta2=E2*tau; eta1=sig/rr if rr>0 else np.inf
sB=np.sqrt(np.diag(cB))
# validación predictiva: ajustar con t<=600 s, predecir el resto
cal=tb<=600; val=tb>600
pBc,_=curve_fit(burgers,tb[cal],eb[cal],p0=pB,sigma=sb[cal],bounds=([0,0,1,0],[0.05,0.05,1e4,1e-3]),maxfev=20000)
pFc,_=curve_fit(findley,tb[cal],eb[cal],p0=pF,sigma=sb[cal],bounds=([0,0,0.01],[0.05,0.05,1]),maxfev=20000)
errB=np.abs(burgers(tb[val],*pBc)-eb[val]); errF=np.abs(findley(tb[val],*pFc)-eb[val])
creep_tot=(burgers(tr.max(),*pB)-burgers(8,*pB))
out['prueba2']=dict(L0_px=L0,L0_sd_px=sL0,L0_mm_regla=L0/PXMM,px_por_mm=PXMM,t_carga_video_s=t_load,
  duracion_bajo_carga_s=float(tr.max()),sigma_MPa=sig/1e6,
  deformacion_8s_pct=float(burgers(8,*pB)*100),deformacion_final_pct=float(np.median(eps[ok&(tr>tr.max()-20)])*100),
  creep_8s_a_final_pct=float(creep_tot*100),
  burgers=dict(eps0=e0,eps2=e2,tau_s=tau,rate_1_s=rr,E1_GPa=E1/1e9,E2_GPa=E2/1e9,eta2_GPa_s=eta2/1e9,eta1_GPa_s=eta1/1e9,
     se=dict(eps0=sB[0],eps2=sB[1],tau=sB[2],rate=sB[3]),R2=r2(eb,burgers(tb,*pB)),RMSE_pct=float(np.sqrt(np.mean((eb-burgers(tb,*pB))**2))*100)),
  findley=dict(eps0=pF[0],A=pF[1],n=pF[2],R2=r2(eb,findley(tb,*pF)),RMSE_pct=float(np.sqrt(np.mean((eb-findley(tb,*pF))**2))*100)),
  prony2=dict(params=list(pP),R2=r2(eb,prony2(tb,*pP)),RMSE_pct=float(np.sqrt(np.mean((eb-prony2(tb,*pP))**2))*100)),
  validacion_600s=dict(burgers_err_max_pct=float(errB.max()*100),burgers_err_medio_pct=float(errB.mean()*100),
     findley_err_max_pct=float(errF.max()*100),findley_err_medio_pct=float(errF.mean()*100),
     burgers_pred_final_pct=float(burgers(tr.max(),*pBc)*100),findley_pred_final_pct=float(findley(tr.max(),*pFc)*100)))

# Figura P2: (a) lineal con datos crudos y modelo, (b) log, (c) residuos
tt=np.logspace(np.log10(8),np.log10(tr.max()),400)
fig,ax=plt.subplots(1,3,figsize=(14,4.2))
m=ok&(tr>=-3)
ax[0].plot(tr[m]/60,eps[m]*100,'.',ms=1.5,color=GRAY,label='Datos (5 fps)',rasterized=True)
ax[0].plot(tb/60,eb*100,'o',ms=4,color=BLUE,label='Mediana por intervalo')
ax[0].plot(tt/60,burgers(tt,*pB)*100,'-',lw=2,color=ORANGE,label=f'Burgers (R² = {out["prueba2"]["burgers"]["R2"]:.3f})')
ax[0].set(xlabel='Tiempo desde la carga (min)',ylabel='Deformación (%)',title='Prueba 2 · 1 kg · σ ≈ 22 MPa'); ax[0].legend(frameon=False,fontsize=8)
ax[0].set_ylim(-0.2, max(eps[ok&(tr>8)])*100+0.15)
ax[1].semilogx(tb,eb*100,'o',ms=4,color=BLUE,label='Datos')
ax[1].semilogx(tt,burgers(tt,*pB)*100,'-',lw=2,color=ORANGE,label='Burgers')
ax[1].semilogx(tt,findley(tt,*pF)*100,'--',lw=1.5,color=MUTED,label=f'Findley (R² = {out["prueba2"]["findley"]["R2"]:.3f})')
ax[1].set(xlabel='Tiempo desde la carga (s, escala log)',ylabel='Deformación (%)',title='Escala logarítmica'); ax[1].legend(frameon=False,fontsize=8)
ax[2].axhline(0,color=MUTED,lw=0.8)
ax[2].semilogx(tb,(eb-burgers(tb,*pB))*100,'o',ms=4,color=ORANGE,label='Burgers')
ax[2].semilogx(tb,(eb-findley(tb,*pF))*100,'s',ms=3.5,color=MUTED,mfc='none',label='Findley')
ax[2].set(xlabel='Tiempo desde la carga (s)',ylabel='Residuo (puntos de %)',title='Residuos del ajuste'); ax[2].legend(frameon=False,fontsize=8)
fig.tight_layout(); fig.savefig(f'{S}/prueba2_creep.png'); plt.close(fig)

# Figura validación
fig,ax=plt.subplots(figsize=(7,4.2))
ax.plot(tb/60,eb*100,'o',ms=4,color=BLUE,label='Datos')
tt2=np.linspace(8,tr.max(),400)
ax.plot(tt2/60,burgers(tt2,*pBc)*100,'-',lw=2,color=ORANGE,label='Burgers ajustado con t ≤ 10 min')
ax.axvspan(0,10,color='#f1f0ec',zorder=0); ax.text(0.3,ax.get_ylim()[1]*0.97 if ax.get_ylim()[1]>0 else 0,'calibración',color=MUTED,va='top',fontsize=8)
ax.set(xlabel='Tiempo desde la carga (min)',ylabel='Deformación (%)',title='Prueba 2 · validación predictiva (10 → 27 min)'); ax.legend(frameon=False,fontsize=8,loc='lower right')
fig.tight_layout(); fig.savefig(f'{S}/prueba2_validacion.png'); plt.close(fig)

# ---------------- PRUEBA 3 ----------------
d=load(f'{DATOS}/prueba3/seguimiento_crudo_px.csv'); t=d['t']; L=d['yl']-d['yu']; ok=~np.isnan(L)
L0=np.nanmean(L[t<=2.8]); PXMM3=2.23; t_load3=6.4; t_break=1223.7
sig3=4.5*g/A
tr=t-t_load3; eps=(L-L0)/L0
rows=[]
for s in np.arange(0,int(t_break-t_load3)):
    m=ok&(tr>=s)&(tr<s+1)
    if m.sum(): rows.append([s+0.5,np.median(L[m]),np.median(L[m])/PXMM3,np.median(eps[m])*100])
rows=np.array(rows)
with open(f'{S}/prueba3_deformacion_1s.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_s','L_px','L_mm','deformacion_%']); w.writerows(np.round(rows,5))
# velocidad de deformación (suavizada)
ts=rows[:,0]; es=rows[:,3]/100
k=15; esm=np.convolve(es,np.ones(k)/k,'same'); rate=np.gradient(esm,ts)
v=(ts>10)&(ts<ts.max()-10); ipk=np.argmax(np.where(v,rate,-1))
out['prueba3']=dict(nota='referencia = filamento ya cargado con 4 kg; t=0 al agregar 0.5 kg; escala en mm de la regla NO valida (regla fuera del plano del filamento)',sigma_4kg_MPa=4*g/A/1e6,L_ref_px=L0,L_ref_mm_regla=L0/PXMM3,px_por_mm=PXMM3,sigma_MPa=sig3/1e6,t_carga_video_s=t_load3,
  t_ruptura_video_s=t_break,tiempo_a_ruptura_s=t_break-t_load3,
  def_15s_pct=float(np.median(es[(ts>12)&(ts<18)])*100),def_60s_pct=float(np.median(es[(ts>55)&(ts<65)])*100),
  def_120s_pct=float(np.median(es[(ts>115)&(ts<125)])*100),def_600s_pct=float(np.median(es[(ts>590)&(ts<610)])*100),
  def_final_pct=float(np.median(es[ts>ts.max()-20])*100),
  vel_max_pct_s=float(rate[ipk]*100),t_vel_max_s=float(ts[ipk]),
  vel_final_pct_min=float(np.mean(rate[(ts>ts.max()-300)&(ts<ts.max()-10)])*100*60))
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
m=ok&(tr>-5)&(t<t_break)
ax[0].plot(tr[m]/60,eps[m]*100,'.',ms=1.5,color=BLUE,rasterized=True)
ax[0].axvline((t_break-t_load3)/60,color=MUTED,ls=':',lw=1); ax[0].text((t_break-t_load3)/60-0.3,1,'ruptura\n(en el nudo)',ha='right',color=MUTED,fontsize=8)
ax[0].set(xlabel='Tiempo desde que se agregaron 0.5 kg (min)',ylabel='Deformación desde el estado con 4 kg (%)',title='Prueba 3 · 4 kg → 4.5 kg (σ ≈ 89 → 100 MPa)')
ax[1].semilogx(ts[v],rate[v]*100*60,'-',lw=2,color=BLUE)
ax[1].set(xlabel='Tiempo desde la carga (s, escala log)',ylabel='Velocidad de deformación (%/min)',title='Velocidad de deformación')
fig.tight_layout(); fig.savefig(f'{S}/prueba3_creep.png'); plt.close(fig)
json.dump(out,open(f'{S}/resultados.json','w'),indent=2,default=float)
print(json.dumps(out,indent=1,default=float))
