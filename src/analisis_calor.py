"""Comparación filamento sin tratar vs calentado/estirado a 90 °C (4.5 kg), ajuste de Burgers y gráficas."""
import numpy as np, json, csv, os, sys
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS=os.path.join(ROOT,'datos'); OUT=os.path.join(ROOT,'resultados','calor'); os.makedirs(OUT,exist_ok=True)
BLUE,GRAY,INK,MUTED,ORANGE='#2a78d6','#9b9a94','#0b0b0b','#52514e','#eb6834'
plt.rcParams.update({'font.size':10,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
  'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e5e0','grid.linewidth':0.6,'figure.dpi':150})
g=9.81; RHO=1.335  # densidad supuesta (PET amorfo) para el área por masa lineal
def area_masa(m_g,L_cm): return (m_g/L_cm)/RHO*1e-4  # m^2  (g/cm / g/cm3 = cm2)
A_sin=area_masa(0.0993,18.9); A_cal=area_masa(0.1260,19.5)
A_sin_v=np.pi*(0.75e-3/2)**2; A_cal_v=np.pi*(0.62e-3/2)**2
def burgers(t,e0,e2,tau,r): return e0+e2*(1-np.exp(-t/tau))+r*t
def prony(t,e0,a1,t1,a2,t2,r): return e0+a1*(1-np.exp(-t/t1))+a2*(1-np.exp(-t/t2))+r*t
PB=([0,0,1,0,30,0],[0.2,0.2,500,0.2,1e5,1e-2])
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
def med_bins(t,y,edges):
    tb=[];yb=[]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(t>=a)&(t<b)&~np.isnan(y)
        if m.sum()>=2: tb.append(np.median(t[m])); yb.append(np.median(y[m]))
    return np.array(tb),np.array(yb)
# ---- calentado (prueba 3 con calor) ----
A=np.genfromtxt(f'{DATOS}/prueba3_calor/seguimiento_1kg_px.csv',delimiter=',',names=True); B=np.genfromtxt(f'{DATOS}/prueba3_calor/seguimiento_px.csv',delimiter=',',names=True)
L1=np.nanmean((A['yl']-A['yu'])[(A['t']>=0.5)&(A['t']<=1.3)])   # con 1 kg
TL=13.0  # fin de la carga a 4.5 kg (s de video)
tc=B['t']-TL; Lc=B['yl']-B['yu']; ec=(Lc-L1)/L1
# ---- sin tratar (prueba 3 original) ----
d=np.genfromtxt(f'{DATOS}/prueba3/seguimiento_crudo_px.csv',delimiter=',',names=True); L3=d['yl']-d['yu']
L4kg=np.nanmean(L3[d['t']<=2.8]); t3=d['t']-6.4; e3=(L3-L4kg)/L4kg; brk3=1223.7-6.4
ok3=(t3>=0)&(t3<brk3)&~np.isnan(e3)
# deformación de creep desde 10 s después de completar la carga
def ref10(t,e): m=(t>=9)&(t<11); return np.nanmedian(e[m])
dc=ec-ref10(tc,ec); d3=e3-ref10(t3,e3)
# ---- ajuste Burgers del calentado (t>=10 s tras la carga) ----
edges=np.r_[10,np.logspace(np.log10(12),np.log10(tc.max()),45)]
tb,eb=med_bins(tc,ec,edges)
p,cov=curve_fit(burgers,tb,eb,p0=[0.025,0.02,60,1e-5],bounds=([0,0,1,0],[0.2,0.2,5000,1e-2]),maxfev=50000)
sig_c=4.5*g/A_cal; sig1=1.0*g/A_cal; dsig=sig_c-sig1
# los parámetros se refieren al incremento de esfuerzo 1 kg -> 4.5 kg (la referencia es el estado con 1 kg)
E1=dsig/p[0]; E2=dsig/p[1]; eta2=E2*p[2]; eta1=dsig/p[3] if p[3]>0 else np.inf
cal=tb<=300; pc,_=curve_fit(burgers,tb[cal],eb[cal],p0=p,bounds=([0,0,1,0],[0.2,0.2,5000,1e-2]),maxfev=50000)
val=tb>300; err=np.abs(burgers(tb[val],*pc)-eb[val])
q,_=curve_fit(prony,tb,eb,p0=[0.02,0.01,20,0.02,300,1e-6],bounds=PB,maxfev=200000)
qc,_=curve_fit(prony,tb[cal],eb[cal],p0=q,bounds=PB,maxfev=200000)
errq=np.abs(prony(tb[val],*qc)-eb[val])
res=dict(
 area_mm2=dict(sin_tratar_masa=A_sin*1e6,calentado_masa=A_cal*1e6,sin_tratar_vernier=A_sin_v*1e6,calentado_vernier=A_cal_v*1e6,rho_supuesta=RHO),
 diametro_equivalente_mm=dict(sin_tratar=2*np.sqrt(A_sin/np.pi)*1e3,calentado=2*np.sqrt(A_cal/np.pi)*1e3),
 esfuerzo_MPa=dict(sin_tratar_4kg=4*g/A_sin/1e6,sin_tratar_4_5kg=4.5*g/A_sin/1e6,calentado_1kg=sig1/1e6,calentado_4_5kg=sig_c/1e6,
                  calentado_4_5kg_zona_0_62mm=4.5*g/A_cal_v/1e6,prueba2_1kg=1*g/A_sin/1e6),
 relacion_estiramiento_por_largo=19.5/19.0,
 calentado=dict(L_1kg_px=L1,def_inmediata_1a4_5kg_pct=float(np.nanmedian(ec[(tc>=0)&(tc<2)])*100),
   def_10s_pct=float(ref10(tc,ec)*100),def_final_pct=float(np.nanmedian(ec[tc>tc.max()-20])*100),
   creep_10s_a_final_pct=float(np.nanmedian(dc[tc>tc.max()-20])*100),duracion_s=float(tc.max()),se_rompio=False,
   burgers=dict(eps0=p[0],eps2=p[1],tau_s=p[2],rate_1_s=p[3],E1_GPa=E1/1e9,E2_GPa=E2/1e9,eta2_Pa_s=eta2,eta1_Pa_s=eta1,
                R2=r2(eb,burgers(tb,*p)),RMSE_pct=float(np.sqrt(np.mean((eb-burgers(tb,*p))**2))*100)),
   validacion_5min_burgers=dict(pred_final_pct=float(burgers(tb[-1],*pc)*100),medido_final_pct=float(eb[-1]*100),
                        err_max_pct=float(err.max()*100),R2_validacion=r2(eb[val],burgers(tb[val],*pc))),
   prony=dict(eps0=q[0],a1=q[1],tau1_s=q[2],a2=q[3],tau2_s=q[4],rate_1_s=q[5],E0_GPa=dsig/q[0]/1e9,E1_GPa=dsig/q[1]/1e9,E2_GPa=dsig/q[3]/1e9,
              eta_flujo_Pa_s=dsig/q[5] if q[5]>0 else None,R2=r2(eb,prony(tb,*q)),RMSE_pct=float(np.sqrt(np.mean((eb-prony(tb,*q))**2))*100)),
   validacion_5min_prony=dict(pred_final_pct=float(prony(tb[-1],*qc)*100),medido_final_pct=float(eb[-1]*100),
                        err_max_pct=float(errq.max()*100),R2_validacion=r2(eb[val],prony(tb[val],*qc)))),
 sin_tratar=dict(creep_10s_a_ruptura_pct=float(np.nanmedian(d3[ok3&(t3>brk3-20)])*100),creep_10s_a_16min_pct=float(np.nanmedian(d3[ok3&(t3>tc.max()-10)&(t3<tc.max()+10)])*100),
   t_ruptura_s=brk3),
 prueba4=dict(carga_ruptura_kg=4.5,t_ruptura_video_s=99.0,lugar='nudo superior'))
json.dump(res,open(f'{OUT}/resultados_calor.json','w'),indent=2,default=float)
# CSV 1 s
with open(f'{OUT}/calentado_deformacion_1s.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_4_5kg_s','L_px','deformacion_desde_1kg_%'])
    for s in range(0,int(tc.max())):
        m=(tc>=s)&(tc<s+1)
        if m.any(): w.writerow([s+0.5,round(float(np.nanmedian(Lc[m])),4),round(float(np.nanmedian(ec[m]))*100,5)])
# ---- Figura 1: comparación de creep ----
fig,ax=plt.subplots(figsize=(7.6,4.4))
m3=ok3&(t3>=10); mc=(tc>=10)
ax.plot(t3[m3]/60,d3[m3]*100,'-',lw=1.6,color=GRAY,label=f'Sin tratar · ≈{4.5*g/A_sin/1e6:.0f} MPa')
ax.plot(tc[mc]/60,dc[mc]*100,'-',lw=2,color=BLUE,label=f'Calentado a 90 °C · ≈{sig_c/1e6:.0f} MPa')
ax.axvline(brk3/60,color=MUTED,ls=':',lw=1); ax.text(brk3/60-0.2,9.5,'ruptura\n(nudo)',ha='right',color=MUTED,fontsize=8)
ax.set(xlabel='Tiempo desde que se completaron los 4.5 kg (min)',ylabel='Deformación por creep (%)',
       title=f'Con 4.5 kg: sin tratar {res["sin_tratar"]["creep_10s_a_16min_pct"]:.1f} % vs calentado {res["calentado"]["creep_10s_a_final_pct"]:.1f} % a los 16 min')
ax.legend(frameon=False,fontsize=9,loc='center',bbox_to_anchor=(0.62,0.45)); fig.tight_layout(); fig.savefig(f'{OUT}/comparacion_calor_vs_sin_tratar.png'); plt.close(fig)
# ---- Figura 2: calentado + Burgers + validación ----
tt=np.linspace(10,tc.max(),400)
fig,ax=plt.subplots(1,2,figsize=(12,4.2))
ax[0].plot(tc[tc>=-2]/60,ec[tc>=-2]*100,'.',ms=1.2,color=GRAY,rasterized=True,label='Datos (5 fps)')
ax[0].plot(tt/60,prony(tt,*q)*100,'-',lw=2,color=BLUE,label=f'Prony, 2 términos (R² = {res["calentado"]["prony"]["R2"]:.3f})')
ax[0].plot(tt/60,burgers(tt,*p)*100,'--',lw=1.3,color=MUTED,label=f'Burgers (R² = {res["calentado"]["burgers"]["R2"]:.3f})')
ax[0].set(xlabel='Tiempo desde la carga de 4.5 kg (min)',ylabel='Deformación desde el estado con 1 kg (%)',title='Filamento calentado, 4.5 kg · ajuste del modelo'); ax[0].legend(frameon=False,fontsize=8,loc='lower right')
ax[1].plot(tb/60,eb*100,'o',ms=4,color=GRAY,label='Datos')
ax[1].plot(tt/60,prony(tt,*qc)*100,'-',lw=2,color=BLUE,label='Prony ajustado con los primeros 5 min')
ax[1].plot(tt/60,burgers(tt,*pc)*100,'--',lw=1.3,color=MUTED,label='Burgers ajustado con los primeros 5 min')
ax[1].axvspan(0,5,color='#f1f0ec',zorder=0); ax[1].text(0.2,eb.max()*100,'calibración',color=MUTED,fontsize=8,va='top')
ax[1].set(xlabel='Tiempo desde la carga de 4.5 kg (min)',ylabel='Deformación (%)',title=f'Validación: Prony predice {res["calentado"]["validacion_5min_prony"]["pred_final_pct"]:.2f} % vs {res["calentado"]["validacion_5min_prony"]["medido_final_pct"]:.2f} % medido (R² = {res["calentado"]["validacion_5min_prony"]["R2_validacion"]:.2f})'); ax[1].legend(frameon=False,fontsize=8,loc='lower right'); ax[1].set_ylim(eb.min()*100-0.3, eb.max()*100+0.6)
fig.tight_layout(); fig.savefig(f'{OUT}/calentado_modelo.png'); plt.close(fig)
print(json.dumps(res,indent=1,default=float))
