"""Figuras finales para póster y reporte (EcoPlastic). Genera PNG 300 dpi y PDF."""
import numpy as np, json, csv, os, cv2
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=TAB=DAT=os.path.join(ROOT,'resultados','finales'); os.makedirs(FIG,exist_ok=True)
VIDEO=os.environ.get('VIDEO_CALOR','')  # ruta opcional al video de la prueba 3 con calor para la figura 1
BLUE,GRAY,INK,MUTED,LIGHT='#2a78d6','#9b9a94','#1a1a19','#52514e','#f1f0ec'
plt.rcParams.update({'font.size':11,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e5e0','grid.linewidth':0.6,'axes.titlesize':12,'axes.titleweight':'bold'})
def save(fig,name):
    fig.tight_layout(); fig.savefig(f'{FIG}/{name}.png',dpi=300); fig.savefig(f'{FIG}/{name}.pdf'); plt.close(fig)
g=9.81; RHO=1.335
A_sin=(0.0993/18.9)/RHO*1e-4; A_cal=(0.1260/19.5)/RHO*1e-4
def burgers(t,e0,e2,tau,r): return e0+e2*(1-np.exp(-t/tau))+r*t
def prony(t,e0,a1,t1,a2,t2,r): return e0+a1*(1-np.exp(-t/t1))+a2*(1-np.exp(-t/t2))+r*t
PB=([0,0,1,0,30,0],[0.2,0.2,500,0.2,1e5,1e-2]); BB=([0,0,1,0],[0.2,0.2,1e4,1e-2])
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
def logbins(t,y,t0,t1,n,mask=None):
    e=np.logspace(np.log10(t0),np.log10(t1),n+1); tb=[];yb=[]
    for a,b in zip(e[:-1],e[1:]):
        m=(t>=a)&(t<b)&~np.isnan(y)
        if mask is not None: m&=mask
        if m.sum()>=2: tb.append(np.median(t[m])); yb.append(np.median(y[m]))
    return np.array(tb),np.array(yb)
def per_second(t,L,e,tmax,mask=None):
    rows=[]
    for s in range(0,int(tmax)):
        m=(t>=s)&(t<s+1)&~np.isnan(e)
        if mask is not None: m&=mask
        if m.any(): rows.append([s+0.5,float(np.median(L[m])),float(np.median(e[m]))*100])
    return rows
R={}
# ===== Prueba 2 =====
d=np.genfromtxt(os.path.join(ROOT,'datos','prueba2','seguimiento_crudo_px.csv'),delimiter=',',names=True); t=d['t']; L=d['yl']-d['yu']
ok2=~np.isnan(L)&(np.abs(d['hu']-11.57)<0.8)&(np.abs(d['hl']-12.75)<0.8)
L0=L[(t>=2)&(t<=4.2)].mean(); t2=t-4.4; e2=(L-L0)/L0; sig2=1*g/A_sin
tb,eb=logbins(t2,e2,8,t2.max(),40,ok2)
p2,_=curve_fit(burgers,tb,eb,p0=[0.004,0.002,60,1e-6],bounds=BB,maxfev=50000)
R['prueba2']=dict(carga_kg=1,esfuerzo_MPa=sig2/1e6,L0_px=L0,L0_mm=L0/2.856,duracion_min=t2.max()/60,
  def_8s_pct=float(burgers(8,*p2)*100),def_final_pct=float(np.median(e2[ok2&(t2>t2.max()-20)])*100),
  burgers=dict(E1_GPa=sig2/p2[0]/1e9,E2_GPa=sig2/p2[1]/1e9,tau_s=p2[2],eta2_Pa_s=sig2/p2[1]*p2[2],eta1_Pa_s=sig2/p2[3],R2=r2(eb,burgers(tb,*p2))))
# ruido / precisión: ventana tranquila 1450-1550 s
mq=ok2&(t2>=1450)&(t2<1550); res=L[mq]-np.polyval(np.polyfit(t2[mq],L[mq],1),t2[mq])
R['precision']=dict(ruido_px_desv=float(res.std()),ruido_def_pct=float(res.std()/L0*100),ventana_s='1450-1550 (prueba 2)',
  mm_por_px=1/2.856, ruido_um=float(res.std()/2.856*1000))
fig,ax=plt.subplots(figsize=(7,3.6))
ax.plot(t2[mq]-1450,res/L0*100,'.',ms=3,color=BLUE)
s=R['precision']['ruido_def_pct']; ax.axhspan(-s,s,color=LIGHT,zorder=0)
ax.axhline(0,color=MUTED,lw=0.8)
ax.set(xlabel='Tiempo (s)',ylabel='Variación de la deformación (%)',title=f'Precisión de la medición: ±{s:.3f} % de deformación (1 desviación estándar)')
ax.set_ylim(-0.05,0.05); save(fig,'fig2_precision_medicion')
# figura prueba 2
fig,ax=plt.subplots(figsize=(7,4.2))
ax.plot(t2[ok2&(t2>=-3)]/60,e2[ok2&(t2>=-3)]*100,'.',ms=1.2,color=GRAY,rasterized=True,label='Datos (5 cuadros/s)')
tt=np.linspace(8,t2.max(),400); ax.plot(tt/60,burgers(tt,*p2)*100,'-',lw=2.2,color=BLUE,label=f'Modelo de Burgers (R² = {R["prueba2"]["burgers"]["R2"]:.2f})')
ax.set_ylim(-0.05,1.0); ax.set(xlabel='Tiempo desde la carga (min)',ylabel='Deformación (%)',title=f'Filamento sin tratar, 1 kg (≈{sig2/1e6:.0f} MPa): {R["prueba2"]["def_final_pct"]:.2f} % en 27 min')
ax.legend(frameon=False,loc='lower right'); save(fig,'fig3_creep_sin_tratar_1kg')
rows=per_second(t2,L,e2,t2.max(),ok2)
with open(f'{DAT}/prueba2_sin_tratar_1kg_deformacion_por_segundo.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_s','L_entre_marcas_px','deformacion_%']); w.writerows([[r[0],round(r[1],4),round(r[2],5)] for r in rows])
# ===== Prueba 3 sin tratar =====
d=np.genfromtxt(os.path.join(ROOT,'datos','prueba3','seguimiento_crudo_px.csv'),delimiter=',',names=True); L3=d['yl']-d['yu']; ok3=~np.isnan(L3)
L4=np.nanmean(L3[d['t']<=2.8]); t3=d['t']-6.4; e3=(L3-L4)/L4; brk=1223.7-6.4
rows3=np.array(per_second(t3,L3,e3,brk,ok3)); ts=rows3[:,0]; es=rows3[:,2]
k=15; esm=np.convolve(es,np.ones(k)/k,'same'); rate=np.gradient(esm,ts)*60
v=(ts>10)&(ts<ts.max()-10); ipk=np.argmax(np.where(v,rate,-1))
R['prueba3_sin_tratar']=dict(carga_kg='4 -> 4.5',esfuerzo_MPa=[4*g/A_sin/1e6,4.5*g/A_sin/1e6],def_15s_pct=float(np.median(es[(ts>12)&(ts<18)])),
  def_2min_pct=float(np.median(es[(ts>115)&(ts<125)])),def_16min_pct=float(np.median(es[(ts>970)&(ts<990)])),def_ruptura_pct=float(np.median(es[ts>ts.max()-20])),
  vel_max_pct_min=float(rate[ipk]),t_vel_max_s=float(ts[ipk]),t_ruptura_min=brk/60,lugar_ruptura='nudo superior')
fig,ax=plt.subplots(1,2,figsize=(11,4))
ax[0].plot(ts/60,es,'-',lw=1.8,color=GRAY); ax[0].axvline(brk/60,color=MUTED,ls=':',lw=1); ax[0].text(brk/60-0.3,4,'ruptura\n(en el nudo)',ha='right',color=MUTED,fontsize=9)
ax[0].set(xlabel='Tiempo desde que se agregaron 0.5 kg (min)',ylabel='Deformación desde el estado con 4 kg (%)',title='Sin tratar, 4 → 4.5 kg: cedencia retardada')
ax[1].semilogx(ts[v],rate[v],'-',lw=1.8,color=GRAY); ax[1].plot(ts[ipk],rate[ipk],'o',color=BLUE,ms=7)
ax[1].annotate(f'máximo: {rate[ipk]:.0f} %/min a los {ts[ipk]:.0f} s',(ts[ipk],rate[ipk]),xytext=(ts[ipk]*2.2,rate[ipk]*0.85),fontsize=9,color=INK,arrowprops=dict(arrowstyle='-',color=MUTED,lw=0.8))
ax[1].set(xlabel='Tiempo (s, escala logarítmica)',ylabel='Velocidad de deformación (%/min)',title='Velocidad de deformación')
save(fig,'fig4_cedencia_sin_tratar_4_5kg')
with open(f'{DAT}/prueba3_sin_tratar_4_5kg_deformacion_por_segundo.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_s','L_entre_marcas_px','deformacion_desde_4kg_%']); w.writerows([[r[0],round(r[1],4),round(r[2],5)] for r in rows3])
# ===== Prueba 3 con calor =====
A_=np.genfromtxt(os.path.join(ROOT,'datos','prueba3_calor','seguimiento_1kg_px.csv'),delimiter=',',names=True); B_=np.genfromtxt(os.path.join(ROOT,'datos','prueba3_calor','seguimiento_px.csv'),delimiter=',',names=True)
L1=np.nanmean((A_['yl']-A_['yu'])[(A_['t']>=0.5)&(A_['t']<=1.3)]); tc=B_['t']-13.0; Lc=B_['yl']-B_['yu']; ec=(Lc-L1)/L1
def ref10(t,e): m=(t>=9)&(t<11); return np.nanmedian(e[m])
dc=ec-ref10(tc,ec); okc=~np.isnan(ec)
d3=e3-ref10(t3,e3); m3=ok3&(t3>=10)&(t3<brk)
sigc=4.5*g/A_cal
# modelos
edges=np.r_[10,np.logspace(np.log10(12),np.log10(tc.max()),45)]; tb=[];eb=[]
for a,b in zip(edges[:-1],edges[1:]):
    m=(tc>=a)&(tc<b)&okc
    if m.sum()>=2: tb.append(np.median(tc[m])); eb.append(np.median(ec[m]))
tb=np.array(tb); eb=np.array(eb)
q,_=curve_fit(prony,tb,eb,p0=[0.02,0.01,20,0.02,300,1e-6],bounds=PB,maxfev=200000)
pb,_=curve_fit(burgers,tb,eb,p0=[0.02,0.02,30,1e-5],bounds=BB,maxfev=200000)
cal=tb<=300; val=~cal
qc,_=curve_fit(prony,tb[cal],eb[cal],p0=q,bounds=PB,maxfev=200000)
pc,_=curve_fit(burgers,tb[cal],eb[cal],p0=pb,bounds=BB,maxfev=200000)
dsig=sigc-1*g/A_cal
R['prueba3_calentado']=dict(carga_kg='1 -> 4.5',esfuerzo_MPa=[1*g/A_cal/1e6,sigc/1e6],temperatura_agua_C=90,largo_cm='19 -> 19.5',lambda_largo=19.5/19,
  def_inmediata_pct=float(ref10(tc,ec)*100),def_final_pct=float(np.median(ec[okc&(tc>tc.max()-20)])*100),
  creep_desde_10s_16min_pct=float(np.median(dc[okc&(tc>tc.max()-20)])*100),duracion_min=tc.max()/60,cedio=False,se_rompio=False,
  prony=dict(a1_pct=q[1]*100,tau1_s=q[2],a2_pct=q[3]*100,tau2_s=q[4],flujo_pct_por_min=q[5]*6000,R2=r2(eb,prony(tb,*q))),
  burgers=dict(R2=r2(eb,burgers(tb,*pb))),
  validacion=dict(calibracion='primeros 5 min',prediccion='5 a 16 min',prony_pred_final_pct=float(prony(tb[-1],*qc)*100),
    burgers_pred_final_pct=float(burgers(tb[-1],*pc)*100),medido_final_pct=float(eb[-1]*100),
    prony_R2_validacion=r2(eb[val],prony(tb[val],*qc)),burgers_R2_validacion=r2(eb[val],burgers(tb[val],*pc)),
    prony_error_max_puntos_pct=float(np.abs(prony(tb[val],*qc)-eb[val]).max()*100)))
R['sin_tratar_creep_desde_10s_16min_pct']=float(np.nanmedian(d3[ok3&(t3>tc.max()-10)&(t3<tc.max()+10)])*100)
R['prueba4_calentado']=dict(carga='progresiva hasta 4.5 kg',ruptura='≈10 s después de completar 4.5 kg, en el nudo superior',esfuerzo_MPa=sigc/1e6)
R['areas']=dict(metodo='A = (masa/largo)/densidad, densidad supuesta 1.335 g/cm3',sin_tratar_mm2=A_sin*1e6,calentado_mm2=A_cal*1e6,
  d_equiv_sin_tratar_mm=2*np.sqrt(A_sin/np.pi)*1e3,d_equiv_calentado_mm=2*np.sqrt(A_cal/np.pi)*1e3,
  vernier_reportado_mm=dict(sin_tratar=0.75,calentado_antes=0.80,calentado_despues=0.62))
# fig comparación
fig,ax=plt.subplots(figsize=(8,4.6))
ax.plot(t3[m3]/60,d3[m3]*100,'-',lw=1.8,color=GRAY,label=f'Sin tratar (≈{4.5*g/A_sin/1e6:.0f} MPa)')
ax.plot(tc[okc&(tc>=10)]/60,dc[okc&(tc>=10)]*100,'-',lw=2.2,color=BLUE,label=f'Calentado a 90 °C (≈{sigc/1e6:.0f} MPa)')
ax.axvline(brk/60,color=MUTED,ls=':',lw=1); ax.text(brk/60-0.2,9.5,'ruptura\n(nudo)',ha='right',color=MUTED,fontsize=9)
ax.set(xlabel='Tiempo con 4.5 kg (min)',ylabel='Deformación por creep (%)',
  title=f'Con 4.5 kg: sin tratar {R["sin_tratar_creep_desde_10s_16min_pct"]:.0f} % vs calentado {R["prueba3_calentado"]["creep_desde_10s_16min_pct"]:.1f} % a los 16 min')
ax.legend(frameon=False,loc='center',bbox_to_anchor=(0.62,0.45)); save(fig,'fig5_comparacion_calentado_vs_sin_tratar')
# fig validación
tt=np.linspace(10,tc.max(),400)
fig,ax=plt.subplots(figsize=(8,4.6))
ax.axvspan(0,5,color=LIGHT,zorder=0); ax.text(0.15,6.45,'calibración\n(5 min)',color=MUTED,fontsize=9,va='top')
ax.text(5.25,2.1,'predicción (5 a 16 min) →',color=MUTED,fontsize=9,va='bottom')
ax.plot(tb/60,eb*100,'o',ms=4.5,color=GRAY,label='Medido',zorder=3)
ax.plot(tt/60,prony(tt,*qc)*100,'-',lw=2.4,color=BLUE,label=f'Prony (R² validación = {R["prueba3_calentado"]["validacion"]["prony_R2_validacion"]:.2f})')
ax.plot(tt/60,burgers(tt,*pc)*100,'--',lw=1.4,color=MUTED,label='Burgers')
ax.set_ylim(2,6.6); ax.set(xlabel='Tiempo con 4.5 kg (min)',ylabel='Deformación (%)',
  title=f'El modelo predice: {R["prueba3_calentado"]["validacion"]["prony_pred_final_pct"]:.2f} % estimado vs {R["prueba3_calentado"]["validacion"]["medido_final_pct"]:.2f} % medido')
ax.legend(frameon=False,loc='lower right'); save(fig,'fig6_validacion_modelo')
rowsc=per_second(tc,Lc,ec,tc.max())
with open(f'{DAT}/prueba3_calentado_4_5kg_deformacion_por_segundo.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t_desde_carga_4_5kg_s','L_entre_marcas_px','deformacion_desde_1kg_%']); w.writerows([[r[0],round(r[1],4),round(r[2],5)] for r in rowsc])
# ===== fig 1: montaje anotado =====
cap=cv2.VideoCapture(VIDEO); cap.set(cv2.CAP_PROP_POS_FRAMES,1500); ok,fr=cap.read()
if not ok:
    json.dump(R,open(f'{TAB}/resultados_completos.json','w'),indent=2,default=float,ensure_ascii=False); raise SystemExit('Figura 1 omitida: define VIDEO_CALOR con la ruta del video')
 fr=cv2.resize(fr,(608,1080),interpolation=cv2.INTER_AREA)
j=np.argmin(np.abs(B_['t']-300)); sx=608/1080
xu,yu,xl,yl=B_['xu'][j]*sx,B_['yu'][j],B_['xl'][j]*sx,B_['yl'][j]
fig,ax=plt.subplots(figsize=(4.2,7.4)); ax.imshow(cv2.cvtColor(fr,cv2.COLOR_BGR2RGB)); ax.axis('off'); ax.grid(False)
for (x,y,txt) in [(xu,yu,'marca superior'),(xl,yl,'marca inferior')]:
    ax.add_patch(plt.Circle((x,y),11,fill=False,ec=BLUE,lw=2)); ax.annotate(txt,(x+12,y),xytext=(x+95,y-20 if 'sup' in txt else y+25),color='white',fontsize=9,fontweight='bold',
       arrowprops=dict(arrowstyle='-',color='white',lw=1),bbox=dict(boxstyle='round,pad=0.25',fc=BLUE,ec='none'))
ax.annotate('',xy=(xu-25,yl),xytext=(xu-25,yu),arrowprops=dict(arrowstyle='<->',color=BLUE,lw=1.8)); ax.text(xu-32,(yu+yl)/2,'L(t)',color='white',fontsize=10,fontweight='bold',ha='right',va='center',bbox=dict(boxstyle='round,pad=0.2',fc=BLUE,ec='none'))
ax.set_title('Banco de creep: filamento, dos marcas y masa colgada',fontsize=10)
save(fig,'fig1_banco_creep')
json.dump(R,open(f'{TAB}/resultados_completos.json','w'),indent=2,default=float,ensure_ascii=False)
print(json.dumps(R,indent=1,default=float,ensure_ascii=False))
