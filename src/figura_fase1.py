"""Figura de la primera fase (2024-2025): alargamiento contra fuerza y rigidez aparente k.

Datos: probetas de 15 cm de tira cortada y de filamento extruido, con masas de 0.05 a 2 kg
(García Hernández y Ortíz de Zárate Granados, 2025). Las rectas se ajustan solo en el
intervalo lineal de cada material; k es el inverso de la pendiente.
"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=os.path.join(ROOT,'resultados','finales'); os.makedirs(FIG,exist_ok=True)
BLUE,GRAY,INK,MUTED='#2a78d6','#9b9a94','#1a1a19','#52514e'
plt.rcParams.update({'font.size':11,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e5e0','grid.linewidth':0.6,'axes.titlesize':12,'axes.titleweight':'bold'})

d=np.genfromtxt(os.path.join(ROOT,'datos','fase1','fase1_hooke.csv'),delimiter=',',names=True)
F=d['fuerza_N']; x1=d['alargamiento_tira_cortada_m']; x2=d['alargamiento_filamento_extruido_m']
lin1=(F<=8.83); lin2=(F>=0.49)&(F<=9.82)          # intervalos lineales
s1,b1=np.polyfit(F[lin1],x1[lin1],1); s2,b2=np.polyfit(F[lin2],x2[lin2],1)
r2=lambda x,y:np.corrcoef(x,y)[0,1]**2
k1,k2=1/s1,1/s2
print(f'Tira cortada: k = {k1:.0f} N/m, R2 = {r2(F[lin1],x1[lin1]):.3f}')
print(f'Filamento extruido: k = {k2:.0f} N/m, R2 = {r2(F[lin2],x2[lin2]):.2f}; razon = {k2/k1:.0f}x')

fig,ax=plt.subplots(figsize=(7,4.2))
ax.plot(F,x1*100,'o',ms=5.5,color=GRAY,label='Tira cortada con navaja')
ax.plot(F,x2*100,'s',ms=5,color=BLUE,label='Filamento extruido')
ff=np.linspace(0,F[lin1].max(),10); ax.plot(ff,(s1*ff+b1)*100,'-',color=GRAY,lw=1.8)
ff=np.linspace(F[lin2].min(),F[lin2].max(),10); ax.plot(ff,(s2*ff+b2)*100,'-',color=BLUE,lw=1.8)
ax.text(9.4,14.5,f'k ≈ {k1:.0f} N/m',color='#4a4945',fontsize=12,fontweight='bold')
ax.text(11,3.8,f'k ≈ {round(k2,-1):,.0f} N/m'.replace(',',' '),color=BLUE,fontsize=12,fontweight='bold')
ax.set(xlabel='Fuerza aplicada (N)',ylabel='Alargamiento (cm)',title='Primera fase: alargamiento contra fuerza (probetas de 15 cm)')
ax.set_ylim(-0.5,22); ax.legend(frameon=False,loc='lower right',bbox_to_anchor=(1,0.33))
fig.tight_layout(); fig.savefig(f'{FIG}/fig0_fase1_hooke.png',dpi=300); fig.savefig(f'{FIG}/fig0_fase1_hooke.pdf')
