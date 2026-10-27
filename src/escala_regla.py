# Escala (px/mm) de la regla a lo largo del video, para corregir cambios de enfoque/zoom del teléfono
import cv2, numpy as np, sys
video,out,x0,x1,y0,y1,mode=sys.argv[1],sys.argv[2],*map(int,sys.argv[3:7]),sys.argv[7]
cap=cv2.VideoCapture(video); fps=cap.get(cv2.CAP_PROP_FPS); n=int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
def period(g):
    ys=np.arange(y0,y1); gf=g.astype(np.float32)
    if mode=='edge':   # muestrear junto al borde izquierdo de la regla
        edge=[]
        for y in ys:
            r=np.convolve(g[y,x0:x1].astype(float),np.ones(3)/3,'same'); idx=np.where(r>180)[0]
            edge.append(x0+idx[0] if len(idx) else np.nan)
        edge=np.array(edge,float); ok=~np.isnan(edge); p=np.polyfit(ys[ok],edge[ok],1)
        xs=[np.polyval(p,y)+5 for y in ys]
    else:              # columna fija
        xs=[x0]*len(ys)
    prof=np.array([cv2.getRectSubPix(gf,(3,1),(float(x),float(y)))[0].mean() for x,y in zip(xs,ys)])
    prof=prof-np.convolve(prof,np.ones(15)/15,'same'); N=len(prof)
    F=np.abs(np.fft.rfft(prof*np.hanning(N),N*64)); fr=np.fft.rfftfreq(N*64); m=(fr>1/3.2)&(fr<1/1.9); k=np.argmax(F*m)
    return 1/fr[k], F[k]/np.median(F[m])
res=[]
for fr in range(0,n,25):
    cap.set(cv2.CAP_PROP_POS_FRAMES,fr); ok,f=cap.read()
    if not ok: break
    p,q=period(cv2.cvtColor(f,cv2.COLOR_BGR2GRAY)); res.append([fr/fps,p,q])
np.savetxt(out,np.array(res),delimiter=',',header='t,px_por_mm,calidad',comments='',fmt='%.5f')
