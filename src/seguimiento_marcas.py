# Seguimiento de banderas blancas: centro vertical = punto medio entre bordes a media altura
# Uso: python seguimiento_marcas.py VIDEO SALIDA.csv x_sup y_sup x_inf y_inf [medio_ancho] [medio_alto]
# Variables opcionales: F0/F1 (cuadro inicial/final), JMAX (salto máximo en px entre cuadros, default 10),
#                       HMAX (alto máximo de la bandera en px, default 25)
import cv2, numpy as np, sys
video, out = sys.argv[1], sys.argv[2]
x_u, y_u, x_l, y_l = map(float, sys.argv[3:7])
HW = int(sys.argv[7]) if len(sys.argv) > 7 else 30
HH = int(sys.argv[8]) if len(sys.argv) > 8 else 20
cap = cv2.VideoCapture(video); fps = cap.get(cv2.CAP_PROP_FPS)
import os
F0 = int(os.environ.get('F0', '0')); F1 = int(os.environ.get('F1', '1000000000'))
cap.set(cv2.CAP_PROP_POS_FRAMES, F0)
pts = {'u': [x_u, y_u], 'l': [x_l, y_l]}
def measure(gray, p):
    x0, y0 = int(round(p[0])), int(round(p[1]))
    g = gray[y0-HH:y0+HH+1, x0-HW:x0+HW+1].astype(float)
    bg = np.median(g)
    prof = np.sort(g, axis=1)[:, -5:].mean(1)        # brillo máx. robusto por fila
    c = HH                                            # buscar el máximo cerca del centro previsto
    lo, hi = max(0, c-10), min(len(prof), c+11)
    k = lo + int(np.argmax(prof[lo:hi])); m = prof[k]
    if m - bg < 25: return None
    h = bg + 0.5*(m - bg)
    a = k
    while a > 0 and prof[a-1] >= h: a -= 1
    b = k
    while b < len(prof)-1 and prof[b+1] >= h: b += 1
    if a == 0 or b == len(prof)-1: return None
    ya = a - (prof[a]-h)/(prof[a]-prof[a-1])
    yb = b + (prof[b]-h)/(prof[b]-prof[b+1])
    rows = g[a:b+1]; mask = rows > h
    cx = np.where(mask)[1].mean() + x0 - HW
    return cx, 0.5*(ya+yb) + y0 - HH, yb-ya
rows = []; i = F0
while i < F1:
    ok, f = cap.read()
    if not ok: break
    gray = cv2.cvtColor(f, cv2.COLOR_BGR2GRAY)
    r = [i, i/fps]
    for k in ('u', 'l'):
        m = measure(gray, pts[k])
        if m and abs(m[1]-pts[k][1]) < float(os.environ.get("JMAX","10")) and 4 < m[2] < float(os.environ.get("HMAX","25")):
            pts[k] = [m[0], m[1]]; r += list(m)
        else: r += [np.nan]*3
    rows.append(r); i += 1
np.savetxt(out, np.array(rows), delimiter=',', header='frame,t,xu,yu,hu,xl,yl,hl', comments='', fmt='%.4f')
print('frames', i)
