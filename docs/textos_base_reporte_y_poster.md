# EcoPlastic: textos base para el reporte y el póster

Borrador para adaptar. Todos los números salen de `02_tablas/` y corresponden a las figuras de `01_figuras/`.

---

## Título

**Medición y predicción de la fluencia lenta (creep) en filamento de PET reciclado con instrumentación de bajo costo**

Alternativa más corta para el póster: *Creep del PET reciclado medido con un celular y predicho con Python*

---

## Resumen (≈ 200 palabras)

El PET reciclado suele perder resistencia tras el reprocesado, y no hay métodos accesibles para evaluar cómo se deforma con el tiempo bajo carga sostenida (creep). En este trabajo desarrollamos un banco de creep de bajo costo: un filamento de PET reciclado con dos marcas, una masa colgada y un teléfono celular que graba el ensayo. Un programa en Python con OpenCV localiza las marcas en cada cuadro y calcula la deformación con una precisión de ±0.007 %. Con este método caracterizamos el filamento sin tratar: con 1 kg (≈ 25 MPa) se deformó 0.91 % en 27 min, y con 4.5 kg (≈ 112 MPa) presentó cedencia retardada, con 14 % de deformación y ruptura a los 20 min. Después evaluamos un tratamiento en agua a 90 °C con estirado. Con la misma masa de 4.5 kg, el filamento tratado se deformó por creep 3.8 % en 16 min, contra 13 % del filamento sin tratar, y no cedió. Finalmente, modelamos el creep con los modelos de Burgers y de series de Prony. Una serie de Prony de dos términos, ajustada solo con los primeros 5 minutos, predijo los siguientes 11 minutos con R² = 0.98. Concluimos que es posible caracterizar y predecir el creep del PET reciclado sin equipo industrial.

**Palabras clave:** PET reciclado, fluencia lenta (creep), viscoelasticidad, serie de Prony, extensometría por video, economía circular.

---

## Planteamiento del problema

El PET reciclado se degrada con el reprocesado y se vuelve más frágil, lo que limita sus usos a productos de bajo valor. Para darle usos estructurales hay que saber cómo se deforma bajo carga sostenida, pero los ensayos de creep normalmente requieren máquinas universales y extensómetros costosos.

## Objetivo

Desarrollar y validar un método de bajo costo (celular y Python) para medir el creep del filamento de PET reciclado, predecirlo con un modelo viscoelástico y usarlo para evaluar un tratamiento térmico.

## Hipótesis

1. Un sistema de video con seguimiento de marcas puede medir el creep con la precisión necesaria para ajustar modelos viscoelásticos.
2. Un modelo viscoelástico calibrado con los primeros minutos de un ensayo puede predecir la deformación posterior.
3. El tratamiento del filamento en agua a 90 °C con estirado reduce su fluencia. El mecanismo propuesto es la orientación de las cadenas y la cristalización por deformación.

---

## Metodología

1. **Material.** Filamento de PET reciclado. Su área efectiva se calculó con la masa por unidad de largo: A = (m/L)/ρ, con ρ = 1.335 g/cm³. Este método es más confiable que el vernier en un filamento que no es perfectamente circular.
2. **Tratamiento.** Calentamiento del filamento en agua a ≈ 90 °C (por encima de la transición vítrea del PET, ≈ 75 °C) y estirado con pinzas. El largo pasó de 19 a 19.5 cm.
3. **Banco de creep.** El filamento cuelga de un soporte y lleva dos banderas blancas como marcas. Se cuelga una masa conocida y un celular fijo graba el ensayo (Figura 1).
4. **Medición por video.** El video se exporta a 5 cuadros/s. Python (OpenCV) localiza el centro de cada bandera en cada cuadro, como el punto medio entre su borde superior y su borde inferior, y calcula la distancia entre marcas L(t). La deformación es ε(t) = (L − L₀)/L₀.
5. **Modelado.** Se ajustaron los modelos de Burgers (resorte + Kelvin-Voigt + amortiguador) y una serie de Prony de dos términos con flujo. Para validarlos, el modelo se ajustó con los primeros 5 minutos y se comparó su predicción con los minutos 5 a 16.

---

## Resultados

**R1. Precisión del método (Figura 2).** Con el filamento ya estable, la medición varía ±0.007 % de deformación (±8 µm). Es una resolución comparable a la de un extensómetro clase B-1 según la norma ASTM E83 (0.01 %), lograda con un celular.

**R2. Creep del filamento sin tratar (Figura 3).** Con 1 kg (≈ 25 MPa), la deformación fue de 0.45 % a los 8 s y de 0.91 % a los 27 min. El modelo de Burgers describe la curva con R² = 0.98.

**R3. Cedencia retardada (Figura 4).** Con 4 kg (≈ 100 MPa) el filamento sin tratar resistía. Al agregar 0.5 kg (≈ 112 MPa) cedió de forma retardada: 1.1 % a los 15 s, 11.4 % a los 2 min y 14.1 % al romperse, a los 20 min. La velocidad máxima de deformación fue de 19 %/min a los 52 s.

**R4. Efecto del tratamiento a 90 °C (Figura 5).** Con la misma masa de 4.5 kg, el filamento tratado se deformó por creep 3.8 % en 16 min, contra 13 % del filamento sin tratar, sin ceder ni romperse.

**R5. Predicción del modelo (Figura 6).** Con solo los primeros 5 minutos, la serie de Prony predijo 6.11 % a los 16 min, contra 6.04 % medido (R² de validación = 0.98, error máximo de 0.08 puntos de %). El modelo de Burgers predijo 7.8 %, porque sobreestima el flujo permanente.

---

## Discusión

- El método por video alcanza una precisión suficiente para estudiar el creep del PET, y cuesta una fracción de un extensómetro comercial.
- El creep del PET reciclado no sigue un flujo lineal, como supone el modelo de Burgers: se desacelera con el tiempo. Por eso una serie de Prony con dos tiempos característicos (≈ 13 s y ≈ 4 min) lo describe y lo predice mejor.
- El filamento tratado resistió mejor la carga sostenida. Parte de la diferencia se explica porque trabajó a menor esfuerzo (≈ 91 contra ≈ 112 MPa), al ser más grueso.
- Las rupturas ocurrieron en los nudos de sujeción. La resistencia medida corresponde al amarre, no al material.

## Conclusiones

1. Es posible medir el creep del PET reciclado con un celular y Python, con una precisión de ±0.007 % de deformación.
2. Una serie de Prony calibrada con 5 minutos predice el creep posterior con R² = 0.98, lo que cumple el criterio de éxito planteado (R² ≥ 0.95). El modelo de Burgers no es adecuado para predecir a largo plazo.
3. El tratamiento en agua a 90 °C con estirado es consistente con una menor fluencia del filamento: 3.8 % contra 13 % con la misma carga.
4. El método permite evaluar la fiabilidad del PET reciclado sin equipo industrial, como aporte a la economía circular.

## Limitaciones

- Número reducido de ensayos, sin réplicas.
- El filamento tratado y el sin tratar no trabajaron al mismo esfuerzo, y el sin tratar había sido cargado antes en la prueba 2.
- El estirado fue pequeño (λ ≈ 1.03).
- No se midió la cristalinidad: la medición por Arquímedes resultó fuera de rango y debe repetirse.
- Las rupturas ocurrieron en los nudos.

## Siguientes pasos

- Medir la densidad por Arquímedes del filamento tratado y del sin tratar para calcular la cristalinidad (Xc = ρc(ρ − ρa)/[ρ(ρc − ρa)]).
- Estiramientos controlados (λ = 2 y 3) al mismo esfuerzo, con réplicas.
- Un amarre sin nudos (por ejemplo, enrollar el filamento en un perno).
- Ensayos más largos y con temperatura controlada.

---

## Frases cortas para el póster

- "Un celular y Python miden la deformación con una precisión de 0.007 %."
- "Con 5 minutos de datos, el modelo predice los siguientes 11 con R² = 0.98."
- "Con 4.5 kg: 13 % sin tratar contra 3.8 % tratado a 90 °C."
- "El PET reciclado sin tratar cede de forma retardada a ≈ 112 MPa."

---

## Preguntas probables de los jueces

**¿Cómo saben que el tratamiento aumentó la cristalinidad?**
Todavía no lo medimos. Lo planteamos como mecanismo probable, y la medición por densidad (Arquímedes) es nuestro siguiente paso. Lo que sí medimos es el efecto en el comportamiento: menos creep y sin cedencia.

**¿La comparación es justa si el esfuerzo fue distinto?**
No del todo: con la misma masa, el filamento tratado trabajó a ≈ 91 MPa y el sin tratar a ≈ 112 MPa. Por eso concluimos que el resultado es "consistente con" una mejora, y proponemos ensayos al mismo esfuerzo.

**¿Por qué Prony y no Burgers?**
Burgers supone un flujo constante que nunca se detiene. Nuestros datos muestran que el creep se desacelera. Prony usa dos tiempos característicos y predice mejor: 6.11 % contra 6.04 % medido, mientras que Burgers predice 7.8 %.

**¿Qué tan confiable es medir con un celular?**
La variación de la medición es de ±0.007 % de deformación, del orden de un extensómetro clase B-1. Además, se mide la distancia entre dos marcas, así que el movimiento del soporte o el deslizamiento del nudo no afectan el resultado.

**¿Por qué se rompió en el nudo?**
El nudo concentra el esfuerzo. Es una limitación del amarre, que proponemos corregir con una sujeción sin nudos.
