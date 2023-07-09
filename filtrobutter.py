import numpy as np
from scipy.signal import butter, iirfilter, butter, sosfilt, sos2tf, bilinear, zpk2tf,freqs,dlti,freqz, lp2hp,lti,lp2bp
import matplotlib.pyplot as plt

ordem = 6  # Ordem do filtro
cutoff_freq = (2/3) * np.pi  # Frequência de corte em radianos / segundo
fc = 44100/3
fs = 44100 # Frequência de amostragem para sinal de audio
wn = fs*cutoff_freq

# Projeto do filtro Butterworth analógico
butter_numerador, butter_denominador = iirfilter(ordem, wn,btype = 'lowpass',analog=True, ftype='butter')
print("Numerador:", butter_numerador)
print("Denominador:", butter_denominador)

# Frequência de amostragem
sampling_freq = 44100

# Transformação bilinear para filtro IIR na forma direta
filtz = dlti(*bilinear(butter_numerador, butter_denominador, fs=sampling_freq))

# Coeficientes do filtro em ponto flutuante
print("Forma direta (ponto flutuante):")
print("b:", filtz.num)
print("a:", filtz.den)

# Projeto do filtro Butterworth analógico
butter_sos = butter(ordem, wn, analog=True, output='sos')
for sos in butter_sos:
  print("Bloco SOS: ", sos)

from scipy.signal import butter, sosfilt, sos2tf, bilinear, zpk2tf
import numpy as np

# Frequência de amostragem normalizada
sampling_freq = 44100

z,p,k = butter(ordem, wn, analog=True, output='zpk')

b,a = zpk2tf(z,p,k)

b_bili,a_bili = bilinear(b,a,fs=sampling_freq)

print(b_bili)
print(a_bili)

# Transformação bilinear para cada seção do filtro em cascata
#for sos in butter_sos:
 # b, a = sos2tf(sos)

#print(b)
#print(a)

# Coeficientes do filtro em ponto flutuante
#print("Forma direta (ponto flutuante):")
#print("b:", bilinear_numerador)
#print("a:", bilinear_denominador)

import numpy as np

# Quantização dos coeficientes
decimal_points = [8, 6, 4, 2]  # Número de casas decimais após a vírgula para quantização
quantized_b_direto = []
quantized_a_direto = []

for dp in decimal_points:
    quantized_b_direto.append(np.round(filtz.num, dp))
    quantized_a_direto.append(np.round(filtz.den, dp))

for i in range(0, len(quantized_b_direto)):
  print("\nForma direta (", decimal_points[i], " casas decimais):")
  print("b:", quantized_b_direto[i])
  print("a:", quantized_a_direto[i])


  for dp in decimal_points:
      quantized_sos = np.round(sos, dp)
      print("\nForma cascata (", dp, " casas decimais):")
      print("a:", quantized_sos)

w, h = freqs(butter_numerador, butter_denominador)

# plt.semilogx(w, 20 * np.log10(abs(h)))

# plt.title('Butterworth filter frequency response')

# plt.xlabel('Frequency [radians / second]')

# plt.ylabel('Amplitude [dB]')

# plt.margins(0, 0.1)

# plt.grid(which='both', axis='both')

# plt.axvline(100, color='green') # cutoff frequency

#plt.show()


w, h = freqs(butter_numerador, butter_denominador)

# plt.semilogx(w, (abs(h)))

# plt.title('Butterworth filter frequency response')

# plt.xlabel('Frequency [radians / second]')

# plt.ylabel('Amplitude [dB]')

# plt.margins(0, 0.1)

# plt.grid(which='both', axis='both')

# plt.axvline(100, color='green') # cutoff frequency

#plt.show()


ws2,hs2 = freqs(butter_numerador, butter_denominador)

wsHz2 = ws2/2*np.pi

# plt.plot(wsHz2*fs/(2*np.pi),20 * np.log10(np.abs(hs2)))
# plt.grid()
# plt.show()


bz,az = freqz(filtz.num,filtz.den)

plt.plot(bz*fs/(2*np.pi), 20*np.log10(np.abs(az)))
plt.axvline(x=fc, color='red', ls='--', label='$f_c$ = {} Hz'.format(round(fc)), alpha=0.7)
plt.legend()
plt.title('Transformação Bilinear')
plt.xlabel('Freq [Hz]')
plt.ylabel('Magnitude [dB]')
plt.grid()
plt.show()


#bz,az = freqz(filtz.num,filtz.den)

hp = lti(*lp2bp(filtz.num,filtz.den))
bzh,azh = freqz(hp.num,hp.den)

plt.plot(bzh*fs/(2*np.pi), 20*np.log10(np.abs(azh)))
plt.axvline(x=fc, color='orange', ls='--', label='$f_c$ = {} Hz'.format(round(fc)), alpha=0.7)
plt.legend()
plt.title('Transformação Bilinear')
plt.xlabel('Freq [Hz]')
plt.ylabel('Magnitude [dB]')
plt.grid()
plt.show()




