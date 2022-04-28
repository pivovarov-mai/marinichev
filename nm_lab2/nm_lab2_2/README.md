# Анализ системы

Посмотрим на график системы:

<img src="https://render.githubusercontent.com/render/math?math=$\begin{cases} x_1 - \cos(x_2) = 3 \\ x_2 - \sin(x_1) = 3 \end{cases}$">


```python
import matplotlib.pyplot as plt
import numpy as np

X1 = np.linspace(-10,5,100)
X2 = np.linspace(-5, 6, 100)

# вариант 12
x1 = np.cos(X2) + 3 
x2 = np.sin(X1) + 3

# область G
a, b = 2, 3.5
c, d = 2.5, 4

plt.plot(x1, X2, color="blue")
plt.plot(X1, x2, color="red")
plt.axvline(a, linestyle='--')
plt.axhline(b, linestyle='--')
plt.axvline(c, linestyle='--')
plt.axhline(d, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y = 0, color='k')
plt.axvline(x = 0, color='k')
plt.axis('equal')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_2/img/output_1_0.png)
    


Посмотрим на решение. Возьмем область <img src="https://render.githubusercontent.com/render/math?math=$G: x_1\in[2, 2.5], x_2\in[3.5, 4]$"> для его поиска


```python
X1 = np.linspace(1.9, 2.6,100)
X2 = np.linspace(3.3, 4.2, 100)

x1 = np.cos(X2) + 3 
x2 = np.sin(X1) + 3 

plt.plot(x1, X2, color="blue")
plt.plot(X1, x2, color="red")
plt.axvline(a, linestyle='--')
plt.axhline(b, linestyle='--')
plt.axvline(c, linestyle='--')
plt.axhline(d, linestyle='--')
plt.grid(True, which='both')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_2/img/output_3_0.png)
    


# Численное решение

Теперь запустите в терминале программу на **С++** при помощи команд
```
make
make run
```
В трех появившихся файлах `answer_NN.txt` можно увидеть полученное численное решение c разной точностью. Проверим его


```python
def f1(x1, x2):
    return x1 - np.cos(x2) - 3

def f2(x1, x2):
    return x2 - np.sin(x1) - 3

x1_0 = 2.21045
x2_0 = 3.80231

print("Значение функций f1(x1, x2) и f2(x1, x2) в найденной точке:")
print(f"f1{x1_0, x2_0} = {f1(x1_0, x2_0):0.2f}")
print(f"f2{x1_0, x2_0} = {f2(x1_0, x2_0):0.2f}")
```

    Значение функций f1(x1, x2) и f2(x1, x2) в найденной точке:
    f1(2.21045, 3.80231) = 0.00
    f2(2.21045, 3.80231) = 0.00
    

Наконец, построим полученную точку на графике


```python
X1 = np.linspace(1.9, 2.6,100)
X2 = np.linspace(3.3, 4.2, 100)

x1 = np.cos(X2) + 3 
x2 = np.sin(X1) + 3 

plt.plot(x1, X2, color="blue")
plt.plot(X1, x2, color="red")
plt.grid(True, which='both')
plt.scatter(x1_0, x2_0, c ="green")
plt.show()
```


    
![png](/nm_lab2/nm_lab2_2/img/output_7_0.png)
    

