# Анализ функции 

Посмотрим на график функции <img src="https://render.githubusercontent.com/render/math?math=$f(x)$">


```python
import matplotlib.pyplot as plt
from scipy.misc import derivative
import numpy as np

x = np.linspace(-2.5, 4.5, 100)
f = pow(3, x) - 5 * x**2 + 1 # вариант 12

plt.plot(x, f,color="red")
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_1_0.png)
    


Посмотрим на первое решение. Возьмем отрезок `[-1, -0.25]` для его поиска


```python
a1, b1 = -1, -0.25

x = np.linspace(-1.5, 0, 100)
f = pow(3, x) - 5 * x**2 + 1

plt.plot(x, f,color="red")
plt.axvline(a1, linestyle='--')
plt.axvline(b1, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_3_0.png)
    


Посмотрим на второе решение. Возьмем отрезок `[0.5, 1]` для его поиска


```python
a2, b2 = 0.5, 1

x = np.linspace(0, 1.5, 100)
f = pow(3, x) - 5 * x**2 + 1

plt.plot(x, f,color="red")
plt.axvline(a2, linestyle='--')
plt.axvline(b2, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_5_0.png)
    


Посмотрим на третье решение. Возьмем отрезок `[3.5, 4.2]` для его поиска


```python
a3, b3 = 3.5, 4.2

x = np.linspace(3, 4.5, 100)
f = pow(3, x) - 5 * x**2 + 1

plt.plot(x, f,color="red")
plt.axvline(a3, linestyle='--')
plt.axvline(b3, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_7_0.png)
    


# Первая производная: <img src="https://render.githubusercontent.com/render/math?math=$f'(x)$">

Для нас важно, чтобы на выбранных выше отрезках знак производной не менялся. Убедимся в этом


```python
def function(x):
    return pow(3, x) - 5 * x**2 + 1

def function_diff(x):
    return derivative(function, x)

x = np.linspace(-2.5, 4.5, 100)

plt.plot(x, function_diff(x), color="red")
plt.axvline(a1, linestyle='--')
plt.axvline(b1, linestyle='--')
plt.axvline(a2, linestyle='--')
plt.axvline(b2, linestyle='--')
plt.axvline(a3, linestyle='--')
plt.axvline(b3, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_9_0.png)
    


# Производная функции <img src="https://render.githubusercontent.com/render/math?math=$\phi(x)$">

Для нас важно, чтобы <img src="https://render.githubusercontent.com/render/math?math=$q =\max_{{x\in[a, b]}} |\phi'(x)| < 1$">


```python
def phi(x):
    return np.sqrt((pow(3, x) + 1) / 5)

def phi_diff(x):
    return derivative(phi, x)

x = np.linspace(-2, 4.5, 100)

plt.plot(x, phi_diff(x), color="red")
plt.axvline(a1, linestyle='--')
plt.axvline(b1, linestyle='--')
plt.axvline(a2, linestyle='--')
plt.axvline(b2, linestyle='--')
plt.axvline(a3, linestyle='--')
plt.axvline(b3, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axhline(y=1, color='b')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_11_0.png)
    


Как можно заметить нам подхдят только два первых отрезка. Тогда попробуем выразить <img src="https://render.githubusercontent.com/render/math?math=$\phi(x)$"> иначе


```python
def log(a, b):
    return np.log(b) / np.log(a)

def phi(x):
    return log(3, 5 * x**2 - 1)

def phi_diff(x):
    return derivative(phi, x)

x = np.linspace(2, 5, 100)

plt.plot(x, phi_diff(x), color="red")
plt.axvline(a3, linestyle='--')
plt.axvline(b3, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axhline(y=1, color='b')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_13_0.png)
    


Теперь условие выполняется для третьего отрезка. Соответсвенно, для первых двух отрезков будем использовать первый вариант выражения <img src="https://render.githubusercontent.com/render/math?math=$\phi(x)$">, а для третьего второй вариант <img src="https://render.githubusercontent.com/render/math?math=$\phi(x)$">

# Вторая производная: <img src="https://render.githubusercontent.com/render/math?math=$f''(x)$">

Для нас важно, чтобы на выбранных выше отрезках знак второй производной не менялся. Убедимся в этом


```python
def function_second_diff(x):
    return derivative(function_diff, x)

x = np.linspace(-2.5, 4.5, 100)

plt.plot(x, function_second_diff(x), color="red")
plt.axvline(a1, linestyle='--')
plt.axvline(b1, linestyle='--')
plt.axvline(a2, linestyle='--')
plt.axvline(b2, linestyle='--')
plt.axvline(a3, linestyle='--')
plt.axvline(b3, linestyle='--')
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_16_0.png)
    


# Численное решение

Теперь запустите в терминале программу на **С++** при помощи команд
```
make
make run
```
В трех появившихся файлах `answer_NN.txt` можно увидеть полученные численные решения. Проверим их


```python
x_1 = -0.555548
x_2 =  0.837941
x_3 =  3.95759

print("Значение функции f(x) в найденных точках:")
print(f"1) f({x_1}) = {function(x_1):0.2f}")
print(f"2) f({x_2}) = {function(x_2):0.2f}")
print(f"3) f({x_3}) = {function(x_3):0.2f}")
```

    Значение функции f(x) в найденных точках:
    1) f(-0.555548) = 0.00
    2) f(0.837941) = 0.00
    3) f(3.95759) = 0.00
    

Наконец, построим полученные точки на графике


```python
x_0 =[x_1, x_2, x_3]
 
y_0 =[function(x_1), function(x_2), function(x_3)]

x = np.linspace(-2.5, 4.5, 100)
f = pow(3, x) - 5 * x**2 + 1

plt.plot(x, f,color="red")
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.scatter(x_0, y_0, c ="blue")
plt.show()
```


    
![png](/nm_lab2/nm_lab2_1/img/output_20_0.png)
    

