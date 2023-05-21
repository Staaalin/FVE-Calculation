# 说明
* BC_Sum
  * 通过 “from BC_Sum import BC_Sum” 将该函数导入至需要的程序中。
  * 使用说明：以积分 $$\int \frac{d^3 p}{(2 \pi)^3} \frac{1}{p^4 + 1}$$ 为例，对其进行Matsubara频率反周期边界条件的求和替换 $$p_x \rightarrow (2n_x+1)\pi/L$$ $$p_y \rightarrow (2n_y+1)\pi/L$$ $$p_z \rightarrow (2n_z+1)\pi/L$$ $$\frac{d^3 p}{(2 \pi)^3} \rightarrow \frac{1}{L^3}\sum_{n_x, n_y, n_z=-\infty}^{\infty}$$ 积分式就变为了 $$\int \frac{d^3 p}{(2 \pi)^3} \frac{1}{p^4 + 1} \rightarrow \frac{1}{L^3}\sum_{n_x, n_y, n_z=-\infty}^{\infty} [((\frac{2n_x+1}{L}\pi)^2 + (\frac{2n_y+1}{L}\pi)^2 + (\frac{2n_z+1}{L}\pi)^2)^2+1]^{-1}$$ 箭头左边的数值计算是容易的，而箭头右边的数值计算肉眼可见得更麻烦。而本程序中"BC_Sum"函数通过优化求和过程的重复项，可以以最快的速度和效率计算箭头右边的这串求和。在导入该函数后，具体操作如下

```python
import numpy as np
import sys
import math

from BC_Sum import BC_Sum # 加载BC_Sum

# 定义被积函数
def f(p):
    return 1/(p**4 + 1)

# 定义参量
L = 1000 # fm
Lm = L/197.32 # 转换为自然单位制，L的单位变为MeV^{-1}
p2_limit = 100 # 该积分虽然在无穷大处是收敛的，但是转化成求和不可能求和至无穷远处，因此我们设置上限为100，积分与求和的空间动量的模长都不超过p2_limit

# 利用BC_Sum计算求和结果（等于上面公式中箭头右边）
r = BC_Sum(f,'aPBC',Lm=Lm,p2_Limit=p2_limit)
print(r)
```

得到结果

```
0.05609531200914537
```

以上是"BC_Sum"函数的使用方法，对于周期边界条件 $$p_x \rightarrow (2n_x)\pi/L$$ $$p_y \rightarrow (2n_y)\pi/L$$ $$p_z \rightarrow (2n_z)\pi/L$$ 可以调整参数“aPBC”为“PBC”实现

```python
r = BC_Sum(f,'aPBC',Lm=Lm,p2_Limit=p2_limit) # 反周期
r = BC_Sum(f,'PBC',Lm=Lm,p2_Limit=p2_limit)  # 周期
```

得到结果

```
0.05393734659017124
```

下面可以比较，使用优化求和过程的重复项的"BC_Sum"函数，和不优化的直接求和的计算速度差异。以周期为例：

```python
import numpy as np
import sys
import math

from scipy.integrate import quad # 不必要，只是在这里用在展示
import time # 不必要，只是在这里用在展示

from BC_Sum import BC_Sum # 加载BC_Sum

# 定义被积函数
def f(p):
    return 1/(p**4 + 1)

# 定义参量
L = 1000 # fm
Lm = L/197.32 # 转换为自然单位制，L的单位变为MeV^{-1}
p2_limit = 100 # 该积分虽然在无穷大处是收敛的，但是转化成求和不可能求和至无穷远处，因此我们设置上限为100，积分与求和的空间动量的模长都不超过p2_limit

# 在做Matsubara频率求和之前的积分结果，这里使用scipy.integrate中的quad函数
r1 = quad(lambda p:4*np.pi*p**2*f(p),0,p2_limit)[0]/(2*np.pi)**3

# 不优化的直接求和，赋值给r2
LL = int((p2_limit*Lm/np.pi)/2)
s1=time.time()
r2 = 0
for i in range(-LL,LL+1):
    for j in range(-LL,LL+1):
        for k in range(-LL,LL+1):
            rs = ((((2*i)*np.pi/Lm)**2 + ((2*j)*np.pi/Lm)**2 + ((2*k)*np.pi/Lm)**2)**0.5)
            if rs < p2_limit:
                r2 += f(rs)
r2 = r2/Lm**3
s2=time.time()

# 使用BC_Sum求和，赋值给r3
r3 = BC_Sum(f,'PBC',Lm=Lm,p2_Limit=p2_limit)
s3=time.time()

print(r1)
print(r2)
print(r3)
print('____________')
print(s2-s1) # 计算r2的速度
print(s3-s2) # 计算r3的速度
```

得到结果：

```
0.05576316384262066
0.053937346590194474
0.05393734659017124
____________
3.542262077331543
0.0845489501953125
```

可见在不损失精度的情况下，BC_Sum函数的计算速度快得多。同时也可以看出在大L的情况下，积分与求和的结果相近。
