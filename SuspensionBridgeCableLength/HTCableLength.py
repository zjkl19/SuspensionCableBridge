# -*- coding: utf-8 -*-
"""
洪塘大桥
矢高f=15.114
l=150
h=80.107-23.732=56.375
Created on Tue Mar 12 10:12:50 2019

@author: Lindinan
"""
import math
import sympy
from scipy.optimize import fsolve
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


h=80.107-23.732    #y|x=l 处高度
f=15    #矢高
l=150    #跨径
n=f/l  #矢跨比

Ec=2.0*10**11   #主缆弹性模量(国际单位）
Ac=0.25*3.14*0.35644**2   #主缆面积（国际单位）
gb=182734.2+142.2*1000+40.0*1000   #加劲梁自重（P58 有解释）。计算过程详见excel表格（联合midas计算）
H=gb*l**2/(8*f)    #主缆拉力的水平分力（P57 3-2-10）
print('成桥阶段水平力H:'+str(H))
gc=78.2*10**3*Ac   #缆索沿纵桥向每延米自重

print (gc/gb)

c=sympy.symbols("c")
#求解空缆线形方程所用的平衡方程
#联立P58 3-2-16~3-2-18所得
def balanceFun(z):
    c=float(z[0])
    return [-c*sympy.cosh(0.5*l/c-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)+c*sympy.cosh(-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)-(h/2+n)]

x = sympy.symbols("x")
#成桥状态线形方程
y=(h/l)*x+(4*f/l**2)*x*(l-x)    #p57 3-2-11
print('y:'+str(y))

coor=[0,10,20,30,40,50,60,70,80,90,100,110,150]

#通过序列项迭代
for i in coor:
    print(y.subs(x, i))

#S:有应力索长
#过溪桥（参考）：S=72.97356
S=float(sympy.integrate((1+sympy.diff(y,x)**2)**0.5, (x, 0, l)))    #P61 3-2-30，积分表达式
print('主缆有应力索长S:'+str(S))

#以下2式仅适用于h=0的结果（p61 3-2-30，积分展开结果）
#S=(l/2)*(1+16*n**2)**0.5+0.125*(l/n)*math.log(4*n+(1+16*n**2)**0.5);
#S=l*(1+8*n**2/3-32*n**4/5);     #有应力索长估算公式（按级数展开结果）以过溪桥为例，误差越1

result=fsolve(balanceFun,[1])    #迭代初值为1
logger.info('平衡方程求解结果:'+str(result))

uy=-result[0]*sympy.cosh(x/result[0]-sympy.asinh(0.5*h/(result[0]*sympy.sinh(0.5*l/result[0])))-0.5*l/result[0])+result[0]*sympy.cosh(-sympy.asinh(0.5*h/(result[0]*sympy.sinh(0.5*l/result[0])))-0.5*l/result[0])   #无应力索线形unstressed y
logger.info('空缆状态下的线形方程uy:\n'+str(uy))

coor=[0,10,20,30,40,50,60,70,80,90,100,110,150]

#通过序列项迭代
for i in coor:
    print(uy.subs(x, i))


#空缆状态下主缆有应力索长

#P61 3-2-34
#感觉书上的公式有误，3-2-31代入3-2-34，好像结果不对
#以下2式仅适用于h=0的结果
#deltaS1=(H*l)*(1+(16*n**2/3))/(Ec*Ac)
#S1=S-deltaS1
#logger.info('S1:'+str(S1))

deltaS1=H*S/(Ec*Ac)
print('主缆在加劲梁自重作用下的弹性伸长值deltaS1:'+str(deltaS1))
S1=S-deltaS1
print('空缆状态下的缆长S1:'+str(S1))

#无应力索长

#适用于h=0情况
#c=l**2/(8*f);
#logger.info(str(c))
#c=69.5058;
#Hc=c*gc;    #缆索在自重状态下的水平拉力
#deltaS2=Hc*(l+c*math.sinh(l/c))/(2*Ec*Ac);

#float(sympy.integrate((1+sympy.diff(y,x)**2)**0.5, (x, 0, l)))
#float(sympy.integrate((1+sympy.diff(uy,x)**2)**0.5, (x, 0, l)))
#logger.info('deltaS2:'+str(deltaS2))

#P62 3-2-36
Hc=result[0]*gc
Hc=1700000
print('空缆阶段水平力Hc:'+str(Hc))
deltaS2=Hc*float(sympy.integrate(1+(sympy.diff(uy,x))**2, (x, 0, l)))/(Ec*Ac)
print('缆索自重作用下的钢缆弹性伸长值deltaS2:'+str(deltaS2))

S0=S-deltaS1-deltaS2;
print('主缆无应力索长S0:'+str(S0))
