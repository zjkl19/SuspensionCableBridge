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

import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sympy
#x = symbols('x')
h=0;f=9;l=70;
n=f/l;  #矢跨比

x = sympy.symbols("x")
y=(h/l)*x+(4*f/l**2)*x*(l-x)    #p57 3-2-11

#S:有应力索长
#过溪桥：S=72.97356
S=float(sympy.integrate((1+sympy.diff(y,x)**2)**0.5, (x, 0, l)))    #P61 3-2-30，积分表达式
logger.info('S:'+str(S))

#p61 3-2-30，积分展开结果
#S=(l/2)*(1+16*n**2)**0.5+0.125*(l/n)*math.log(4*n+(1+16*n**2)**0.5);
#S=l*(1+8*n**2/3-32*n**4/5);     #有应力索长估算公式（按级数展开结果）以过溪桥为例，误差越1

from scipy.optimize import fsolve

c=sympy.symbols("c")
#s =sympy.solve(-c*sympy.cosh(0.5*l/c-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)+c*sympy.cosh(-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)-(h/2+n), c)
#print(float(s))
def balanceFun(z):
    c=float(z[0])
    return [-c*sympy.cosh(0.5*l/c-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)+c*sympy.cosh(-sympy.asinh(0.5*h/(c*sympy.sinh(0.5*l/c)))-0.5*l/c)-(h/2+n)]

result=fsolve(balanceFun,[1])
logger.info('result:\n'+str(result))

uy=-result[0]*sympy.cosh(x/result[0]-sympy.asinh(0.5*h/(result[0]*sympy.sinh(0.5*l/result[0])))-0.5*l/result[0])+result[0]*sympy.cosh(-sympy.asinh(0.5*h/(result[0]*sympy.sinh(0.5*l/result[0])))-0.5*l/result[0])   #无应力索线形unstressed y
logger.info('uy:\n'+str(uy))
#def balanceFun(z):
#    c=float(z[0])
#    c1=float(z[1])
#    c2=float(z[2])
#    return [-c*math.cosh((70/2)/c+c1)+c2-(0/2+9/70),
#            c*math.cosh(c1)-c2,
#            -math.asinh(0*(2*c*math.sinh(0.5*70/c))**-1)-0.5*70/c-c1]
    #return [-c*math.cosh((l/2)/c+c1)+c2-(h/2+n),
    #        c*math.cosh(c1)-c2,
    #        -math.asinh(h*(2*c*math.sinh(0.5*l/c))**-1)-0.5*l/c-c1]
#            -math.asinh(h*(2*c*math.sinh(0.5*l/c))**-1)-0.5*l/c-c1]
#    return [-c*math.cosh((l/2)/c+c1)+c2-(h/2+n),
#result=fsolve(balanceFun,[7.17,24677,504.246])
#logger.info('result:\n'+str(result))

#空缆状态下主缆有应力索长
Ec=2.0*10**11;   #主缆弹性模量(国际单位）
Ac=0.25*3.14*0.10695**2;   #主缆面积（国际单位）

gb=168000;   #加劲梁自重

H=gb*l**2/(8*f);

deltaS1=(H*l)*(1+(16*n**2/3))/(Ec*Ac);
S1=S-deltaS1;

logger.info('S1:'+str(S1))

deltaS11=H*S/(Ec*Ac);
S11=S-deltaS11;
logger.info('S11:'+str(S11))

#无应力索长
gc=8.348*10**4;   #缆索自重
c=l**2/(8*f);
logger.info(str(c))
c=69.5058;
Hc=c*gc;    #缆索在自重状态下的水平拉力

deltaS2=Hc*(l+c*math.sinh(l/c))/(2*Ec*Ac);

#def fff(x):
#    return (1-x**2)**(1/2)

#w, err = integrate.quad(f,-1,1)
#print(w,err)

#float(sympy.integrate((1+sympy.diff(y,x)**2)**0.5, (x, 0, l)))
#float(sympy.integrate((1+sympy.diff(uy,x)**2)**0.5, (x, 0, l)))
logger.info('deltaS2:'+str(deltaS2))
deltaS22=Hc*float(sympy.integrate(1+(sympy.diff(uy,x))**2, (x, 0, l)))/(Ec*Ac)
logger.info('deltaS22:'+str(deltaS22))

from sympy import  integrate ,cos,sin
from sympy.abc import  a,x,y
r=integrate(sin(x)/x,(x,-float("inf"),float("inf")))

r1=integrate(sympy.cosh(x))
logger.info('r1:'+str(r1))
S0=S-deltaS1-deltaS2;
