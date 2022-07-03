# RISC-V Simulator
with Tomasulo Algorithm and 2-bit branch predictor

设计图：
![](https://raw.githubusercontent.com/Jerx2y/RISC-V-Simulator-Assignment-2022/master/design.jpg)

| 数据点 | 预测直接跳转 | 二位饱和预测 | 二级预测 |
| -------- | -------- | -------- | -------- |
| array_test1* | 0.454545 | 0.454545 | 0.5 |
| array_test2* | 0.5 | 0.461538 | 0.5 | 
| basicopt1 | 0.410226 | 0.696085 | 0.801101 | 
| bulgarian | 0.493601 | 0.90262 | 0.906201 | 
| expr* | 0.621622 | 0.81982 | 0.72973 | 
| gcd* | 0.375 | 0.625 | 0.658333 | 
| hanoi | 0.5002 | 0.73535 | 0.617403 | 
| lvalue2* | 0.333333 | 0.666667 | 0.666667 | 
| magic | 0.468564 | 0.712623 | 0.707215 | 
| manyarguments* | 0.2 | 0.8 | 0.8 | 
| multiarray* | 0.814815 | 0.740741 | 0.617284 | 
| naive* | / | / | / | 
| pi | 0.577298 | 0.783479 | 0.79273 | 
| qsort | 0.699933 | 0.833028 | 0.86584 | 
| queens | 0.367784 | 0.603649 | 0.663001 | 
| statement_test* | 0.405941 | 0.59901 | 0.613861 | 
| superloop | 0.127054 | 0.911658 | 0.892184 | 
| tak | 0.250004 | 0.800112 | 0.740332 | 