算法过程：
1.设置各参数

2.随机产生初始群体——pop=initpop(popsize,chromlength)

3.故障类型编码，每一行为一种！code(1,:)，正常；code(2,:)，50％；code(3,:)，150％。实际故障测得数据编码，这里Unnoralcode,188%

4.开始迭代（M次）：
     1）计算目标函数值：欧氏距离[objvalue]=calobjvalue(pop,i)
     2）计算群体中每个个体的适应度fitvalue=calfitvalue(objvalue)
     3）选择newpop=selection(pop,fitvalue); objvalue=calobjvalue(newpop,i); %
        交叉newpop=crossover(newpop,pc,k);  objvalue=calobjvalue(newpop,i); %
        变异newpop=mutation(newpop,pm);     objvalue=calobjvalue(newpop,i); %

5.求出群体中适应值最大的个体及其适应值

6.迭代停止判断。