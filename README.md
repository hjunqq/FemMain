FemMain
=======
采用的并行模式如下：
do
{
	iiter++;
	
	//交换数据
	InteractValue=Fem.ExchangeData();
	
	//将交换的数据放入迭代过程
	Fem.SetInteractResult(InteractValue);
	
	//求解
	Fem.Solve();
	
	//检查收敛性
	Converge=Fem.ConvergeCheck();
	cout << "Iter=    " << iiter << endl;
	//后处理
	Fem.ComputeElementStress();
	Fem.CountElement();
	Fem.SendResultToNode();
	Fem.GIDOutResult(iiter);
	
}

//所有的进程达到收敛或者达到最大的迭代步之后退出迭代循环
while (Converge == false && iiter<MaxIter);

多个进程一直在迭代，知道所有的进程都达到收敛才退出迭代过程。
