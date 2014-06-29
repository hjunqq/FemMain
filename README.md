FemMain
=======
把Converge变量存储为一个含有两个元素的数组，
Converge[0]存储的是本地的收敛条件
Converge[1]存储的是远程的收敛条件
do 
		{
			iiter++;
			
			InteractValue=Fem.ExchangeData();
			Fem.SetInteractResult(InteractValue);
			如果本地收敛，就不需要再迭代了。
			if (Converge[0] == false)
			{
				Fem.Solve();
			}

			Converge = Fem.ConvergeCheck();

			cout << "Iter=    " << iiter << endl;

			Fem.ComputeElementStress();

			Fem.CountElement();

			Fem.SendResultToNode();

			Fem.GIDOutResult(iiter);


		} while (Converge[1] == false && iiter<MaxIter);

多个进程一直在迭代，直到所有的进程都达到收敛才退出迭代过程。

