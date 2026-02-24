#ifndef TOOLS_H_
#define TOOLS_H_

double myRound(double value) {
	unsigned decimals = 6;
	double factor = pow(10,decimals);
	return round(value*factor)/factor;
}

double myRound(double value, int decimals) {
	double factor = pow(10,decimals);
	return round(value*factor)/factor;
}

/*
 * Fisher-Yates 洗牌算法
 */
void Generate_Rand_List1(int *randlist, int len)
{
    assert(randlist != NULL);

    // 洗牌算法，从后往前打乱
    for (int i = len - 1; i > 0; i--) {
        int randid = rand() % (i + 1); // 保证随机索引范围是 [0, i]
        int tmp = randlist[i];
        randlist[i] = randlist[randid];
        randlist[randid] = tmp;
    }
}

/*
 * 生成[0,len-1]的随机数列
 */
void Generate_Rand_List(int *randlist, int len)
{
	assert(randlist != NULL);

	for (int i = 0; i < len; i++)
		randlist[i] = i;
//	printf("打乱前");
//	for(int i = 0; i < len; i++)
//	{
//		printf("%d ", randlist[i]);
//		fflush(stdout);
//	}
//	printf("\n");
//	fflush(stdout);

	// swap
	for (int i = 0; i < len; i++)
	{
		int randid = rand() % len;
		int tmp = randlist[i];
		randlist[i] = randlist[randid];
		randlist[randid] = tmp;
	}

//	printf("打乱后");
//	for(int i = 0; i < len; i++)
//	{
//		printf("%d ", randlist[i]);
//		fflush(stdout);
//	}
//	printf("\n");
//	fflush(stdout);
}

void Quick_Sort_up(int *idlist, int *objlist, int l, int r)
{
	if (l < r)
	{
		int id = l, jd = r;
		int x = objlist[l];							// MDFD 原本是double类型数据
		int y = idlist[l];
		while (id < jd)
		{
			// Find first > x, from right to left
			while (id < jd && objlist[jd] >= x)
				jd--;
			if (id < jd)
			{
				objlist[id] = objlist[jd];
				idlist[id++] = idlist[jd];
			}

			// Find first <= x, from left to right
			while (id < jd && objlist[id] < x)
				id++;
			if (id < jd)
			{
				objlist[jd] = objlist[id];
				idlist[jd--] = idlist[id];
			}
		}
		objlist[id] = x;
		idlist[id] = y;
		Quick_Sort_up(idlist, objlist, l, id - 1);
		Quick_Sort_up(idlist, objlist, id + 1, r);
	}
}


#endif /* TOOLS_H_ */
