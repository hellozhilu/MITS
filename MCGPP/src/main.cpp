//============================================================================
// Name        : MCGP.cpp
// Author      : Lin
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
//#include "test_ILS_v1.3.2.h"
//#include "test_infeasible_v1.1.h"
//#include "test_ML_v1.1.0.h"
//#include "test_ILSWT_v2.3.0.h"
//#include "test_MLWT_v2.2.3.h"
//#include "test_MLWT_v2.2.3.h"
//#include "MLWT.h"
//#include "test_MLWT_v2.2.2.3.h"
//#include "MLWT_1222_modified.h"			// 目前效果最好的版本
//#include "DGCWT_v1.0.h"

#include "1-MLWT_convergence.h"

using namespace std;

int main(int argc, char** argv)
{
//	bool tuning = false;
	bool tuning = true;

	if(tuning)
	{
		Read_Parameters(argc, argv);
		srand(seed);
		string rltfile = "./results/result_MLILS";
		execute_each_instance(filename, rltfile);
	}
	else
	{
		srand(14);					//14
		string dir = ".//Instances//";
		string instance_list = "instances_list.txt";

		Read_Parameters(argc, argv);
		string rltfile = to_string(seed) + "result_ILS";
		execute_each_instance(dir, instance_list, rltfile);
		cout << "\n!!!Hello World!!!\n" << endl; // prints !!!Hello World!!!
	}

	return 0;
}
