#include <iostream>
#include "MITS.h"

using namespace std;

int main(int argc, char **argv)
{
#ifdef _WIN32
    // 设置控制台代码页为 UTF-8
    SetConsoleOutputCP(CP_UTF8);
    // 或者使用中文代码页
    // SetConsoleOutputCP(936);
#endif
    // bool tuning = false;
    	bool tuning = true;

    if (tuning)
    {
        root_path = "./";
        Read_Parameters(argc, argv);
        srand(seed);
        string rltfile = "./results/result_MLILS";
        execute_each_instance(filename, rltfile);
    }
    else
    {
        root_path = "../";
#ifdef _WIN32
        // 设置控制台代码页为 UTF-8
        SetConsoleOutputCP(CP_UTF8);
        // 或者使用中文代码页
        // SetConsoleOutputCP(936);
#endif

        srand(14); //14
        string dir = "..//Instances//";
        string instance_list = "instances_list.txt";

        Read_Parameters(argc, argv);
        string rltfile = to_string(seed) + "result_ILS";
        execute_each_instance(dir, instance_list, rltfile);
        cout << "\n!!!Hello World!!!\n" << endl; // prints !!!Hello World!!!
    }

    return 0;
}
