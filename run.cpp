#include <mpi.h>
#include <stdio.h>
#include <bits/stdc++.h> 
using namespace std;
int main()
{
	string str = "ctest -V -R TestColonicCryptSimulation"; 
	const char *command = str.c_str(); 
	system(command);

	return 0;
}
