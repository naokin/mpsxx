#ifndef _PARSING_FCIDUMP_H
#define _PARSING_FCIDUMP_H

//
// parsing molpro FCIDUMP file
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <legacy/DENSE/TArray.h>

std::vector<std::string> gettoken(std::ifstream& fin);

void parsing_reorder
(std::ifstream& frord, std::vector<int>& reorder);

void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::TArray<double,2>& oneint);
void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::TArray<double,2>& oneint, const std::vector<int>& reorder);

void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::TArray<double,4>& twoint);
void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::TArray<double,4>& twoint, const std::vector<int>& reorder);


void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::TArray<double,2>& oneint, btas::TArray<double,4>& twoint);
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::TArray<double,2>& oneint, btas::TArray<double,4>& twoint, const std::vector<int>& reorder);

void writing_fcidump
(std::ofstream& fdump, const int& norbs, const int& nelec, const double& ecore, const btas::TArray<double,2>& oneint, const btas::TArray<double,4>& twoint);

#endif // _PARSING_FCIDUMP_H
