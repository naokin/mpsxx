#ifndef _PARSING_FCIDUMP_H
#define _PARSING_FCIDUMP_H

//
// parsing molpro FCIDUMP file
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <btas/DENSE/DArray.h>

std::vector<std::string> gettoken(std::ifstream& fin);

void parsing_reorder
(std::ifstream& frord, std::vector<int>& reorder);

void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::DArray<2>& oneint);
void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::DArray<2>& oneint, const std::vector<int>& reorder);

void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::DArray<4>& twoint);
void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::DArray<4>& twoint, const std::vector<int>& reorder);


void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::DArray<2>& oneint, btas::DArray<4>& twoint);
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::DArray<2>& oneint, btas::DArray<4>& twoint, const std::vector<int>& reorder);

void writing_fcidump
(std::ofstream& fdump, const int& norbs, const int& nelec, const double& ecore, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint);

#endif // _PARSING_FCIDUMP_H
