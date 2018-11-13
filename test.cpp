/********************************************************************
 ********************************************************************
 ** O(n^3) implementation derived from libhungarian by Cyrill Stachniss, 2004
 **
 ** C++ class implementation of the Hungarian algorithm by David Schwarz, 2012.
 **
 ** Minor modifications by Justin Buchanan
 **
 ** Solving the Minimum Assignment Problem using the
 ** Hungarian Method.
 **
 ** ** This file may be freely copied and distributed! **
 **
 **
 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 ** PURPOSE.
 **
 ********************************************************************
 ********************************************************************/

#include <iostream>
#include "hungarian.hpp"

using namespace std;

const Hungarian::Matrix EXAMPLE1 = {
    {100, 100, 1},
    {100, 2, 21512},
    {1, 4, 9852},
    {6, 30252, 400},
};

const Hungarian::Matrix SOLUTION1 = {
    {0, 0, 1, 0},
    {0, 1, 0, 0},
    {1, 0, 0, 0},
    {0, 0, 0, 1},
};

// clang-format off
const Hungarian::Matrix EXAMPLE2 = {
  {100, 1},
  {100, 12},
  {1, 4},
  {6, 30252}
};

// TODO: this isn't right
const Hungarian::Matrix SOLUTION2 = {
  {0,1,0,0},
  {0,0,1,0},
  {1,0,0,0},
  {0,0,0,1},
};
// clang-format on

bool runExample(const Hungarian::Matrix& cost,
                const Hungarian::Matrix& solution) {
  Hungarian::Result r = Hungarian::Solve(cost, Hungarian::MODE_MINIMIZE_COST);

  cerr << "cost-matrix:";
  Hungarian::PrintMatrix(r.cost);

  if (!r.success) {
    cerr << "Failed to find solution :(" << endl;
    return false;
  }

  cerr << "assignment:";
  Hungarian::PrintMatrix(r.assignment);

  const bool correct = (r.assignment == solution);
  cerr << "Solution correct? " << correct << endl;

  return correct;
}

int main() {
  bool success = true;

  if (!runExample(EXAMPLE1, SOLUTION1)) success = false;
  cerr << "--------------------" << endl;
  if (!runExample(EXAMPLE2, SOLUTION2)) success = false;

  return success ? 0 : 1;
}
