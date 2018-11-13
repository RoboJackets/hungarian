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

#pragma once

#include <iostream>
#include <stdio.h>
#include <vector>

typedef enum {
  HUNGARIAN_MODE_MINIMIZE_COST,
  HUNGARIAN_MODE_MAXIMIZE_UTIL,
} MODE;

typedef enum {
  HUNGARIAN_NOT_ASSIGNED,
  HUNGARIAN_ASSIGNED,
} ASSIGN;

class Hungarian {

public:
  /** This method initialize the hungarian_problem structure and init
   *  the  cost matrices (missing lines or columns are filled with 0).
   *  It returns the size of the quadratic(!) assignment matrix. **/

  Hungarian();
  Hungarian(const std::vector<std::vector<int>> &, int, int, MODE);

  int init(const std::vector<std::vector<int>> &input_matrix, int rows,
           int cols, MODE mode);

  /** This method computes the optimal assignment. **/
  bool solve();

  /** Accessor for the cost **/
  int cost() const;

  /** Reference accessor for assignment **/
  const std::vector<std::vector<int>> &assignment() const;

  /** Print the computed optimal assignment. **/
  void print_assignment();

  /** Print the cost matrix. **/
  void print_cost();

  /** Print cost matrix and assignment matrix. **/
  void print_status();

protected:
  bool check_solution(const std::vector<int> &row_dec,
                      const std::vector<int> &col_inc,
                      const std::vector<int> &col_vertex);
  bool assign_solution(const std::vector<int> &row_dec,
                       const std::vector<int> &col_inc,
                       const std::vector<int> &col_vertex);

private:
  int m_cost;
  int m_rows;
  int m_cols;
  std::vector<std::vector<int>> m_costmatrix;
  std::vector<std::vector<int>> m_assignment;
};
