/*************************************************************
* @author Triston Ruiseco
* @file LT.h
* @date 03/22/2021
* @brief Loop tracker class to track the completion of
         nested loops and output completion iteratively.
*************************************************************/

#ifndef LT_H
#define LT_H

#include <iostream>
#include <string>

class LT
{
public:

  /**
   * @pre none
   * @param d: the total number of steps in the loop
   * @param m: tracking message in output
   * @post none
   */
  LT(double d, std::string str);

  /**
   * @pre none
   * @post m_n iterated
   * @brief outputs current progress to console
   */
  void track();

private:
  int m_n;
  double m_p;
  double m_d;
  std::string m_m;
};
#endif
