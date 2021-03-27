/*************************************************************
* @author Triston Ruiseco
* @file LT.cpp
* @date 03/22/2021
* @brief Implementation file for LT.h.
*************************************************************/

#include "LT.h"

LT::LT(double d, std::string str){
  std::cout << "\n";
  m_n = 0;
  m_p = 0.0;
  m_d = d;
  m_d /= 100;
  m_m = str;
}

void LT::track(){
  m_n++;
  m_p = double(m_n)/m_d;
  if(int(m_p) == m_p){
    std::cout << "\r" << m_m << int(m_p) << "%" << std::flush;
  }
}
