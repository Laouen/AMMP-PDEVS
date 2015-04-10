#ifndef BOOST_SIMULATION_PDEVS_PARSER_ITERATOR_H
#define BOOST_SIMULATION_PDEVS_PARSER_ITERATOR_H

#include <iostream>
#include <iterator>
#include "parser.h"

using namespace std;
namespace Parser {

class iterator : public std::iterator<bidirectional_iterator_tag, TiXmlElement>{

private:
  TiXmlElement* p;

public:
  iterator(TiXmlElement* x) :p(x) {}
  iterator(const iterator& mit) 
  : p(mit.p) {}
  iterator& operator++() {
  	++p;
  	return *this;
  }
  iterator operator++(int) {
  	iterator tmp(*this); 
  	operator++(); 
  	return tmp;
  }
  bool operator==(const iterator& rhs) {
  	return p==rhs.p;
  }
  bool operator!=(const iterator& rhs) {
  	return p!=rhs.p;
  }
  int& operator*() {
  	return *p;
  }
};

}

#endif // BOOST_SIMULATION_PDEVS_PARSER_ITERATOR_H