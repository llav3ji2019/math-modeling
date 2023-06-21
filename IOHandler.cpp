#include <string>
#include "IOHandler.h"

namespace {
  const std::string IS_SUSTAINABLE_LABEL = " yes";
  const std::string IS_UNSUSTAINABLE_LABEL = " no";
  const std::string IS_BIFURCATION_LABEL = " buf";
}

std::ostream& print::printSustainable(std::ostream& out) {
  out << IS_SUSTAINABLE_LABEL;
  return out;
}

std::ostream& print::printUnsustainable(std::ostream& out) {
  out << IS_UNSUSTAINABLE_LABEL;
  return out;
}

std::ostream& print::printBifurcation(std::ostream& out) {
  out << IS_BIFURCATION_LABEL;
  return out;
}
