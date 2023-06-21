#ifndef COURS_IOHANDLER_H
#define COURS_IOHANDLER_H

#include <iosfwd>

namespace print {
  std::ostream& printSustainable(std::ostream& out);

  std::ostream& printUnsustainable(std::ostream& out);

  std::ostream& printBifurcation(std::ostream& out);
};

#endif
