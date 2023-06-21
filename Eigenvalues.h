#ifndef COURS_EIGENVALUES_H
#define COURS_EIGENVALUES_H

#include <vector>
#include <complex>
#include "Flags.h"

namespace ev {
  std::vector<std::complex<double>> getEigenvalues(const std::vector<std::vector<double>> &matrix);

  std::ostream& printEigenvalues(std::ostream& out, std::vector<std::complex<double>>& data);

  fl::Flags checkFlags(const std::vector<std::complex<double>>& eigenvalues);
}

#endif
