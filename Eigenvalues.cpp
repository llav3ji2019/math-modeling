#include "Eigenvalues.h"

std::vector<std::complex<double>> ev::getEigenvalues(const std::vector<std::vector<double>> &matrix) {
  std::vector<std::complex<double>> eigenvalues(4);
  for (int i = 0; i < 4; i += 2) {
    double b = -(matrix[i + 1][i + 1] + matrix[i][i]);
    double c = matrix[i][i] * matrix[i + 1][i + 1] - matrix[i][i + 1] * matrix[i + 1][i];
    double discriminant = std::pow(b, 2) - 4 * c;
    if (discriminant >= 0) {
      eigenvalues[i] = std::complex<double>(0.5 * (-b + std::sqrt(discriminant)), 0);
      eigenvalues[i + 1] = std::complex<double>(0.5 * (-b - std::sqrt(discriminant)), 0);
    } else {
      discriminant *= -1;
      eigenvalues[i] = std::complex<double>(0.5 * -b, 0.5 * std::sqrt(discriminant));
      eigenvalues[i + 1] = std::complex<double>(0.5 * -b, -0.5 * std::sqrt(discriminant));
    }
  }
  return eigenvalues;
}

std::ostream& ev::printEigenvalues(std::ostream& out, std::vector<std::complex<double>>& data) {
  int counter = 1;
  for (auto val: data) {
    out << "a" << counter << " =";
    if (val.imag() == 0) {
      out << val.real() << " ";
    } else {
      out << val.real();
      out.setf(std::ios::showpos);
      out << val.imag();
      out.unsetf(std::ios::showpos);
      out << " i ";
    }
    counter++;
  }
  return out;
}

fl::Flags ev::checkFlags(const std::vector<std::complex<double>>& eigenvalues) {
  fl::Flags eigenvaluesFlags;
  for (auto val: eigenvalues) {
    if (std::abs(val.real()) == 0) {
      eigenvaluesFlags.isBifurcation = true;
      break;
    }
    if (val.real() > 0) {
      eigenvaluesFlags.isSustainable = false;
    }
  }
  return eigenvaluesFlags;
}
