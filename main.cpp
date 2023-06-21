#include <iostream>
#include <cmath>
#include "zeroin.c"
#include "IOHandler.h"
#include "Flags.h"
#include "Eigenvalues.h"
#include <vector>
#include <fstream>
#include <complex>

namespace {
  const double X2_STEP = 0.1;
  const double X4_STEP = 0.01;
  const double X2_BEGIN = 0;
  const double X2_END = 8;
  const double X4_BEGIN = 0.01;
  const double X4_END = 8;
  const double EPS = 0.000001;
  const int P4 = 2;
  const int P8 = 1000;
  const int P9 = 22;
  const std::string X1_TABLE_PATH = "../x1_table.txt";
  const std::string X2_TABLE_PATH = "../x2_table.txt";
  const std::string X3_TABLE_PATH = "../x3_table.txt";
  const std::string X4_TABLE_PATH = "../x4_table.txt";
  const std::string COMMON_TABLE_PATH = "../common_table.txt";

  double x2 = 0;

  double getX1() {
    return x2 / P9 * (P4 + 1);
  }

  double getP2() {
    return x2 * (1 + P4) / (P9 * (1 - getX1()) * std::exp(x2 / (1 + x2 / P8)));
  }

  double getP3() {
    return getP2();
  }

  double getX3(double x4) {
    return (3 * x4 + 2 * x2) / 22;
  }

  double funcForZeroin(double x4) {
    return getX1() - getX3(x4) + getP3() * (1 - getX3(x4)) * std::exp(x4 / (1 + x4 / P8));
  }

  std::vector<std::vector<double>> getJacobiMatrix(double p2, double x1, double x3, double x4) {
    std::vector<std::vector<double>> jacobiMatrix;
    jacobiMatrix.resize(4, std::vector<double>(4));
    jacobiMatrix[0][0] = -p2 * std::exp(x2 / (1 + x2 / 1000)) - 1;
    jacobiMatrix[0][1] = p2 * (1 - x1) * std::pow(1000 / (1000 + x2), 2) * std::exp(x2 / (1 + x2 / 1000));
    jacobiMatrix[0][2] = 0;
    jacobiMatrix[0][3] = 0;
    jacobiMatrix[1][0] = -22 * p2 * std::exp(x2 / (1 + x2 / 1000));
    jacobiMatrix[1][1] = 22 * p2 * (1 - x1) * std::pow(1000 / (1000 + x2), 2) * std::exp(x2 / (1 + x2 / 1000)) - 3;
    jacobiMatrix[1][2] = 0;
    jacobiMatrix[1][3] = 0;
    jacobiMatrix[2][0] = 1;
    jacobiMatrix[2][1] = 0;
    jacobiMatrix[2][2] = -p2 * std::exp(x4 / (1 + x4 / 1000)) - 1;
    jacobiMatrix[2][3] = p2 * (1 - x3) * std::pow(1000 / (1000 + x4), 2) * std::exp(x4 / (1 + x4 / 1000));
    jacobiMatrix[3][0] = 0;
    jacobiMatrix[3][1] = 1;
    jacobiMatrix[3][2] = -22 * p2 * std::exp(x4 / (1 + x4 / P8));
    jacobiMatrix[3][3] = 22 * p2 * (1 - x3) * std::pow(1000 / (1000 + x4), 2) * std::exp(x4 / (1 + x4 / 1000)) - 3;
    return jacobiMatrix;
  }

  std::ostream& printParameters(std::ostream& out, double x1, double p2, double x3, double x4) {
    return out << "x2 = " << x2 << " x1 = " << x1 << " p2 = " << p2 << " x3 = " << x3 << " x4 = " << x4;
  }
  
  void printDataToTables(std::ostream& out1, std::ostream& out2, std::ostream& out3, std::ostream& out4, double x4, const std::vector<std::complex<double>>& eigenvalues) {
    fl::Flags eigenvaluesFlags = ev::checkFlags(eigenvalues);
    if (getP2() > 0) {
      if (getX1() > 0) {
        out1 << getX1() << " " << getP2();
        if (eigenvaluesFlags.isBifurcation) {
          print::printBifurcation(out1) << "\n";
        } else if (eigenvaluesFlags.isSustainable) {
          print::printSustainable(out1) << "\n";
        } else {
          print::printUnsustainable(out1) << "\n";
        }
      }
      if (getX3(x4) > 0) {
        out3 << getX3(x4) << " " << getP2();
        if (eigenvaluesFlags.isBifurcation) {
          print::printBifurcation(out3) << "\n";
        } else if (eigenvaluesFlags.isSustainable) {
          print::printSustainable(out3) << "\n";
        } else {
          print::printUnsustainable(out3) << "\n";
        }
      }
      out2 << x2 << " " << getP2();
      out4 << x4 << " " << getP2();
      if (eigenvaluesFlags.isBifurcation) {
        print::printBifurcation(out2) << "\n";
        print::printBifurcation(out4) << "\n";
      } else if (eigenvaluesFlags.isSustainable) {
        print::printSustainable(out2) << "\n";
        print::printSustainable(out4) << "\n";
      } else {
        print::printUnsustainable(out2) << "\n";
        print::printUnsustainable(out4) << "\n";
      }
    }
  }

  std::ostream& printParametersToTable(std::ostream& out, double x4) {
    out << getX1() << " " << x2 << " " << getX3(x4) << " " << x4 << " " << getP2();
    return out;
  }
}

int main() {
  std::ofstream x1TableOutput(X1_TABLE_PATH);
  std::ofstream x2TableOutput(X2_TABLE_PATH);
  std::ofstream x3TableOutput(X3_TABLE_PATH);
  std::ofstream x4TableOutput(X4_TABLE_PATH);
  std::ofstream commonTableOutput(COMMON_TABLE_PATH);
  for (x2 = X2_BEGIN; x2 <= X2_END; x2 += X2_STEP) {
    for (double x4 = X4_BEGIN; x4 <= X4_END; x4 += X4_STEP) {
      int flag;
      if (funcForZeroin(x4 - X4_STEP) * funcForZeroin(x4) < 0) {
        double currentX4 = zeroin(x4 - X4_STEP, x4, funcForZeroin, EPS , &flag);
        printParameters(std::cout, getX1(), getP2(), getX3(currentX4), currentX4) << "\n";
        printParametersToTable(commonTableOutput, currentX4) << "\n";
        std::vector<std::vector<double>> jacobiMatrix = getJacobiMatrix(getP2(), getX1(), getX3(currentX4), currentX4);
        std::vector<std::complex<double>> eigenvalues = ev::getEigenvalues(jacobiMatrix);
        ev::printEigenvalues(std::cout, eigenvalues) << "\n\n";
        printDataToTables(x1TableOutput, x2TableOutput, x3TableOutput, x4TableOutput, currentX4, eigenvalues);
      }
    }
  }
  x1TableOutput.close();
  x2TableOutput.close();
  x3TableOutput.close();
  x4TableOutput.close();
  commonTableOutput.close();
  return 0;
}
