
#pragma once

#include "Kernel.h"

/**
 * Computes the residual and Jacobian contribution for the weak form
 * of the BiharmonicPolar equation:
 *
 * \int Laplacian(u) * Laplacian(v) dx
 */
class BiharmonicPolar : public Kernel
{
public:
  static InputParameters validParams();

  BiharmonicPolar(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const VariableSecond & _second_u;
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  // const unsigned int _dim;
  const MaterialProperty<Real> & _lambda110;
  const MaterialProperty<Real> & _lambda11;
};
