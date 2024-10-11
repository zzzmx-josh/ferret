#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
// #include "ComputeSpontaneousEigenStrain.h"
#include "ComputeEigenstrainBase.h"
/**
 * ComputeSpontaneousEigenStrain the base class for computing spontaneous polar strain contributions (cubic)
 */
class ComputeSpontaneousEigenStrain : public ComputeEigenstrainBase
{
public:
    ComputeSpontaneousEigenStrain(const InputParameters & parameters);
  static InputParameters validParams();
  void computeQpEigenstrain();
private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;

  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q44;
  std::vector<Real> _vals;
  RankTwoTensor _polar_strain;
};
