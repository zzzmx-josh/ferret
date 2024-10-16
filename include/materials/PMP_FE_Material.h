/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#pragma once

#include "Material.h"
#include "RotationMatrix.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"
#include "Kernel.h"

class PMP_FE_Material;

class PMP_FE_Material:public Material
{
public:
    static InputParameters validParams();
    PMP_FE_Material(const InputParameters &parameters);

protected:
    virtual void computeQpProperties() override;

    const Real &_a0,&_a1,&_a11,&_a12,&_a111,&_a112,&_a123,&_a1111,&_a1112,&_a1122,&_a1123;
    const Real &_theta;
    const std::vector<Real> &_Permittivity;
    const std::vector<Real> &_ElecotrictConstants;    // Electrostrictive constants Q11 Q12 Q44
    const std::vector<Real> &_FlexoConstants;
    const std::vector<Real> &_InterfaceConstants;
    unsigned int _npolar;

    const VariableValue &_couple_P1, &_couple_P2, &_couple_P3;
    const VariableGradient &_grad_couple_varphi, &_grad_couple_P1, &_grad_couple_P2, &_grad_couple_P3;

    const VariableValue &_E00,&_E01,&_E02;

    MaterialProperty<RankTwoTensor> &_K;         // dielectric constants
    MaterialProperty<RankFourTensor> &_Q;        // electrostrictive constants
    MaterialProperty<RankFourTensor> &_g;        // interface constants
    MaterialProperty<RankFourTensor> &_f;        // flexoelectricity constants
    MaterialProperty<RealVectorValue> &_E;
    MaterialProperty<RealVectorValue> &_D;
    MaterialProperty<RankTwoTensor> &_GradientP;
    MaterialProperty<RealVectorValue> &_dpsidP;
    MaterialProperty<RankTwoTensor> &_d2psidPdP;
    MaterialProperty<Real> &_f_bulk;
    MaterialProperty<Real> &_f_electric;
    MaterialProperty<Real> &_f_gradient;

    RankTwoTensor _KCoe;
    RankFourTensor _QCoe;
    RankFourTensor _gCoe;
    RankFourTensor _FlexoCoe;
};
