/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 2.0.0        */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#pragma once

#include "Material.h"
#include "RotationMatrix.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"
#include "Kernel.h"

class PMP_Ferroelectrics_SP_Material;

class PMP_Ferroelectrics_SP_Material:public Material
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_SP_Material(const InputParameters &parameters);

protected:
    virtual void computeQpProperties() override;

    const Real &_P0,&_eps0;
    const Real &_a0,&_a1,&_a11,&_a12,&_a111,&_a112,&_a123,&_a1111,&_a1112,&_a1122,&_a1123;
    const Real &_theta;
    const std::vector<Real> &_PermittivityCoe;
    const std::vector<Real> &_PiezoelectricityCoe;
    const std::vector<Real> &_FlexoelectricityCoe;
    const std::vector<Real> &_GradientEnergyCoe;
    unsigned int _npolar;

    const VariableValue &_couple_P1, &_couple_P2, &_couple_P3;
    const VariableGradient &_grad_couple_varphi, &_grad_couple_P1, &_grad_couple_P2, &_grad_couple_P3;

    const VariableValue &_E00,&_E01,&_E02;

    MaterialProperty<RankThreeTensor> &_e_const;  // piezoelectric constants
    MaterialProperty<RankThreeTensor> &_e;
    MaterialProperty<RealVectorValue> &_E;        // electric field
    MaterialProperty<RankTwoTensor> &_K;          // permittivity
    MaterialProperty<RankFourTensor> &_g;         // gradient energy
    MaterialProperty<RankFourTensor> &_flexo;     // flexoelectricity
    MaterialProperty<RankFourTensor> &_dedP;
    MaterialProperty<RankTwoTensor> &_d2NPdPdP;
    MaterialProperty<RealVectorValue> &_vec_P;
    MaterialProperty<RankTwoTensor> &_GradientP;
    MaterialProperty<RealVectorValue> &_dpsidP;
    MaterialProperty<RankTwoTensor> &_d2psidPdP;
    MaterialProperty<Real> &_f_bulk;
    MaterialProperty<Real> &_f_electric;
    MaterialProperty<Real> &_f_gradient;

    RankTwoTensor _K_Rot;
    RankThreeTensor _PiezoCoe;
    RankThreeTensor _PiezoCoe_Rot;
    RankFourTensor _gCoe;
    RankFourTensor _FlexoCoe;
    RankFourTensor _FlexoCoe_Rot;
    RankTwoTensor RotG2L;
    RealVectorValue vec_P;
    RealVectorValue dNormPdP;
    RankTwoTensor d2NormPdPdP;
};
