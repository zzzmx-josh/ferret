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

#include "ComputeStressBase.h"
#include "RotationMatrix.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

class PMP_Ferroelectrics_ComputeStress;

class PMP_Ferroelectrics_ComputeStress:public ComputeStressBase
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_ComputeStress(const InputParameters &parameters);

protected:
    virtual void computeQpStress() override;
    const Real &_strain_gradient;
    unsigned int _ndisp;
    unsigned int _npolar;

    const VariableSecond &_second_grad_couple_u1;
    const VariableSecond &_second_grad_couple_u2;
    const VariableSecond &_second_grad_couple_u3;

    MaterialProperty<Real> &_f_elastic;

    // define derivatives for residual and stiffness matrix
    MaterialProperty<RankTwoTensor> &_stress;
    MaterialProperty<RankFourTensor> &_dstressdeps;
    MaterialProperty<RankThreeTensor> &_dstressddphi;
    MaterialProperty<RankThreeTensor> &_dstressdP;
    MaterialProperty<RankFourTensor> &_dstressddP;

    MaterialProperty<RankThreeTensor> &_dDddu;
    MaterialProperty<RankTwoTensor> &_dDddphi;
    MaterialProperty<RankTwoTensor> &_dDdP;

    MaterialProperty<RealVectorValue> &_dfdP;
    MaterialProperty<RankTwoTensor> &_dfddP;
    MaterialProperty<RankThreeTensor> &_d2fdPdeps;
    MaterialProperty<RankTwoTensor> &_d2fdPddphi;
    MaterialProperty<RankTwoTensor> &_d2fdPdP;
    MaterialProperty<RankFourTensor> &_d2fddPddP;
    // contributed from Flexoelectricity
    MaterialProperty<RankThreeTensor> &_deps;
    MaterialProperty<RankFourTensor> &_d2fddPdeps;   // contributed from flexoelectricity, method with reduced order
    MaterialProperty<RankFourTensor> &_d2fdPddeps;   // contributed from flexoelectricity, method with strain gradient
    MaterialProperty<RankThreeTensor> &_t;          // double stress tensor, method with strain gradient
    MaterialProperty<RankFourTensor> &_dtdP;         // method with strain gradient

    // get material properties
    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankTwoTensor> &_total_strain;
    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
    const MaterialProperty<RealVectorValue> &_E;
    const MaterialProperty<RankTwoTensor> &_K;          // permittivity
    const MaterialProperty<RankFourTensor> &_g;         // interface constants
    const MaterialProperty<RankFourTensor> &_flexo;     // flexoelectric constants
    const MaterialProperty<Real> &_flexo_status;
    const MaterialProperty<RealVectorValue> &_vec_P;
    const MaterialProperty<RankTwoTensor> &_GradientP;
    const MaterialProperty<RealVectorValue> &_dpsidP;
    const MaterialProperty<RankTwoTensor> &_d2psidPdP;
    const MaterialProperty<RankThreeTensor> &_deps0dP;
    const MaterialProperty<RankFourTensor> &_d2eps0dPdP;

};
