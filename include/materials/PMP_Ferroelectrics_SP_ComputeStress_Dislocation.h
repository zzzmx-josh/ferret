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

class PMP_Ferroelectrics_SP_ComputeStress_Dislocation;

class PMP_Ferroelectrics_SP_ComputeStress_Dislocation:public ComputeStressBase
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_SP_ComputeStress_Dislocation(const InputParameters &parameters);

protected:
    virtual void computeQpStress() override;
    unsigned int _ndisp;
    unsigned int _npolar;

    MaterialProperty<Real> &_f_elastic;

    // define derivatives for residual and stiffness matrix
    MaterialProperty<RankTwoTensor> &_stress;
    MaterialProperty<RankFourTensor> &_dstressdeps;
    MaterialProperty<RankThreeTensor> &_dstressddphi;
    MaterialProperty<RankThreeTensor> &_dstressdP;

    MaterialProperty<RealVectorValue> &_D;
    MaterialProperty<RankThreeTensor> &_dDddu;
    MaterialProperty<RankTwoTensor> &_dDddphi;
    MaterialProperty<RankTwoTensor> &_dDdP;

    MaterialProperty<RealVectorValue> &_dfdP;
    MaterialProperty<RankTwoTensor> &_dfddP;
    MaterialProperty<RankThreeTensor> &_d2fdPdeps;
    MaterialProperty<RankTwoTensor> &_d2fdPddphi;
    MaterialProperty<RankTwoTensor> &_d2fdPdP;
    MaterialProperty<RankFourTensor> &_d2fddPddP;

    // get material properties
    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankTwoTensor> &_total_strain;
    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
    const MaterialProperty<RankThreeTensor> &_e;          // piezoelectric constants
    const MaterialProperty<RealVectorValue> &_E;
    const MaterialProperty<RankTwoTensor> &_K;          // permittivity
    const MaterialProperty<RankFourTensor> &_g;         // gradient energy constants
    const MaterialProperty<RankFourTensor> &_flexo;     // flexoelectric constants
    const MaterialProperty<RankFourTensor> &_dedP;
    const MaterialProperty<RankTwoTensor> &_d2NPdPdP;
    const MaterialProperty<RankThreeTensor> &_e_const;  // piezoelectric constants
    const MaterialProperty<RealVectorValue> &_vec_P;
    const MaterialProperty<RankTwoTensor> &_GradientP;
    const MaterialProperty<RealVectorValue> &_dpsidP;
    const MaterialProperty<RankTwoTensor> &_d2psidPdP;
    const MaterialProperty<RankThreeTensor> &_deps0dP;
    const MaterialProperty<RankFourTensor> &_d2eps0dPdP;
    const MaterialProperty<RankTwoTensor> &_stress_D;
    const MaterialProperty<RankTwoTensor> &_beta_eD;

};
