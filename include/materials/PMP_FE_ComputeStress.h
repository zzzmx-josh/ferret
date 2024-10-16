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

#include "ComputeStressBase.h"
#include "RotationMatrix.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"
#include "Kernel.h"

class PMP_FE_ComputeStress;

class PMP_FE_ComputeStress:public ComputeStressBase
{
public:
    static InputParameters validParams();
    PMP_FE_ComputeStress(const InputParameters &parameters);

protected:
    virtual void computeQpStress() override;
    unsigned int _ndisp;
    unsigned int _npolar;

    MaterialProperty<RankTwoTensor> &_stress;
    MaterialProperty<RankTwoTensor> &_GradientG;
    MaterialProperty<RankTwoTensor> &_GradientF;
    MaterialProperty<Real> &_f_elastic;

    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankTwoTensor> &_total_strain;
    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
    const MaterialProperty<RankFourTensor> &_g;     // interface constants
    const MaterialProperty<RankFourTensor> &_f;     // flexoelectric constants
    const MaterialProperty<RankTwoTensor> &_GradientP;

};
