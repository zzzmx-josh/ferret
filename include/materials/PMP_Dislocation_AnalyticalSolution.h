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

class PMP_Dislocation_AnalyticalSolution;

class PMP_Dislocation_AnalyticalSolution:public Material
{
public:
    static InputParameters validParams();
    PMP_Dislocation_AnalyticalSolution(const InputParameters &parameters);

protected:
    virtual void computeQpProperties() override;

    const std::vector<Real> &_b_n;
    const std::vector<Real> &_r0;
    const Real &_h;
    const Real &_angle;

    unsigned int _ndisp;

    MaterialProperty<RankTwoTensor> &_stress_D;
    MaterialProperty<RankTwoTensor> &_beta_eD;

    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
};
