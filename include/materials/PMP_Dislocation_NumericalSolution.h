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

class PMP_Dislocation_NumericalSolution;

class PMP_Dislocation_NumericalSolution:public Material
{
public:
    static InputParameters validParams();
    PMP_Dislocation_NumericalSolution(const InputParameters &parameters);

protected:
    virtual void computeQpProperties() override;

    const VariableValue &_sigma_D_11, &_sigma_D_22, &_sigma_D_33;
    const VariableValue &_sigma_D_12, &_sigma_D_13, &_sigma_D_23;

    const VariableValue &_beta_eD_11, &_beta_eD_12, &_beta_eD_13;
    const VariableValue &_beta_eD_21, &_beta_eD_22, &_beta_eD_23;
    const VariableValue &_beta_eD_31, &_beta_eD_32, &_beta_eD_33;

    MaterialProperty<RankTwoTensor> &_stress_D;
    MaterialProperty<RankTwoTensor> &_beta_eD;

};
