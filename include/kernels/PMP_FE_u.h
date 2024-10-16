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

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

class PMP_FE_u;

class PMP_FE_u:public Kernel
{
public:
    static InputParameters validParams();
    PMP_FE_u(const InputParameters &parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    unsigned int _component;

    unsigned int _ndisp;
    std::vector<unsigned int> _u_var;
    const unsigned int _varphi_var;
    unsigned int _npolar;
    std::vector<unsigned int> _P_var;

    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
    const MaterialProperty<RankTwoTensor> &_stress;
    const MaterialProperty<RealVectorValue> &_E;
    const MaterialProperty<RankThreeTensor> &_deps0dP;
    const MaterialProperty<RankFourTensor> &_f;
    const MaterialProperty<RankTwoTensor> &_GradientP;

};

