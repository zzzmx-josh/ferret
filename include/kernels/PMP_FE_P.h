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

#include "TimeKernel.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"
#include "ComputeStressBase.h"

class PMP_FE_P;

class PMP_FE_P:public TimeKernel
{
public:
    static InputParameters validParams();
    PMP_FE_P(const InputParameters & parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    const Real &_inv_alpha;
    unsigned int _component;

    unsigned int _ndisp;
    std::vector<unsigned int> _u_var;
    const unsigned int _varphi_var;
    unsigned int _npolar;
    std::vector<unsigned int> _P_var;

    const MaterialProperty<RankTwoTensor> &_K;
    const MaterialProperty<RankFourTensor> &_g;
    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankFourTensor> &_elasticity_tensor;
    const MaterialProperty<RankTwoTensor> &_stress;
    const MaterialProperty<RealVectorValue> &_E;
    const MaterialProperty<RankTwoTensor> &_GradientP;
    const MaterialProperty<RankTwoTensor> &_GradientG;
    const MaterialProperty<RankTwoTensor> &_GradientF;
    const MaterialProperty<RankThreeTensor> &_deps0dP;
    const MaterialProperty<RankFourTensor> &_d2eps0dPdP;
    const MaterialProperty<RealVectorValue> &_dpsidP;
    const MaterialProperty<RankTwoTensor> &_d2psidPdP;
    const MaterialProperty<RankFourTensor> &_f;

    RankTwoTensor I;

};
