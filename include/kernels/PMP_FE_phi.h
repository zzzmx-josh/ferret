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
#include "ComputeStressBase.h"

class PMP_FE_phi;

class PMP_FE_phi:public Kernel
{
public:
    static InputParameters validParams();
    PMP_FE_phi(const InputParameters & parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    unsigned int _npolar;
    std::vector<unsigned int> _P_var;

    const MaterialProperty<RankTwoTensor> &_mechanical_strain;
    const MaterialProperty<RankTwoTensor> &_K;    // Permitivity
    const MaterialProperty<RealVectorValue> &_D;
    const MaterialProperty<RankThreeTensor> &_deps0dP;

};
