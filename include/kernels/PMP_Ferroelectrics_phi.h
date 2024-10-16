/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 1.1          */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#pragma once

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

class PMP_Ferroelectrics_phi;

class PMP_Ferroelectrics_phi:public Kernel
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_phi(const InputParameters & parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    unsigned int _ndisp;
    std::vector<unsigned int> _u_var;
    const unsigned int _varphi_var;
    unsigned int _npolar;
    std::vector<unsigned int> _P_var;

    const MaterialProperty<RealVectorValue> &_D;
    const MaterialProperty<RankThreeTensor> &_dDddu;
    const MaterialProperty<RankTwoTensor> &_dDddphi;
    const MaterialProperty<RankTwoTensor> &_dDdP;

};
