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

#include "TimeKernel.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

class PMP_Ferroelectrics_P;

class PMP_Ferroelectrics_P:public TimeKernel
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_P(const InputParameters & parameters);

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

    // Variables for second order derivatives
    const VariableSecond & _second_u;
    const VariableTestSecond & _second_test;
    const VariablePhiSecond & _second_phi;

    const MaterialProperty<RealVectorValue> &_dfdP;
    const MaterialProperty<RankTwoTensor> &_dfddP;
    const MaterialProperty<RankThreeTensor> &_d2fdPdeps;
    const MaterialProperty<RankTwoTensor> &_d2fdPddphi;
    const MaterialProperty<RankTwoTensor> &_d2fdPdP;
    const MaterialProperty<RankFourTensor> &_d2fddPddP;
    const MaterialProperty<Real> &_flexo_status;
    const MaterialProperty<RankFourTensor> &_d2fddPdeps;
    const MaterialProperty<RankFourTensor> &_d2fdPddeps;

};
