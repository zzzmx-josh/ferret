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

class PMP_Ferroelectrics_u;

class PMP_Ferroelectrics_u:public Kernel
{
public:
    static InputParameters validParams();
    PMP_Ferroelectrics_u(const InputParameters &parameters);

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

    // Variables for second order derivatives
    const VariableSecond & _second_u;
    const VariableTestSecond & _second_test;
    const VariablePhiSecond & _second_phi;

    const MaterialProperty<RankTwoTensor> &_stress;
    const MaterialProperty<RankFourTensor> &_dstressdeps;
    const MaterialProperty<RankThreeTensor> &_dstressddphi;
    const MaterialProperty<RankThreeTensor> &_dstressdP;
    const MaterialProperty<RankFourTensor> &_dstressddP;
    const MaterialProperty<Real> &_flexo_status;
    const MaterialProperty<RankThreeTensor> &_t;
    const MaterialProperty<RankFourTensor> &_dtdP;

};
