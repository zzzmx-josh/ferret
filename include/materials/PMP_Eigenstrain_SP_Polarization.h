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

#include "ComputeEigenstrainBase.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

class PMP_Eigenstrain_SP_Polarization;

class PMP_Eigenstrain_SP_Polarization : public ComputeEigenstrainBase
{
public:
    static InputParameters validParams();
    PMP_Eigenstrain_SP_Polarization (const InputParameters & parameters);

protected:
    virtual void computeQpEigenstrain();

    const Real &_P0,&_eps0;
    const VariableValue &_couple_P1, &_couple_P2, &_couple_P3;
    unsigned int _npolar;

    MaterialProperty<RankTwoTensor> &_epsilon0;
    MaterialProperty<RankThreeTensor> &_deps0dP;
    MaterialProperty<RankFourTensor> &_d2eps0dPdP;

};
