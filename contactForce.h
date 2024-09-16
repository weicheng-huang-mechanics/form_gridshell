#ifndef CONTACTFORCE_H
#define CONTACTFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class contactForce
{
public:
	contactForce(elasticRod &m_rod, timeStepper &m_stepper, double m_stiffTimes);
	~contactForce();

	void computeFc(MatrixXd m_contactNodes, VectorXd m_contactRefLen);
	void computeJc(MatrixXd m_contactNodes, VectorXd m_contactRefLen);

private:

    elasticRod *rod;
    timeStepper *stepper;

    double refLen;
    double edgeLen;
    double epsX;

    double stiffness;

    Vector3d tangent;

    Vector3d currecntNode1;
    Vector3d currecntNode2;
    int local_nv;

    Vector3d Fc;
    Matrix3d Jc;

    Vector3d u;
    Matrix<double,1,3> v;

    Matrix3d Id3;

    int currentDof, currentDof1, currentDof2;
};

#endif
