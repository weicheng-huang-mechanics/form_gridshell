#include "timeStepper.h"

timeStepper::timeStepper(elasticRod &m_rod)
{
	rod = &m_rod;
	kl = 10; // lower diagonals
	ku = 10; // upper diagonals
	freeDOF = rod->uncons;
	ldb = freeDOF;
	NUMROWS = 2 * kl + ku + 1;
	totalForce = new double[freeDOF];
	jacobianLen = (2 * kl + ku + 1) * freeDOF;
	jacobian = new double [jacobianLen];
	nrhs = 1;
    ipiv = new int[freeDOF];
    info = 0;

    /*

    ForceVector.setZero(freeDOF, 1);
    VectorMotion.setZero(freeDOF, 1);

    JacobianMatrix.resize(freeDOF, freeDOF);

    */
}

timeStepper::~timeStepper()
{
	;
}

double* timeStepper::getForce()
{
	return totalForce;
}

double* timeStepper::getJacobian()
{
	return jacobian;
}

void timeStepper::addForce(int ind, double p)
{
	if (rod->getIfConstrained(ind) == 0) // free dof
	{
		mappedInd = rod->fullToUnconsMap[ind];
		totalForce[mappedInd] = totalForce[mappedInd] + p; // subtracting elastic force

		//ForceVector(mappedInd) = ForceVector(mappedInd) + p;
	}
}

void timeStepper::addJacobian(int ind1, int ind2, double p)
{
	mappedInd1 = rod->fullToUnconsMap[ind1];
	mappedInd2 = rod->fullToUnconsMap[ind2];
	if (rod->getIfConstrained(ind1) == 0 && rod->getIfConstrained(ind2) == 0 && p != 0) // both are free
	{
		row = kl + ku + mappedInd2 - mappedInd1;
        col = mappedInd1;
        offset = row + col * NUMROWS;
        jacobian[offset] = jacobian[offset] + p;

        //JacobianMatrix.coeffRef(mappedInd2, mappedInd1) += p;
	}
}

void timeStepper::setZero()
{
	for (int i=0; i < freeDOF; i++)
		totalForce[i] = 0;
	for (int i=0; i < jacobianLen; i++)
		jacobian[i] = 0;

	/*

	ForceVector.setZero(freeDOF, 1);
    VectorMotion.setZero(freeDOF, 1);

    JacobianMatrix.setZero();

    */
}

void timeStepper::integrator()
{

	dgbsv_(&freeDOF, &kl, &ku, &nrhs, jacobian, &NUMROWS, ipiv, totalForce, &ldb, &info);


	/*

	ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(JacobianMatrix);
    VectorMotion = cg.solve(ForceVector);

    for (int i=0; i < freeDOF; i++)
    {
		totalForce[i] = VectorMotion(i);
    }

    */
}
