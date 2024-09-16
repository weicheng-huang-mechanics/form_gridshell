#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"
#include "contactForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	void updateEachRod(int n);
	int simulationRunning();
	int getnumRods();
	int numPoints(int n);
	double getScaledCoordinate(int n, int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);
		
private:

	// Physical parameters
	double rodRadius;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	int numRods;
	double deltaLen;
	double stiffTimes;
	double speed;
    
	double tol, stol;
	int maxIter;
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	VectorXd theta;
	
	// Rod
	std::vector<elasticRod*> rodsVector;
	
	// set up the time stepper
	std::vector<timeStepper*> stepperVector;
	std::vector<double*> totalForceVector;
	double currentTime;
	
	// declare the forces
	std::vector<elasticStretchingForce*> v_stretchForce;
	std::vector<elasticBendingForce*> v_bendingForce;
	std::vector<elasticTwistingForce*> v_twistingForce;
	std::vector<inertialForce*> v_inertialForce;
	std::vector<externalGravityForce*> v_gravityForce;
	std::vector<contactForce*> v_contactForce;

	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry(int i);
	void rodBoundaryCondition(int i);
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	MatrixXd rodStart;
	MatrixXd rodEnd;

	double xStart;
	double yStart;
	double xEnd;
	double yEnd;

	double rodLength;
	int ne_local;
	int nv_local;
	Vector3d direction;
	double delta_local;

	double ComputeDistance(double x1, double y1, double x2, double y2);
	int getIntersection(double p0_x, double p0_y, double p1_x, double p1_y, 
    double p2_x, double p2_y, double p3_x, double p3_y, double &i_x, double &i_y);

    void computeContactPair();

    MatrixXd contactPair;
    int contact_num;

    MatrixXi contactMapping;

    void prepareContact(int n);

    MatrixXd contactNodes;

    void computeContactRefLen();

    VectorXd contactRefLen;

    Vector3d xCurrent;
	Vector3d target;
	Vector3d towards;
};

#endif
