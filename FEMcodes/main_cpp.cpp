/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Interactive deformable object simulator" driver application,         *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin, Daniel Schroeder          *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
//#include <windows.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <float.h>
using namespace std;

#ifdef WIN32
  #include <windows.h>
#endif

#ifdef __APPLE__
  #include "TargetConditionals.h"
#endif

#include "getopts.h"
//#include "initGraphics.h"
//#include "sceneObjectDeformable.h"
#include "performanceCounter.h"
#include "tetMesh.h"
#include "StVKCubeABCD.h"
#include "StVKTetABCD.h"
#include "StVKTetHighMemoryABCD.h"
#include "implicitBackwardEulerSparse.h"
#include "eulerSparse.h"
//#include "centralDifferencesSparse.h"
#include "StVKInternalForces.h"
#include "StVKStiffnessMatrix.h"
//#include "StVKInternalForcesMT.h"
//#include "StVKStiffnessMatrixMT.h"
#include "StVKForceModel.h"
#include "massSpringSystemForceModel.h"
#include "corotationalLinearFEM.h"
//#include "corotationalLinearFEMMT.h"
#include "corotationalLinearFEMForceModel.h"
#include "linearFEMForceModel.h"
#include "isotropicHyperelasticFEM.h"
//#include "isotropicHyperelasticFEMMT.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "isotropicMaterial.h"
#include "StVKIsotropicMaterial.h"
#include "neoHookeanIsotropicMaterial.h"
#include "MooneyRivlinIsotropicMaterial.h"
#include "getIntegratorSolver.h"
#include "volumetricMeshLoader.h"
#include "StVKElementABCDLoader.h"
#include "generateMeshGraph.h"
#include "generateMassMatrix.h"
#include "massSpringSystem.h"
//#include "massSpringSystemMT.h"
#include "massSpringSystemFromObjMeshConfigFile.h"
#include "massSpringSystemFromTetMeshConfigFile.h"
#include "massSpringSystemFromCubicMeshConfigFile.h"
//#include "graph.h"
//#include "renderSprings.h"
#include "configFile.h"
//#include "GL/glui.h"
//#include "lighting.h"
#include "loadList.h"
#include "matrixIO.h"

// graphics
char windowTitleBase[4096] = "Real-time sim";
void displayFunction(void);
int windowID;
int windowWidth = 800;
int windowHeight = 600;
double zNear=0.01;
double zFar=10.0;
double cameraRadius;
double focusPositionX, focusPositionY, focusPositionZ;
double cameraLongitude, cameraLattitude;
//SphericalCamera * camera = NULL;
int g_iLeftMouseButton=0, g_iMiddleMouseButton=0, g_iRightMouseButton=0;
int g_vMousePos[2] = {0,0};
int shiftPressed=0;
int altPressed=0;
int ctrlPressed=0;
int renderWireframe=1;
int renderAxes=0;
int renderDeformableObject=1;
int renderSecondaryDeformableObject=1;
int useRealTimeNormals = 0;
int renderGroundPlane = 0;
int renderFixedVertices = 1;
int renderSprings = 0;
int renderVertices = 0;
int lockScene=0;
int pauseSimulation=0;
int singleStepMode=0;
char groundPlaneString[128];
double groundPlaneHeight;
double groundPlaneLightHeight = 10.0;
double groundPlaneSize = 15.0;
// config file
string configFilename;
char renderingMeshFilename[4096];
char secondaryRenderingMeshFilename[4096];
char secondaryRenderingMeshInterpolationFilename[4096];
char volumetricMeshFilename[4096];
char customMassSpringSystem[4096];
char deformableObjectMethod[4096];
char fixedVerticesFilename[4096];
char massSpringSystemObjConfigFilename[4096];
char massSpringSystemTetMeshConfigFilename[4096];
char massSpringSystemCubicMeshConfigFilename[4096];
char invertibleMaterialString[4096] = "__none";
char initialPositionFilename[4096];
char initialVelocityFilename[4096];
char forceLoadsFilename[4096];
char outputFilename[4096];
int corotationalLinearFEM_warp = 1;
const int max_corotationalLinearFEM_warp = 2;
char implicitSolverMethod[4096];
char solverMethod[4096];
char extraSceneGeometryFilename[4096];
char lightingConfigFilename[4096];
float dampingMassCoef; // Rayleigh mass damping
float dampingStiffnessCoef; // Rayleigh stiffness damping
float dampingLaplacianCoef = 0.0; // Laplacian damping (rarely used)
float deformableObjectCompliance = 1.0; // scales all user forces by the provided factor
// adjusts the stiffness of the object to cause all frequencies scale by the provided factor:
// keep it to 1.0 (except for experts)
float frequencyScaling = 1.0;
int maxIterations; // for implicit integration
double epsilon; // for implicit integration
char backgroundColorString[4096] = "255 255 255";
int numInternalForceThreads;
int numSolverThreads;
// simulation
int syncTimestepWithGraphics=1;
/**Important Parameters**************/
float timeStep = 1.0 / 30;
int timestepCounter = 1,numOutputSteps,bodyType; // bodyType=3 for solid body, =2 for membrane
double heavingAmplitude,heavingFrequency,forceRatio,*flowGrid[3];
double *bodyNormalVectorCW,*markerForce,markerArea,flowRe;
int nFlowGrid[3],*surfaceTriangleMesh,nTriangleMesh;
/************************************/
int subTimestepCounter = 0;
float newmarkBeta = 0.25;
float newmarkGamma = 0.5;
int use1DNewmarkParameterFamily = 1;
int substepsPerTimeStep = 1;
double inversionThreshold;
double fps = 0.0;
const int fpsBufferSize = 5;
int fpsHead = 0;
double fpsBuffer[fpsBufferSize];
double cpuLoad = 0;
double forceAssemblyTime = 0.0;
double forceAssemblyLocalTime = 0.0;
const int forceAssemblyBufferSize = 50;
int forceAssemblyHead = 0;
double forceAssemblyBuffer[forceAssemblyBufferSize];
double systemSolveTime = 0.0;
double systemSolveLocalTime = 0.0;
const int systemSolveBufferSize = 50;
int systemSolveHead = 0;
double systemSolveBuffer[systemSolveBufferSize];
int enableTextures = 0;
int staticSolver = 0;
int graphicFrame = 0;
int lockAt30Hz = 0;
int pulledVertex = -1;
int forceNeighborhoodSize = 5;
int dragStartX, dragStartY;
int explosionFlag = 0;
PerformanceCounter titleBarCounter;
PerformanceCounter explosionCounter;
PerformanceCounter cpuLoadCounter;
int numFixedVertices;
int * fixedVertices;
int numForceLoads = 0;
double * forceLoads = NULL;
IntegratorBase * integratorBase = NULL;
ImplicitNewmarkSparse * implicitNewmarkSparse = NULL;
IntegratorBaseSparse * integratorBaseSparse = NULL;
ForceModel * forceModel = NULL;
StVKInternalForces * stVKInternalForces = NULL;
StVKStiffnessMatrix * stVKStiffnessMatrix = NULL;
StVKForceModel * stVKForceModel = NULL;
MassSpringSystemForceModel * massSpringSystemForceModel = NULL;
CorotationalLinearFEMForceModel * corotationalLinearFEMForceModel = NULL;
int enableCompressionResistance = 1;
double compressionResistance = 500;
int centralDifferencesTangentialDampingUpdateMode = 1;
int positiveDefinite = 0;
int addGravity=0;
double g=9.81;
VolumetricMesh * volumetricMesh = NULL;
TetMesh * tetMesh = NULL;
Graph * meshGraph = NULL;
enum massSpringSystemSourceType { OBJ, TETMESH, CUBICMESH, CHAIN, NONE } massSpringSystemSource = NONE;
enum deformableObjectType { STVK, COROTLINFEM, LINFEM, MASSSPRING, INVERTIBLEFEM, UNSPECIFIED } deformableObject = UNSPECIFIED;
enum invertibleMaterialType { INV_STVK, INV_NEOHOOKEAN, INV_MOONEYRIVLIN, INV_NONE } invertibleMaterial = INV_NONE;
enum solverType { IMPLICITNEWMARK, IMPLICITBACKWARDEULER, EULER, SYMPLECTICEULER, CENTRALDIFFERENCES, UNKNOWN } solver = UNKNOWN;
MassSpringSystem * massSpringSystem = NULL;
//RenderSprings * renderMassSprings = NULL;
SparseMatrix * massMatrix = NULL;
SparseMatrix * LaplacianDampingMatrix = NULL;
int n; // the number of vertices
double * u = NULL, *u_base=NULL;  // u_base is used to store the u of the last timeStep.
double * uvel = NULL, *uvel_base=NULL;
double * uaccel = NULL, *uaccel_base=NULL;
double * f_ext = NULL;
double * f_extBase = NULL;
double * f_int =NULL, *f_int_base=NULL;
double * uSecondary = NULL;
double * uInitial = NULL;
double * velInitial = NULL;
int * verticesRelation_vega_IB; // vega Vertice, IB vertice,

// set up the configuration file
void initConfigurations()
{
  printf("Parsing configuration file %s...\n", configFilename.c_str());
  ConfigFile configFile;
  // specify the entries of the config file, at least one of the following must be present:
  configFile.addOptionOptional("volumetricMeshFilename", volumetricMeshFilename, "__none");
  configFile.addOptionOptional("customMassSpringSystem", customMassSpringSystem, "__none");
  configFile.addOptionOptional("deformableObjectMethod", deformableObjectMethod, "StVK");
  configFile.addOptionOptional("massSpringSystemObjConfigFilename", massSpringSystemObjConfigFilename, "__none");
  configFile.addOptionOptional("massSpringSystemTetMeshConfigFilename", massSpringSystemTetMeshConfigFilename, "__none");
  configFile.addOptionOptional("massSpringSystemCubicMeshConfigFilename", massSpringSystemCubicMeshConfigFilename, "__none");
  // option for corotational linear FEM: if warp is disabled, one gets purely linear FEM
  configFile.addOptionOptional("corotationalLinearFEM_warp", &corotationalLinearFEM_warp, corotationalLinearFEM_warp);
  configFile.addOptionOptional("implicitSolverMethod", implicitSolverMethod, "none"); // this is now obsolete, but preserved for backward compatibility, use "solver" below
  configFile.addOptionOptional("solver", solverMethod, "implicitNewmark");
  configFile.addOptionOptional("centralDifferencesTangentialDampingUpdateMode", &centralDifferencesTangentialDampingUpdateMode, centralDifferencesTangentialDampingUpdateMode);
  configFile.addOptionOptional("initialPositionFilename", initialPositionFilename, "__none");
  configFile.addOptionOptional("initialVelocityFilename", initialVelocityFilename, "__none");
  configFile.addOptionOptional("outputFilename", outputFilename, "__none");
  configFile.addOptionOptional("addGravity", &addGravity, addGravity);
  configFile.addOptionOptional("g", &g, g);
  configFile.addOptionOptional("renderingMeshFilename", renderingMeshFilename, "__none");
  configFile.addOptionOptional("secondaryRenderingMeshFilename", secondaryRenderingMeshFilename, "__none");
  configFile.addOptionOptional("secondaryRenderingMeshInterpolationFilename", secondaryRenderingMeshInterpolationFilename, "__none");
  configFile.addOptionOptional("useRealTimeNormals", &useRealTimeNormals, 0);
  configFile.addOptionOptional("fixedVerticesFilename", fixedVerticesFilename, "__none");
  configFile.addOptionOptional("enableCompressionResistance", &enableCompressionResistance, enableCompressionResistance);
  configFile.addOptionOptional("compressionResistance", &compressionResistance, compressionResistance);
  configFile.addOption("timestep", &timeStep);
  configFile.addOptionOptional("substepsPerTimeStep", &substepsPerTimeStep, substepsPerTimeStep);
  configFile.addOptionOptional("syncTimestepWithGraphics", &syncTimestepWithGraphics, syncTimestepWithGraphics);
  configFile.addOption("dampingMassCoef", &dampingMassCoef);
  configFile.addOption("dampingStiffnessCoef", &dampingStiffnessCoef);
  configFile.addOptionOptional("dampingLaplacianCoef", &dampingLaplacianCoef, dampingLaplacianCoef);
  configFile.addOptionOptional("newmarkBeta", &newmarkBeta, newmarkBeta);
  configFile.addOptionOptional("newmarkGamma", &newmarkGamma, newmarkGamma);
  configFile.addOption("deformableObjectCompliance", &deformableObjectCompliance);
  configFile.addOption("frequencyScaling", &frequencyScaling);
  configFile.addOptionOptional("forceNeighborhoodSize", &forceNeighborhoodSize, forceNeighborhoodSize);
  configFile.addOptionOptional("maxIterations", &maxIterations, 1);
  configFile.addOptionOptional("epsilon", &epsilon, 1E-6);
  configFile.addOptionOptional("numInternalForceThreads", &numInternalForceThreads, 0);
  configFile.addOptionOptional("numSolverThreads", &numSolverThreads, 1);
  configFile.addOptionOptional("inversionThreshold", &inversionThreshold, -DBL_MAX);
  configFile.addOptionOptional("forceLoadsFilename", forceLoadsFilename, "__none");
  configFile.addOptionOptional("windowWidth", &windowWidth, windowWidth);
  configFile.addOptionOptional("windowHeight", &windowHeight, windowHeight);
  configFile.addOptionOptional("cameraRadius", &cameraRadius, 17.5);
  configFile.addOptionOptional("focusPositionX", &focusPositionX, 0.0);
  configFile.addOptionOptional("focusPositionY", &focusPositionY, 0.0);
  configFile.addOptionOptional("focusPositionZ", &focusPositionZ, 0.0);
  configFile.addOptionOptional("cameraLongitude", &cameraLongitude, -60.0);
  configFile.addOptionOptional("cameraLattitude", &cameraLattitude, 20.0);
  configFile.addOptionOptional("renderWireframe", &renderWireframe, 1);
  configFile.addOptionOptional("renderSecondaryDeformableObject", &renderSecondaryDeformableObject, renderSecondaryDeformableObject);
  configFile.addOptionOptional("renderAxes", &renderAxes, renderAxes);
  configFile.addOptionOptional("extraSceneGeometry", extraSceneGeometryFilename, "__none");
  configFile.addOptionOptional("enableTextures", &enableTextures, enableTextures);
  configFile.addOptionOptional("backgroundColor", backgroundColorString, backgroundColorString);
  configFile.addOption("lightingConfigFilename", lightingConfigFilename);
  configFile.addOptionOptional("groundPlane", groundPlaneString, "__none");
  configFile.addOptionOptional("singleStepMode", &singleStepMode, singleStepMode);
  configFile.addOptionOptional("pauseSimulation", &pauseSimulation, pauseSimulation);
  configFile.addOptionOptional("lockAt30Hz", &lockAt30Hz, lockAt30Hz);
  configFile.addOptionOptional("invertibleMaterial", invertibleMaterialString, invertibleMaterialString);
  configFile.addOption("heavingAmplitude",&heavingAmplitude);  // Add by CJ Yuan
  configFile.addOption("heavingFrequency",&heavingFrequency);  // Add by CJ Yuan
  configFile.addOptionOptional("numOutputSteps",&numOutputSteps,1);  // Add by CJ Yuan
  configFile.addOptionOptional("forceRatio",&forceRatio,100.0); // Add by CJ Yuan
  configFile.addOptionOptional("bodyType",&bodyType,2); // Add by CJ Yuan. default bodyType is membrane.
  configFile.addOption("nFlowGridX",&nFlowGrid[0]); // Add by CJ Yuan
  configFile.addOption("nFlowGridY",&nFlowGrid[1]); // Add by CJ Yuan
  configFile.addOption("nFlowGridZ",&nFlowGrid[2]); // Add by CJ Yuan
  configFile.addOption("markerArea",&markerArea); // Add by CJ Yuan
  configFile.addOptionOptional("flowRe",&flowRe,200.0); // Add by CJ Yuan

  // parse the configuration file
  if (configFile.parseOptions((char*)configFilename.c_str()) != 0)
  {  printf("Error parsing options.\n"); exit(1); }

  // the config variables have now been loaded with their specified values
  // informatively print the variables (with assigned values) that were just parsed
  configFile.printOptions();
  // set the solver based on config file input
  solver = UNKNOWN;
  if (strcmp(implicitSolverMethod, "implicitNewmark") == 0)
    solver = IMPLICITNEWMARK;
  if (strcmp(implicitSolverMethod, "implicitBackwardEuler") == 0)
    solver = IMPLICITBACKWARDEULER;

  if (strcmp(solverMethod, "implicitNewmark") == 0)
    solver = IMPLICITNEWMARK;
  if (strcmp(solverMethod, "implicitBackwardEuler") == 0)
    solver = IMPLICITBACKWARDEULER;
  if (strcmp(solverMethod, "Euler") == 0)
    solver = EULER;
  if (strcmp(solverMethod, "symplecticEuler") == 0)
    solver = SYMPLECTICEULER;
  if (strcmp(solverMethod, "centralDifferences") == 0)
    solver = CENTRALDIFFERENCES;

  if (solver == UNKNOWN)
  {
    printf("Error: unknown implicit solver specified.\n");
    exit(1);
  }
}

// program initialization
void initSimulation()
{
  volumetricMesh = NULL;
  massSpringSystem = NULL;

  // set deformable material type
  if (strcmp(volumetricMeshFilename, "__none") != 0)
  {
    if (strcmp(deformableObjectMethod, "StVK") == 0)
      deformableObject = STVK;
    if (strcmp(deformableObjectMethod, "CLFEM") == 0)
      deformableObject = COROTLINFEM;
    if (strcmp(deformableObjectMethod, "LinearFEM") == 0)
      deformableObject = LINFEM;
    if (strcmp(deformableObjectMethod, "InvertibleFEM") == 0)
      deformableObject = INVERTIBLEFEM;
  }

  if (strcmp(massSpringSystemObjConfigFilename, "__none") != 0)
    massSpringSystemSource = OBJ;
  else if (strcmp(massSpringSystemTetMeshConfigFilename, "__none") != 0)
    massSpringSystemSource = TETMESH;
  else if (strcmp(massSpringSystemCubicMeshConfigFilename, "__none") != 0)
    massSpringSystemSource = CUBICMESH;
  else if (strncmp(customMassSpringSystem, "chain", 5) == 0)
    massSpringSystemSource = CHAIN;

  if ((massSpringSystemSource == OBJ) || (massSpringSystemSource == TETMESH) || (massSpringSystemSource == CUBICMESH) || (massSpringSystemSource == CHAIN))
    deformableObject = MASSSPRING;

  if (deformableObject == UNSPECIFIED)
  { printf("Error: no deformable model specified.\n"); exit(1); }

  // load volumetric mesh
  if ((deformableObject == STVK) || (deformableObject == COROTLINFEM) || (deformableObject == LINFEM) || (deformableObject == INVERTIBLEFEM))
  {
    printf("Loading volumetric mesh from file %s...\n", volumetricMeshFilename);

    VolumetricMesh::fileFormatType fileFormat = VolumetricMesh::ASCII;
    int verbose = 0;
    volumetricMesh = VolumetricMeshLoader::load(volumetricMeshFilename, fileFormat, verbose);
    if (volumetricMesh == NULL)
    {
      printf("Error: unable to load the volumetric mesh from %s.\n", volumetricMeshFilename);
      exit(1);
    }

    n = volumetricMesh->getNumVertices();
    printf("Num vertices: %d. Num elements: %d\n", n, volumetricMesh->getNumElements());
    meshGraph = GenerateMeshGraph::Generate(volumetricMesh);

    // create the mass matrix
    int inflate3Dim = true; // necessary so that the returned matrix is 3n x 3n
    GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, inflate3Dim);

    // create the internal forces for STVK and linear FEM materials
    if (deformableObject == STVK || deformableObject == LINFEM)  // LINFEM is constructed from stVKInternalForces
    {
      unsigned int loadingFlag = 0; // 0 = use the low-memory version, 1 = use the high-memory version
      StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh, loadingFlag);
      if (precomputedIntegrals == NULL)
      {
        printf("Error: unable to load the StVK integrals.\n");
        exit(1);
      }

      printf("Generating internal forces and stiffness matrix models...\n"); fflush(NULL);
      if (numInternalForceThreads == 0)
        stVKInternalForces = new StVKInternalForces(volumetricMesh, precomputedIntegrals, addGravity, g);
      else
      {printf("Invalid numInternalForceThreads %d \n",numInternalForceThreads);exit(0);}
        //stVKInternalForces = new StVKInternalForcesMT(volumetricMesh, precomputedIntegrals, addGravity, g, numInternalForceThreads);

      if (numInternalForceThreads == 0)
        stVKStiffnessMatrix = new StVKStiffnessMatrix(stVKInternalForces);
      else
        {printf("Invalid numInternalForceThreads %d \n",numInternalForceThreads);exit(0);}
        //stVKStiffnessMatrix = new StVKStiffnessMatrixMT(stVKInternalForces, numInternalForceThreads);
    }
  }
  // load mass spring system (if any)
  // codes deleted
  int scaleRows =1;
  meshGraph->GetLaplacian(&LaplacianDampingMatrix, scaleRows);
  LaplacianDampingMatrix->ScalarMultiply(dampingLaplacianCoef);

  if (!((deformableObject == MASSSPRING) && (massSpringSystemSource == CHAIN)))
  {
    // read the fixed vertices
    // 1-indexed notation
    if (strcmp(fixedVerticesFilename, "__none") == 0)
    {
      numFixedVertices = 0;
      fixedVertices = NULL;
    }
    else
    {
      if (LoadList::load(fixedVerticesFilename, &numFixedVertices,&fixedVertices) != 0)
      {
        printf("Error reading fixed vertices.\n");
        exit(1);
      }
      LoadList::sort(numFixedVertices, fixedVertices);
    }
  }
  else
  {
    numFixedVertices = 1;
    fixedVertices = (int*) malloc (sizeof(int) * numFixedVertices);
    fixedVertices[0] = massSpringSystem->GetNumParticles();
  }
  printf("Loaded %d fixed vertices. They are:\n",numFixedVertices);
  LoadList::print(numFixedVertices,fixedVertices);
  // create 0-indexed fixed DOFs
  int numFixedDOFs = 3 * numFixedVertices;
  int * fixedDOFs = (int*) malloc (sizeof(int) * numFixedDOFs);
  for(int i=0; i<numFixedVertices; i++)
  {
    fixedDOFs[3*i+0] = 3*fixedVertices[i]-3;
    fixedDOFs[3*i+1] = 3*fixedVertices[i]-2;
    fixedDOFs[3*i+2] = 3*fixedVertices[i]-1;
  }
  for(int i=0; i<numFixedVertices; i++)
    fixedVertices[i]--;
  printf("Boundary vertices processed.\n");

  // make room for deformation and force vectors
  u = (double*) calloc (3*n, sizeof(double)); u_base = (double*) calloc(3*n, sizeof(double));
  uvel = (double*) calloc (3*n, sizeof(double)); uvel_base = (double*) calloc(3*n, sizeof(double));
  uaccel = (double*) calloc (3*n, sizeof(double)); uaccel_base = (double*) calloc(3*n, sizeof(double));
  f_ext = (double*) calloc (3*n, sizeof(double));
  f_extBase = (double*) calloc (3*n, sizeof(double));
  f_int = (double*) calloc (3*n, sizeof(double)); f_int_base = (double*) calloc(3*n, sizeof(double));
  verticesRelation_vega_IB=(int *)calloc(3*n,sizeof(int));

  // load initial condition
  if (strcmp(initialPositionFilename, "__none") != 0)
  {
    int m1, n1;
    ReadMatrixFromDisk_(initialPositionFilename, &m1, &n1, &uInitial);
    if ((m1 != 3*n) || (n1 != 1))
    {
      printf("Error: initial position matrix size mismatch.\n");
      exit(1);
    }
  }
  else if ((deformableObject == MASSSPRING) && (massSpringSystemSource == CHAIN))
  {
    uInitial = (double*) calloc (3*n, sizeof(double));
    int numParticles = massSpringSystem->GetNumParticles();
    for(int i=0; i<numParticles; i++)
    {
      uInitial[3*i+0] = 1.0 - ((numParticles == 1) ? 1.0 : 1.0 * i / (numParticles - 1));
      uInitial[3*i+1] = 1.0 - ((numParticles == 1) ? 0.0 : 1.0 * i / (numParticles-1));
      uInitial[3*i+2] = 0.0;
    }
  }
  else
    uInitial = (double*) calloc (3*n, sizeof(double));

  // load initial velocity
  if (strcmp(initialVelocityFilename, "__none") != 0)
  {
    int m1, n1;
    ReadMatrixFromDisk_(initialVelocityFilename, &m1, &n1, &velInitial);
    if ((m1 != 3*n) || (n1 != 1))
    {
      printf("Error: initial position matrix size mismatch.\n");
      exit(1);
    }
  }

  // load force loads
  if (strcmp(forceLoadsFilename, "__none") != 0)
  {
    int m1;
    ReadMatrixFromDisk_(forceLoadsFilename, &m1, &numForceLoads, &forceLoads);
    if (m1 != 3*n)
    {
      printf("Mismatch in the dimension of the force load matrix.\n");
      exit(1);
    }
  }

  // create force model, to be used by the integrator
  printf("Creating force model...\n");
  if (deformableObject == STVK)
  {
    printf("Force model: STVK\n");
    stVKForceModel = new StVKForceModel(stVKInternalForces, stVKStiffnessMatrix);
    forceModel = stVKForceModel;
  }

  if (deformableObject == COROTLINFEM)
  {
    printf("Force model: COROTLINFEM\n");
    TetMesh * tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
    if (tetMesh == NULL)
    { printf("Error: the input mesh is not a tet mesh (CLFEM deformable model).\n");exit(1);}
    CorotationalLinearFEM * corotationalLinearFEM;
    if (numInternalForceThreads == 0)
      corotationalLinearFEM = new CorotationalLinearFEM(tetMesh);
    else
      {printf("Invalid numInternalForceThreads %d \n",numInternalForceThreads);exit(0);}

      //corotationalLinearFEM = new CorotationalLinearFEMMT(tetMesh, numInternalForceThreads);

    corotationalLinearFEMForceModel = new CorotationalLinearFEMForceModel(corotationalLinearFEM, corotationalLinearFEM_warp);
    forceModel = corotationalLinearFEMForceModel;
  }

  if (deformableObject == LINFEM)
  {
    printf("Force model: LINFEM\n");
    LinearFEMForceModel * linearFEMForceModel = new LinearFEMForceModel(stVKInternalForces);
    forceModel = linearFEMForceModel;
  }

  if (deformableObject == INVERTIBLEFEM)
  {
    printf("Force model: INVERTIBLEFEM\n");
    TetMesh * tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);

    if (tetMesh == NULL)
    {
      printf("Error: the input mesh is not a tet mesh (Invertible FEM deformable model).\n");
      exit(1);
    }
    IsotropicMaterial * isotropicMaterial = NULL;
    // create the invertible material model
    if (strcmp(invertibleMaterialString, "StVK") == 0)
      invertibleMaterial = INV_STVK;
    if (strcmp(invertibleMaterialString, "neoHookean") == 0)
      invertibleMaterial = INV_NEOHOOKEAN;
    if (strcmp(invertibleMaterialString, "MooneyRivlin") == 0)
      invertibleMaterial = INV_MOONEYRIVLIN;

    switch (invertibleMaterial)
    {
      case INV_STVK:
      {
        isotropicMaterial = new StVKIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
        printf("Invertible material: StVK.\n");
        break;
      }

      case INV_NEOHOOKEAN:
        isotropicMaterial = new NeoHookeanIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
        printf("Invertible material: neo-Hookean.\n");
        break;

      case INV_MOONEYRIVLIN:
        isotropicMaterial = new MooneyRivlinIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
        printf("Invertible material: Mooney-Rivlin.\n");
        break;

      default:
        printf("Error: invalid invertible material type.\n");
        exit(1);
        break;
    }

    // create the invertible FEM deformable model
    IsotropicHyperelasticFEM * isotropicHyperelasticFEM;
    if (numInternalForceThreads == 0)
      isotropicHyperelasticFEM = new IsotropicHyperelasticFEM(tetMesh, isotropicMaterial, inversionThreshold, addGravity, g);
    else
      {printf("Invalid numInternalForceThreads %d \n",numInternalForceThreads);exit(0);}
      //isotropicHyperelasticFEM = new IsotropicHyperelasticFEMMT(tetMesh, isotropicMaterial, inversionThreshold, addGravity, g, numInternalForceThreads);

    // create force model for the invertible FEM class
    IsotropicHyperelasticFEMForceModel * isotropicHyperelasticFEMForceModel = new IsotropicHyperelasticFEMForceModel(isotropicHyperelasticFEM);
    forceModel = isotropicHyperelasticFEMForceModel;
  }

  // initialize the integrator
  printf("Initializing the integrator, n = %d...\n", n);
  printf("Solver type: %s\n", solverMethod);

  integratorBaseSparse = NULL;
  if(solver == IMPLICITNEWMARK)
  {
    implicitNewmarkSparse = new ImplicitNewmarkSparse(3*n, timeStep, massMatrix, forceModel, positiveDefinite, numFixedDOFs, fixedDOFs,
       dampingMassCoef, dampingStiffnessCoef, maxIterations, epsilon, newmarkBeta, newmarkGamma, numSolverThreads);
    integratorBaseSparse = implicitNewmarkSparse;
  }
  else if(solver == IMPLICITBACKWARDEULER)
  {
    implicitNewmarkSparse = new ImplicitBackwardEulerSparse(3*n, timeStep, massMatrix, forceModel, positiveDefinite, numFixedDOFs, fixedDOFs,
       dampingMassCoef, dampingStiffnessCoef, maxIterations, epsilon, numSolverThreads);
    integratorBaseSparse = implicitNewmarkSparse;
  }
  else if(solver == EULER)
  {
    int symplectic = 0;
    integratorBaseSparse = new EulerSparse(3*n, timeStep, massMatrix, forceModel, symplectic, numFixedDOFs, fixedDOFs, dampingMassCoef);
  }
  else if(solver == SYMPLECTICEULER)
  {
    int symplectic = 1;
    integratorBaseSparse = new EulerSparse(3*n, timeStep, massMatrix, forceModel, symplectic, numFixedDOFs, fixedDOFs, dampingMassCoef);
  }
  /*else if (solver == CENTRALDIFFERENCES)
  {
    integratorBaseSparse = new CentralDifferencesSparse(3*n, timeStep, massMatrix, forceModel, numFixedDOFs, fixedDOFs, dampingMassCoef, dampingStiffnessCoef, centralDifferencesTangentialDampingUpdateMode, numSolverThreads);
  }*/
  integratorBase = integratorBaseSparse;

  if(integratorBase == NULL)
  {  printf("Error: failed to initialize numerical integrator.\n"); exit(1);}

  // set integration parameters
  integratorBaseSparse->SetDampingMatrix(LaplacianDampingMatrix);
  integratorBase->ResetToRest();
  integratorBase->SetState(uInitial, velInitial);
  integratorBase->SetTimestep(timeStep/substepsPerTimeStep);

  if(implicitNewmarkSparse!=NULL)
  {
    implicitNewmarkSparse->UseStaticSolver(staticSolver);
    if (velInitial != NULL)
      implicitNewmarkSparse->SetState(implicitNewmarkSparse->Getq(), velInitial);
  }
  // clear fps buffer
  for(int i=0; i<fpsBufferSize; i++)
    fpsBuffer[i] = 0.0;
  for(int i=0; i<forceAssemblyBufferSize; i++)
    forceAssemblyBuffer[i] = 0.0;
  for(int i=0; i<systemSolveBufferSize; i++)
    systemSolveBuffer[i] = 0.0;
  printf("heavingAmplitude = %lf   heaving Frequency = %lf \n",heavingAmplitude,heavingFrequency);
  titleBarCounter.StartCounter();
  printf(" Vega initiation finished \n");
}

void restrictOneDirectionMoving(char direction_)
{
    int i,j;
    integratorBase->GetqState(u,uvel,uaccel);
    switch(direction_)
    {
        case 'x': j=0; break;
        case 'y': j=1; break;
        case 'z': j=2; break;
        default : printf("ERROR: Invalid OneDirectionMoving %c exit.. \n",direction_); exit(0);
    }
    for(i=0;i<n;i++)
    {   u[i*3+j]=0.0; uvel[i*3+j]=0.0; uaccel[i*3+j]=0.0;}
    integratorBase->SetqState(u,uvel,uaccel);
}

void readVerticesRelation_vega_IB(char *verticesRelationFileName_)
{
    FILE *fp; int i,j;
    fp=fopen(verticesRelationFileName_,"r"); if(fp==NULL){printf("Error:cannot open '%s'",verticesRelationFileName_);exit(EXIT_FAILURE);}
    fscanf(fp,"%d \n",&i);
    while(i>0)
    {
       fscanf(fp,"%d ",&j);
       fscanf(fp,"%d %d %d \n",&verticesRelation_vega_IB[3*j-3],&verticesRelation_vega_IB[3*j-2],&verticesRelation_vega_IB[3*j-1]);
       verticesRelation_vega_IB[3*j-3]--; verticesRelation_vega_IB[3*j-2]--; verticesRelation_vega_IB[3*j-1]--;
       i--;
    }
    printf("read VerticesRelation_vega_IB Finished ! \n");
}

void readMovingVertices(char *movingVerticesFilename_,int *movingVertices_,int *numMovingVertices_)
{
    FILE *fp; int i=0,j;
    fp=fopen(movingVerticesFilename_,"r"); if(fp==NULL){printf("Error:cannot open '%s'",movingVerticesFilename_);exit(EXIT_FAILURE);}
    while(!feof(fp)){ fscanf(fp,"%d,",&movingVertices_[i]); i++;if(i>=1000){printf("");exit(0);}} i--;
    *numMovingVertices_=i; printf("Moving Vertices # %d: \n",i);
    for(j=0;j<i;j++)printf("%d ",movingVertices_[j]); printf("\n");
}

void readLoadingVertices(char *loadingVerticesFileName_,int *loadingVertices_,int *numLoadingVertices_)
{
    FILE *fp; int i=0,j;
    fp=fopen(loadingVerticesFileName_,"r"); if(fp==NULL){printf("Error:cannot open '%s'",loadingVerticesFileName_);exit(EXIT_FAILURE);}
    while(!feof(fp)){ fscanf(fp,"%d,",&loadingVertices_[i]); i++;if(i>=1000){printf("");exit(0);}} i--;
    *numLoadingVertices_=i; printf("External Force Loading Vertices # %d: \n",i);
    for(j=0;j<i;j++)printf("%d ",loadingVertices_[j]); printf("\n");
}

void setVerticesMoving(double timeStep_, int numTimeSteps_, int *movingVertices_, int numMovingVertices_)
{
    int i,temp_; double amplitude_=0.15,vibrateRatio_=0.4; FILE *fp; char appendLogoFileName[]="logo.file";
    //int movingVertices[]={49,50,51,52,101,102,103,104,153,154,155,156,205,206,207,208}, numMovingVertices=16;
    //movingVertices[]={49,51,101,103,153,155,205,207}, numMovingVertices=8;
    integratorBase->GetqState(u,uvel,uaccel);
    //fp=fopen(appendLogoFileName,"at");
    //fprintf(fp,"%d %lf %lf %lf \n",numTimeSteps_,timeStep_*numTimeSteps_,u[2],amplitude_*sin(timeStep_*numTimeSteps_*vibrateRatio_*3.1415926));
    //fclose(fp);

    for(i=0;i<numMovingVertices_;i++)
    {
        temp_=movingVertices_[i]-1;
        u[temp_*3]=0.0; u[temp_*3+1]=0.0;
        u[temp_*3+2]=amplitude_*sin(timeStep_*numTimeSteps_*vibrateRatio_*3.1415926);
        uvel[temp_*3]=0.0; uvel[temp_*3+1]=0.0; //uvel[temp_*3+2]=0.0;
        uvel[temp_*3+2]=amplitude_*cos(timeStep_*numTimeSteps_*vibrateRatio_*3.1415926)*vibrateRatio_*3.1415926;
        uaccel[temp_*3]=0.0; uaccel[temp_*3+1]=0.0; uaccel[temp_*3+2]=0.0;
        uaccel[temp_*3+2]=-1.0*amplitude_*sin(timeStep_*numTimeSteps_*vibrateRatio_*3.1415926)*vibrateRatio_*3.1415926*vibrateRatio_*3.1415926;
    }
    integratorBase->SetqState(u,uvel,uaccel);
}

void setExternalForces(double timeStep_,int numTimeSteps_,int *externalForcesVertices_,int numExternalForcesVertices_)
{
    int i,temp_; double correctionFactor_=1.23;//memset(f_ext,0,3*n*sizeof(double));
    for(i=0;i<numExternalForcesVertices_;i++)
    {
        temp_=externalForcesVertices_[i]-1; f_ext[temp_*3]=0.0; f_ext[temp_*3+1]=0.0;
        //f_ext[temp_*3+2]=7.0*sin(timeStep_*numTimeSteps_*1.0*3.1415926);
        f_ext[temp_*3+2]=1.0/numExternalForcesVertices_*correctionFactor_;
    }
    integratorBase->SetExternalForces(f_ext);
}

void filterForce(int bodyType_,double *force_)
{
    int i,j,temp,before_,after_; double ratio_=6.0;
    if(bodyType_==2)temp=n/2; else temp=n;

    for(j=0;j<3;j++)
    {
        if(fabs(force_[j])/fabs(force_[j+3])>=ratio_)force_[j]=force_[j+3];
    }
    for(i=1;i<temp;i++)
    {
        before_=3*i-3; if(i==temp-1)after_=3*i-3; else after_=3*i+3;
        for(j=0;j<3;j++)
        {
            if((fabs(force_[i*3+j])/fabs(force_[before_+j])>=ratio_)&&(fabs(force_[i*3+j])/fabs(force_[after_+j])>=ratio_/2.0))
                force_[i*3+j]=0.5*(force_[before_+j]+force_[after_+j]);
        }
    }
}

void readFlowGridXYZ(void)
{
    FILE *fp; char *flowGridName[3]={"xgrid.dat","ygrid.dat","zgrid.dat"}; int i,j,temp_; double previous_,current_;
    for(i=0;i<3;i++)
    {
        flowGrid[i]=(double*)calloc(nFlowGrid[i],sizeof(double));
        fp=fopen(flowGridName[i],"r");
        fscanf(fp,"%d %lf \n",&temp_,&previous_);
        for(j=0;j<nFlowGrid[i];j++)
        {
            fscanf(fp,"%d %lf \n",&temp_,&current_);
            flowGrid[i][j]=0.5*(previous_+current_);
            previous_=current_;
        }
        fclose(fp);
        printf(" Vega read flow Grid %d #%d  %s finished \n",i,nFlowGrid[i],flowGridName[i]);
    }
}

void readSurfaceTriangleMesh(void)
{
    FILE *fp; int i;
    fp=fopen("vegaSurfaceMesh.dat","r");
    fscanf(fp,"%d \n",&nTriangleMesh);
    surfaceTriangleMesh=(int*)calloc(3*nTriangleMesh,sizeof(int));
    for(i=0;i<nTriangleMesh;i++)
        fscanf(fp,"%d %d %d \n",&surfaceTriangleMesh[i*3],&surfaceTriangleMesh[i*3+1],&surfaceTriangleMesh[i*3+2]);
    fclose(fp); printf("read SurfaceTriangleMesh finished \n");
}

void obtainBodyNormalVectorCW(void)
{
    int nVertices_,i,j,k,before_,after_; int *verticesMark_; double vec1[3],vec2[3],CrossProductVec12[3],vecNorm_;
    if(bodyType==2)nVertices_=n/2; else if(bodyType==3) nVertices_=n; else {printf("ERROR invalid bodyType \n");exit(0);}
    verticesMark_=(int*)calloc(nVertices_,sizeof(int)); integratorBase->GetqState(u);
    for(i=0;i<nTriangleMesh;i++)
    {
        for(j=0;j<3;j++)
        {
            if(verticesMark_[surfaceTriangleMesh[i*3+j]]!=1)
            {
                if(j==0)before_=2; else before_=j-1;
                if(j==2)after_=0; else after_=j+1;
                Vec3d v1 = *(volumetricMesh->getVertex(surfaceTriangleMesh[i*3+before_]));
                Vec3d v2 = *(volumetricMesh->getVertex(surfaceTriangleMesh[i*3+j]));
                Vec3d v3 = *(volumetricMesh->getVertex(surfaceTriangleMesh[i*3+after_]));
                for(k=0;k<3;k++)
                {
                    v1[k]=v1[k]+u[surfaceTriangleMesh[i*3+before_]*3+k];
                    v2[k]=v2[k]+u[surfaceTriangleMesh[i*3+j]*3+k];
                    v3[k]=v3[k]+u[surfaceTriangleMesh[i*3+after_]*3+k];
                }
                for(k=0;k<3;k++)
                {
                    vec1[k]=v3[k]-v2[k];
                    vec2[k]=v1[k]-v2[k];
                }
                // CrossProductVec12 = vec1 X vec2
                CrossProductVec12[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
                CrossProductVec12[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
                CrossProductVec12[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
                // normalize CrossProductVec12
                vecNorm_=sqrt(pow(CrossProductVec12[0],2.0)+pow(CrossProductVec12[1],2.0)+pow(CrossProductVec12[2],2.0));
                for(k=0;k<3;k++)CrossProductVec12[k]=CrossProductVec12[k]/vecNorm_;
                for(k=0;k<3;k++)bodyNormalVectorCW[surfaceTriangleMesh[i*3+j]*3+k]=CrossProductVec12[k];
                verticesMark_[surfaceTriangleMesh[i*3+j]]=1;
            }
        }
    }
    free(verticesMark_);
    printf("obtain Body Normal Vector CW Finished \n");
}

int binarySearchIndex(int start_,int end_,double *array_,double data_)
{
    int middle_=(start_+end_)/2;
    if(start_==end_||start_==(end_-1))return start_;
    if(data_<array_[middle_])return binarySearchIndex(start_,middle_,array_,data_);
    else if(data_==array_[middle_])return middle_;
    else return binarySearchIndex(middle_,end_,array_,data_);
}

void obtainMarkerInterpolateIndexRatio(int *markerInterpolateIndex_,double *markerInterpolateRatio_)
{
    double vectorNorm_,ratio_; int nVertices_,i,j,index_;
    if(bodyType==2)nVertices_=n/2; else if(bodyType==3) nVertices_=n; else;
    printf("bodyType= %d nVertices = %d \n",bodyType,nVertices_);
    integratorBase->GetqState(u); vectorNorm_=1.8*fabs(flowGrid[0][nFlowGrid[0]/2]-flowGrid[0][nFlowGrid[0]/2-1]);
    for(i=0;i<nVertices_;i++)
    {
        Vec3d v=*(volumetricMesh->getVertex(i));
        for(j=0;j<3;j++)
        {
            v[j]=v[j]+u[i*3+j]+bodyNormalVectorCW[i*3+j]*vectorNorm_;
            index_=binarySearchIndex(nFlowGrid[j]/5,nFlowGrid[j]*4/5,flowGrid[j],v[j]);
            ratio_=(flowGrid[j][index_+1]-v[j])/(flowGrid[j][index_+1]-flowGrid[j][index_]);
            //printf("Index=%d  ratio=%lf ",index_,ratio_);
            markerInterpolateIndex_[i*6+j]=index_; markerInterpolateRatio_[i*6+j]=ratio_;
            v[j]=v[j]-2.0*bodyNormalVectorCW[i*3+j]*vectorNorm_;
            index_=binarySearchIndex(nFlowGrid[j]/5,nFlowGrid[j]*4/5,flowGrid[j],v[j]);
            ratio_=(flowGrid[j][index_+1]-v[j])/(flowGrid[j][index_+1]-flowGrid[j][index_]);
            markerInterpolateIndex_[i*6+3+j]=index_; markerInterpolateRatio_[i*6+3+j]=ratio_;
        }
    }
    printf("obtain Marker Interpolate Index & Ratio Finished \n");
}

void obtainMarkerForceFromPressureAndVelocity(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_)
{
    int nVertices_,i,j,k; double vectorNorm_,shearForce_[3],surfaceGradient[3][3];
    vectorNorm_=1.8*fabs(flowGrid[0][nFlowGrid[0]/2]-flowGrid[0][nFlowGrid[0]/2-1]);
    if(bodyType==2)nVertices_=n/2; else if(bodyType==3) nVertices_=n; else;
    for(i=0;i<nVertices_;i++)
    {
        // calculate pressure induced force
        for(j=0;j<3;j++)
            markerForce[i*3+j]=(markerPressure_[i*2+1]-markerPressure_[i*2])*markerArea*bodyNormalVectorCW[i*3+j];
        // calculate viscous shear flow induced force
        for(j=0;j<3;j++)
            shearForce_[j]=(markerInterpolateVelocity_[i*6+j]+markerInterpolateVelocity_[i*6+3+j]-2.0*bodyMotion_[i*3+j])/vectorNorm_/flowRe*markerArea;
        // surface gradient operator = I - nn
        for(j=0;j<3;j++)
        for(k=0;k<3;k++)
        {
            surfaceGradient[j][k]=-1.0*bodyNormalVectorCW[i*3+j]*bodyNormalVectorCW[i*3+k];
            if(j==k)surfaceGradient[j][k]=surfaceGradient[j][k]+1.0;
        }
        // (I-nn).shearForce
        for(j=0;j<3;j++)
            markerForce[i*3+j]=markerForce[i*3+j]+surfaceGradient[j][0]*shearForce_[0]
                        +surfaceGradient[j][1]*shearForce_[1]+surfaceGradient[j][2]*shearForce_[2];
    }
}

void saveVegaMarkerForce(char *fileMarker,double *markerForce_)
{
    int numVertices_,numElements_; FILE *fp; char bodyForceFileName[30];
    numVertices_=volumetricMesh->getNumVertices();
    numElements_=volumetricMesh->getNumElements();
    integratorBase->GetqState(u);
    sprintf(bodyForceFileName,"%s_%d.dat",fileMarker,timestepCounter);
    fp=fopen(bodyForceFileName,"w"); printf("*** Vega write %s \n",bodyForceFileName);
    fprintf(fp, "TITLE = \" 3 Dimensional Tetrahedral FE mesh data \" \n");
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"fx\", \"fy\", \"fz\" \n");
    fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON \n",numVertices_,numElements_);
    if(bodyType==2) // membrane
    {
        for(int i=0; i < numVertices_; i++)
        {
            Vec3d v = *(volumetricMesh->getVertex(i));
            if(i<numVertices_/2)
                fprintf(fp,"%.15G %.15G %.15G %.15G %.15G %.15G \n",v[0]+u[i*3],
                    v[1]+u[i*3+1]+heavingAmplitude*sin(timestepCounter*timeStep/heavingFrequency*2.0*3.14159265359),
                    v[2]+u[i*3+2],markerForce_[i*3],markerForce_[i*3+1],markerForce_[i*3+2]);
            else fprintf(fp,"%.15G %.15G %.15G %.15G %.15G %.15G \n",v[0]+u[i*3],
                    v[1]+u[i*3+1]+heavingAmplitude*sin(timestepCounter*timeStep/heavingFrequency*2.0*3.14159265359),
                    v[2]+u[i*3+2],0.0,0.0,0.0);
        }
    }
    else if(bodyType==3) // 3D solid body
    {
        for(int i=0; i < numVertices_; i++)
        {
            Vec3d v = *(volumetricMesh->getVertex(i));
            fprintf(fp,"%.15G %.15G %.15G %.15G %.15G %.15G \n",v[0]+u[i*3],
                v[1]+u[i*3+1]+heavingAmplitude*sin(timestepCounter*timeStep/heavingFrequency*2.0*3.14159265359),
                v[2]+u[i*3+2],markerForce_[i*3],markerForce_[i*3+1],markerForce_[i*3+2]);
        }
    }
    else {printf("**bodyType ERROR %d \n",bodyType); exit(0);}

    for(int el=0; el < numElements_; el++)
    {
        for(int j=0; j <=3; j++)
            fprintf(fp, "%d ", volumetricMesh->getVertexIndex(el,j)+1);
        fprintf(fp,"\n");
    }
    fclose(fp);
}

// main function
//extern "C" int vega_main_cpp(void);
extern "C" void vega_FEM_initiate_cpp(void);
extern "C" void vega_interpolateIndexRatio_cpp(int *markerInterpolateIndex_,double *markerInterpolateRatio_);
extern "C" void vega_deformation_cpp(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_);
extern "C" void vega_reNewBodyPosition_cpp(void);

void vega_FEM_initiate_cpp(void)
{
  char fileNameWrite[20]; int nVertices_;
  char argv[]="vega.config"; char movingVerticesFileName[]="PlateUpDown.moving.bou";
  char loadingVerticesFileName[]="PlateUpDown.forceLoadingVertices.bou"; char * configFilenameC = argv;

  configFilename = string(configFilenameC); printf("Loading scene configuration from %s.\n", configFilename.c_str());
  initConfigurations(); // parse the config file
  initSimulation(); // init the simulation
  //readVerticesRelation_vega_IB("verticesRelation_vega_IB.dat");
  sprintf(fileNameWrite,"TET_original.obj"); volumetricMesh->saveToFE_obj(fileNameWrite);
  readFlowGridXYZ(); readSurfaceTriangleMesh(); if(bodyType==2)nVertices_=n/2; else nVertices_=n;
  bodyNormalVectorCW=(double*)calloc(nVertices_*3,sizeof(double)); markerForce=(double*)calloc(nVertices_*3,sizeof(double));
  printf("********** Vega FSI Initiation finished !  timeStep %lf, \n",timeStep);
}

void vega_interpolateIndexRatio_cpp(int *markerInterpolateIndex_,double *markerInterpolateRatio_)
{
    obtainBodyNormalVectorCW(); obtainMarkerInterpolateIndexRatio(markerInterpolateIndex_,markerInterpolateRatio_);
}

void vega_deformation_cpp(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_)
{
  int i,temp;
  obtainMarkerForceFromPressureAndVelocity(markerPressure_,markerInterpolateVelocity_,bodyMotion_);
  if(timestepCounter%numOutputSteps==0)saveVegaMarkerForce("Vega_Force",markerForce);

  integratorBase->SetqState(u_base,uvel_base,uaccel_base); integratorBase->SetInternalForces(f_int_base);

  if(bodyType==3)
  {
      for(i=0;i<3*n;i++)f_ext[i]=markerForce[i]*forceRatio;
  }
  else if(bodyType==2)
  {
      for(i=0;i<3*n/2;i++)f_ext[i]=markerForce[i]*forceRatio;
      for(i=3*n/2;i<3*n;i++)f_ext[i]=0.0;
  }
  else ;
  filterForce(bodyType,f_ext);
  if(timestepCounter%numOutputSteps==0)saveVegaMarkerForce("Vega_Force_filtered",f_ext);

  integratorBase->SetExternalForces(f_ext);
  integratorBase->DoTimestep();
  integratorBase->GetqState(u,uvel,uaccel);
  printf("**** d_heavingAmplitude =%lf \n",heavingAmplitude*(sin(timestepCounter*timeStep/heavingFrequency*2.0*3.14159265359)-
         sin((timestepCounter-1)*timeStep/heavingFrequency*2.0*3.14159265359)));

  if(bodyType==3) temp=n; else if(bodyType==2)temp=n/2; else;
  for(i=0;i<temp;i++)
  {
      bodyMotion_[i*3]=(u[i*3]-u_base[i*3])/timeStep;
      bodyMotion_[i*3+1]=(u[i*3+1]+heavingAmplitude*sin(timestepCounter*timeStep/heavingFrequency*2.0*3.14159265359)
                             -u_base[i*3+1]-heavingAmplitude*sin((timestepCounter-1)*timeStep/heavingFrequency*2.0*3.14159265359))/timeStep;
      bodyMotion_[i*3+2]=(u[i*3+2]-u_base[i*3+2])/timeStep;
  }
  //obtainBodyNormalVectorCW(); obtainMarkerInterpolateIndexRatio(MarkerInterpolateIndex_,MarkerInterpolateRatio_);
  if(timestepCounter%numOutputSteps==0)saveVegaMarkerForce("Vega_markerVelocity",bodyMotion_);
}

void vega_reNewBodyPosition_cpp(void)
{
  integratorBase->GetqState(u_base,uvel_base,uaccel_base);
  integratorBase->GetInternalForces(f_int_base);
  timestepCounter++;
}

int vega_main_cpp(void)
{
  int numTimeSteps,i,j;  char fileNameWrite[20]; std::string num_; std::stringstream out;
  int movingVertices[1000],numMovingVertices=0,loadingVertices[1000],numLoadingVertices=0; FILE *fp;
  char appendLogoFileName[]="logo.dat"; double temp1_,temp2_;

  //char argv[]="C:/Dropbox/FSI/FEM/VegaFEM-v2.1/examples/asianDragon/asianDragon_MooneyRivlin.config";
  char argv[]="vega.config";
  char movingVerticesFileName[]="PlateUpDown.moving.bou";
  char loadingVerticesFileName[]="PlateUpDown.forceLoadingVertices.bou";
  char * configFilenameC = argv;

  printf("Starting application.\n");
  configFilename = string(configFilenameC);
  printf("Loading scene configuration from %s.\n", configFilename.c_str());

  initConfigurations(); // parse the config file
  initSimulation(); // init the simulation
  //readMovingVertices(movingVerticesFileName,movingVertices,&numMovingVertices);
  //readLoadingVertices(loadingVerticesFileName,loadingVertices,&numLoadingVertices);

  sprintf(fileNameWrite,"TET_original.obj"); volumetricMesh->saveToFE_obj(fileNameWrite);

  //numTimeSteps = 2000; numOutputSteps=100;
  printf("timeStep %lf, totalNumSteps %d, numOutput %d. Start time steps... \n",timeStep,numTimeSteps,numOutputSteps);
  //setExternalForces(timeStep,i,loadingVertices,numLoadingVertices);
  // Notice : cannot set load force and moving displacement before time steps !!!
  for(i=0;i<numTimeSteps;i++)
  {
      //setExternalForces(timeStep,i,loadingVertices,numLoadingVertices);
      //setVerticesMoving(timeStep,i,movingVertices,numMovingVertices);
      printf("\r steps #%5d ",i);
      integratorBase->DoTimestep();
      //restrictOneDirectionMoving('x');
      if(i%numOutputSteps==0)
      {
        integratorBase->GetqState(u); integratorBase->GetInternalForces(f_int);
        //sprintf(fileNameWrite,"TET_mesh_%d.obj",i); volumetricMesh->saveToFE_obj(fileNameWrite);
        sprintf(fileNameWrite,"TET_deformation_%d.obj",i); volumetricMesh->saveDeformationToFE_obj(fileNameWrite,u,f_int);
        printf("  write %s \n",fileNameWrite);
        /*fp=fopen(appendLogoFileName,"at"); fprintf(fp,"%lf ",i*timeStep);
        for(j=0;j<numLoadingVertices/2;j++)
        {
            temp1_=u[3*loadingVertices[1]-1]; temp2_=sqrt(pow(temp1_,2.0)+pow(u[3*loadingVertices[1]-2],2.0)+pow(u[3*loadingVertices[1]-3],2.0));
            fprintf(fp,"%lf %lf ",temp1_,temp2_);
        } fprintf(fp,"\n");fclose(fp);
        */

      }
  }
  return 0;
}

