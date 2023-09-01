
#include "mglet_cdef.h"

#include <iostream>
#include <chrono>
#include <cmath>

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkRectilinearGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPolyData.h"
#include "vtkDataSet.h"

#include "vtkDecimatePro.h"
#include "vtkContourFilter.h"
#include "vtkPlane.h"
#include "vtkPlaneCutter.h"
#include "vtkGlyph3D.h"
#include "vtkArrowSource.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkExtractGrid.h"

#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPPolyDataWriter.h"

#include "vtkMPIController.h"
#include "vtkProcessGroup.h"
#include <vtkMPI.h>


// just picking a tag which is available
static const int ISO_VALUE_RMI_TAG = 300;
static const int ISO_OUTPUT_TAG = 301;


// defining a struct of all information
struct ParallelIsoArgs_tmp{
    vtkMPIController* subcontroller;    // a seperate controller for I/O operations
    int* procPerIO;
    int argc; char** argv; int* retVal;
#ifdef _DOUBLE_PRECISION_
    double* xi[3]; double* dxi[3]; double* ddxi[3];     // for x,y,z, dx,dy,dz, ddx,ddy,ddz (1D arrays)
    double* array[21];                                  // for u,v,w,p etc. (3D arrays)
#else
    float* xi[3]; float* dxi[3]; float* ddxi[3];        // for x,y,z, dx,dy,dz, ddx,ddy,ddz (1D arrays)
    float* array[21];                                   // for u,v,w,p etc. (3D arrays)
#endif
    void (*cp_mgdims)(int*,int*,int*,const int*);       // function pointers
    void (*cp_get_ip1)(int*,const int*);
    void (*cp_get_ip3)(int*,const int*);
    void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*);
    void (*cp_iterate_grids_lvl)(const int*,int*,int*);
    void (*cp_compute_q)(const int*);
    void (*cp_compute_lambda2)(const int*);
    void (*cp_compute_p_extended)(const int*);
    void (*cp_compute_max_finecell)(const int*,int*);
    void (*cp_compute_event)(const int*);
    int* myid; int* numprocs; int* numgrids; int* istep; int* ilevel;
    int* nscal; int* itscal; int* lvlmin; int* lvlmax;
    double* isosurface;
    double* compression;
};


int get_ngrids_lvl( ParallelIsoArgs_tmp* args, const int ilevel ){
    int lvlcounter = 1;
    int igrid = -1;
    while ( true ) {
        args->cp_iterate_grids_lvl( &ilevel, &lvlcounter, &igrid );
        if ( igrid > 0 ){ lvlcounter++; } else { break; }
    }
    return (lvlcounter-1);
}


int get_ngrids_tot( ParallelIsoArgs_tmp* args ){
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);
    int counter = 0;
    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){
        int lvlcounter = 1;
        int igrid = -1;
        while ( true ) {
            args->cp_iterate_grids_lvl( &ilvl, &lvlcounter, &igrid );
            if ( igrid > 0 ){ lvlcounter++; } else { break; }
        }
        counter = counter + (lvlcounter-1);
    }
    return counter;
}


void get_grid_reduction( ParallelIsoArgs_tmp* args, const int igrid,
int& rfro, int& rbac, int& rrgt, int& rlft, int& rbot, int& rtop, int defval=1 ){
    int nfro = -1; int nbac = -1; int nrgt = -1;
    int nlft = -1; int nbot = -1; int ntop = -1;
    const int nbl = defval;
    rfro = nbl; rbac = nbl; rrgt = nbl; rlft = nbl; rbot = nbl; rtop = nbl;
    args->cp_mgbasb( &nfro, &nbac, &nrgt, &nlft, &nbot, &ntop, &igrid );
    // avoid getting to close the parent boundary
    // if ( nfro == 8 ){ rfro = 4; }
    // if ( nbac == 8 ){ rbac = 4; }
    // if ( nrgt == 8 ){ rrgt = 4; }
    // if ( nlft == 8 ){ rlft = 4; }
    // if ( nbot == 8 ){ rbot = 4; }
    // if ( ntop == 8 ){ rtop = 4; }
    return;
}


// This will be called by all processes
void lambda2Routine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkDecimatePro* deci;
    vtkContourFilter* iso;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);

    // determine file name
    const char* scalarName = "lambda2";
    char nmbr[9]; char fileName[26];  // 4+8 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName,"VTK/lambda2_");
    strcat(fileName,nmbr);
    strcat(fileName,".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){
        // results for Lambda2 are stored in field "hilf" (args->array[20])
        args->cp_compute_lambda2( &ilvl );

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );
                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    float val = (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    float val = (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[20];  // 20 = hilf (now containing lambda2 values)

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                // setting the scalar field that defines the isosurface
                // the loops realize the conversion from Fortran to C array
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj);
                            float val = (float) scalarInput[ip3-1+add];
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

    // computing iso-surface on multi-block data set
    const float isoValue = *(args->isosurface);

    iso = vtkContourFilter::New();
    iso->SetInputData(pds);
    iso->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
    iso->ComputeScalarsOn();

    iso->SetValue(0,  -1.0*isoValue);
    iso->SetValue(1,  -2.0*isoValue);
    iso->SetValue(2,  -4.0*isoValue);
    iso->SetValue(3, -16.0*isoValue);

    iso->GenerateTrianglesOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();

    // composing the data set to one polydata object (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(iso->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    const int ntri = comp->GetOutput()->GetNumberOfCells();
    deci = vtkDecimatePro::New();
    if ( ntri > 0 ) {
        // decimating to reduce data size before gathering and writing
        deci->SetInputData(comp->GetOutput());
        deci->SetTargetReduction(*(args->compression));
        deci->SetErrorIsAbsolute(1);
        deci->SetAbsoluteError(0.1);
        deci->PreserveTopologyOn();
        deci->Update();
    }

    if ( !isIOrank ) {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        if ( ntri > 0 ) {
            controller->Send(deci->GetOutput(), destRank, ISO_OUTPUT_TAG);
        } else {
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
    } else {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );
        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        if ( ntri > 0 ) {
            outputCopy->ShallowCopy(deci->GetOutput());
        } else {
            outputCopy->ShallowCopy(comp->GetOutput());
        }
        app->AddInputData(outputCopy);
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    iso->Delete();
    comp->Delete();
    deci->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
void eventRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkDecimatePro* deci;
    vtkContourFilter* iso;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);

    // determine file name
    const char* scalarName = "event";
    char nmbr[9]; char fileName[25];  // 4+7 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName,"VTK/event_");
    strcat(fileName,nmbr);
    strcat(fileName,".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){
        // computes a field that allows to extract sweep and eject events
        // (see Voermans 2018: "Hydr. Response to Coherent Turb. Motion")
        args->cp_compute_event( &ilvl );

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );
                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    float val = (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    float val = (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[20];  // 20 = hilf (now containing event-indicating field values)

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                // setting the scalar field that defines the isosurface
                // the loops realize the conversion from Fortran to C array
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj);
                            float val = (float) scalarInput[ip3-1+add];
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

    // computing iso-surface on multi-block data set
    const float isoValue = *(args->isosurface);

    iso = vtkContourFilter::New();
    iso->SetInputData(pds);
    iso->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
    iso->ComputeScalarsOn();

    iso->SetValue(0, -2.0*isoValue);
    iso->SetValue(1, -1.0*isoValue);
    iso->SetValue(2, +1.0*isoValue);
    iso->SetValue(3, +2.0*isoValue);

    iso->GenerateTrianglesOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();

    // composing the data set to one polydata object (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(iso->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    const int ntri = comp->GetOutput()->GetNumberOfCells();
    deci = vtkDecimatePro::New();
    if ( ntri > 0 ) {
        // decimating to reduce data size before gathering and writing
        deci->SetInputData(comp->GetOutput());
        deci->SetTargetReduction(*(args->compression));
        deci->SetErrorIsAbsolute(1);
        deci->SetAbsoluteError(0.1);
        deci->PreserveTopologyOn();
        deci->Update();
    }

    if ( !isIOrank ) {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        if ( ntri > 0 ) {
            controller->Send(deci->GetOutput(), destRank, ISO_OUTPUT_TAG);
        } else {
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
    } else {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );
        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        if ( ntri > 0 ) {
            outputCopy->ShallowCopy(deci->GetOutput());
        } else {
            outputCopy->ShallowCopy(comp->GetOutput());
        }
        app->AddInputData(outputCopy);
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes.
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    iso->Delete();
    comp->Delete();
    deci->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}


// This will be called by all processes
void uRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkDecimatePro* deci;
    vtkContourFilter* iso;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);

    // determine file name
    const char* scalarName = "u";
    char nmbr[9]; char fileName[20];  // 4+2 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/u_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );
                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    // staggering [ 1 0 0 ] is considered
                    float val = (float) ( 0.5*px[i] + 0.5*px[i+1] );
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    float val = (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[0];  // 0 = velocity U
                auto* scalarAvgInput = args->array[10];  // 10 = velocity U_AVG

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                // setting the scalar field that defines the isosurface
                // the loops realize the conversion from Fortran to C array
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj);
                            float val = (float) ( scalarInput[ip3-1+add] - scalarAvgInput[ip3-1+add] );
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

    // computing iso-surface on multi-block data set
    const float isoValue = *(args->isosurface);
    iso = vtkContourFilter::New();
    iso->SetInputData(pds);
    iso->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
    iso->ComputeScalarsOn();

    iso->SetValue(0, -2.0*isoValue);
    iso->SetValue(1, -1.0*isoValue);
    iso->SetValue(2,  1.0*isoValue);
    iso->SetValue(3,  2.0*isoValue);

    iso->GenerateTrianglesOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();

    // composing the data set (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(iso->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    const int ntri = comp->GetOutput()->GetNumberOfCells();
    deci = vtkDecimatePro::New();
    if ( ntri > 0 ) {
        // decimating to reduce data size before gathering and writing
        deci->SetInputData(comp->GetOutput());
        deci->SetTargetReduction(*(args->compression));
        deci->SetErrorIsAbsolute(1);
        deci->SetAbsoluteError(0.1);
        deci->PreserveTopologyOn();
        deci->Update();
    }

    if ( !isIOrank ) {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        if ( ntri > 0 ) {
            controller->Send(deci->GetOutput(), destRank, ISO_OUTPUT_TAG);
        } else {
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
    } else {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );
        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        if ( ntri > 0 ) {
            outputCopy->ShallowCopy(deci->GetOutput());
        } else {
            outputCopy->ShallowCopy(comp->GetOutput());
        }
        app->AddInputData(outputCopy);
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes.
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    iso->Delete();
    comp->Delete();
    deci->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
void wRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkDecimatePro* deci;
    vtkContourFilter* iso;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);

    // determine file name
    const char* scalarName = "w";
    char nmbr[9]; char fileName[20];  // 4+2 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/w_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );
                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    float val = (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    // staggering [ 0 0 1 ] is considered
                    float val = (float) ( 0.5*pz[k] + 0.5*pz[k+1] );
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[2];  // 2 = velocity W
                auto* scalarAvgInput = args->array[12];  // 12 = velocity W_AVG

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                // setting the scalar field that defines the isosurface
                // the loops realize the conversion from Fortran to C array
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj);
                            float val = (float) ( scalarInput[ip3-1+add] - scalarAvgInput[ip3-1+add] );
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

    // computing iso-surface on multi-block data set
    const float isoValue = *(args->isosurface);
    iso = vtkContourFilter::New();
    iso->SetInputData(pds);
    iso->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
    iso->ComputeScalarsOn();

    iso->SetValue(0, -2.0*isoValue);
    iso->SetValue(1, -1.0*isoValue);
    iso->SetValue(2,  1.0*isoValue);
    iso->SetValue(3,  2.0*isoValue);

    iso->GenerateTrianglesOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();

    // composing the data set (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(iso->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    const int ntri = comp->GetOutput()->GetNumberOfCells();
    deci = vtkDecimatePro::New();
    if ( ntri > 0 ) {
        // decimating to reduce data size before gathering and writing
        deci->SetInputData(comp->GetOutput());
        deci->SetTargetReduction(*(args->compression));
        deci->SetErrorIsAbsolute(1);
        deci->SetAbsoluteError(0.1);
        deci->PreserveTopologyOn();
        deci->Update();
    }

    if ( !isIOrank ) {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        if ( ntri > 0 ) {
            controller->Send(deci->GetOutput(), destRank, ISO_OUTPUT_TAG);
        } else {
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
    } else {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );
        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        if ( ntri > 0 ) {
            outputCopy->ShallowCopy(deci->GetOutput());
        } else {
            outputCopy->ShallowCopy(comp->GetOutput());
        }
        app->AddInputData(outputCopy);
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes.
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    iso->Delete();
    comp->Delete();
    deci->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
void tiRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkDecimatePro* deci;
    vtkContourFilter* iso;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);
    const int nscal = *(args->nscal);

    // ID of the scalar
    const int iTscal = *(args->itscal);
    assert( iTscal > 0 && "ERROR: Argument iTscal starts with 1 (Fortran / MGLET style)" );

    // determine file name
    const char* scalarName = "t";
    char nmbr[9]; char fileName[21];  // 4+3 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    if ( iTscal == 1 ){
        strcpy(fileName, "VTK/t1_");
    } else if ( iTscal == 2 ){
        strcpy(fileName, "VTK/t2_");
    } else if ( iTscal == 3 ){
        strcpy(fileName, "VTK/t3_");
    } else if ( iTscal == 4 ){
        strcpy(fileName, "VTK/t4_");
    } else {
        std::cout << "Invalid number for scalar has been detected." << std::endl;
    }
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );

                // computing the scalar pointer
                int ip3t = nscal * ip3 - (nscal-1);

                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    float val = (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    float val = (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[4];  // 4 = scalar T
                auto* scalarAvgInput = args->array[14];  // 14 = scalar T_AVG

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                // setting the scalar field that defines the isosurface
                // the loops realize the conversion from Fortran to C array
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj) + (kk*jj*ii) * (iTscal-1);
                            float val = (float) ( scalarInput[ip3t-1+add] - scalarAvgInput[ip3t-1+add] );
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

    // computing iso-surface on multi-block data set
    const float isoValue = *(args->isosurface);
    iso = vtkContourFilter::New();
    iso->SetInputData(pds);
    iso->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
    iso->ComputeScalarsOn();

    iso->SetValue(0, -2.0*isoValue);
    iso->SetValue(1, -1.0*isoValue);
    iso->SetValue(2,  1.0*isoValue);
    iso->SetValue(3,  2.0*isoValue);

    iso->GenerateTrianglesOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();

    // composing the data set (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(iso->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    const int ntri = comp->GetOutput()->GetNumberOfCells();
    deci = vtkDecimatePro::New();
    if ( ntri > 0 ) {
        // decimating to reduce data size before gathering and writing
        deci->SetInputData(comp->GetOutput());
        deci->SetTargetReduction(*(args->compression));
        deci->SetErrorIsAbsolute(1);
        deci->SetAbsoluteError(0.1);
        deci->PreserveTopologyOn();
        deci->Update();
    }

    if ( !isIOrank ) {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        if ( ntri > 0 ) {
            controller->Send(deci->GetOutput(), destRank, ISO_OUTPUT_TAG);
        } else {
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
    } else {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );
        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        if ( ntri > 0 ) {
            outputCopy->ShallowCopy(deci->GetOutput());
        } else {
            outputCopy->ShallowCopy(comp->GetOutput());
        }
        app->AddInputData(outputCopy);
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes.
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    iso->Delete();
    comp->Delete();
    deci->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
// This function is to plot a field value interpolated on a vertical x-y slice
void pHorSliceRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkPlaneCutter* cutter;
    vtkPlane* plane;
    vtkRectilinearGrid* grid;
    vtkMultiBlockDataSet* mbds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int nbl = 1;
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = get_ngrids_lvl( args, *(args->ilevel) );
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );

    const char* scalarName = "p";

    // determine file name
    char nmbr[9]; char fileName[21];  // 4+3 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/pH_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // results for extended pressure field are stored in field "hilf" (args->array[20])
    // args->cp_compute_p_extended( args->ilevel );

    // Inserting grids into multi-block data set
    mbds = vtkMultiBlockDataSet::New();
    mbds->SetNumberOfBlocks( numGrids );

    for ( int iter = 1; iter <= numGrids; iter++ ){

        // getting the value of igrid (possibly also -1)
        args->cp_iterate_grids_lvl( args->ilevel, &iter, &igrid );
        assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");

        int ii; int jj; int kk; int ip1; int ip3;
        // getting the grid dimension and the pointers
        args->cp_mgdims( &kk, &jj, &ii, &igrid );
        args->cp_get_ip1( &ip1, &igrid );
        args->cp_get_ip3( &ip3, &igrid );

        vtkNew<vtkFloatArray> xCoords;
        auto* px = args->xi[0];
        int istart = ip1 - 1;
        int iend = ip1 + ii;
        for (int i = istart+nbl; i < iend-nbl; i++){
            float val = (float) px[i];
            xCoords->InsertNextValue( val );
        }

        vtkNew<vtkFloatArray> yCoords;
        auto* py = args->xi[1];
        int jstart = ip1 - 1;
        int jend = ip1 + jj;
        for (int j = jstart+nbl; j < jend-nbl; j++){
            float val = (float) py[j];
            yCoords->InsertNextValue( val );
        }

        vtkNew<vtkFloatArray> zCoords;
        auto* pz = args->xi[2];
        int kstart = ip1 - 1;
        int kend = ip1 + kk;
        for (int k = kstart+nbl; k < kend-nbl; k++){
            float val = (float) pz[k];
            zCoords->InsertNextValue( val );
        }

        auto* scalarInput = args->array[3];  // 3 = P
        auto* scalarAvgInput = args->array[13];  // 13 = scalar P_AVG

        vtkNew<vtkFloatArray> scalars;
        scalars->SetName(scalarName);

        // setting the scalar field that defines the isosurface
        // the loops realize the conversion from Fortran to C array
        for (int k = nbl; k < kk-nbl; k++){
            for (int j = nbl; j < jj-nbl; j++){
                for (int i = nbl; i < ii-nbl; i++){
                    int add = k + j*(kk) + i*(kk*jj);
                    float val = (float) ( scalarInput[ip3-1+add] - scalarAvgInput[ip3-1+add] );
                    scalars->InsertNextValue(val);
                }
            }
        }

        // configuring the rectilinear grids
        grid = vtkRectilinearGrid::New();
        grid->SetDimensions( ii-2*nbl, jj-2*nbl, kk-2*nbl );
        grid->SetXCoordinates(xCoords);
        grid->SetYCoordinates(yCoords);
        grid->SetZCoordinates(zCoords);
        grid->GetPointData()->SetScalars(scalars);

        mbds->SetBlock( iter, grid );
    }

    // this ranks is an I/O ranks and gathers the data to write a file
    vtkAppendPolyData* app = vtkAppendPolyData::New();

    for ( int iPlane = 0; iPlane < 3; iPlane++){

        // computing the parameters
        double origin_vec[3] = {0.0, 0.0, 0.0105 - ((float) iPlane) };
        double norm_vec[3] = {0.0, 0.0, 1.0};
        // create a plane instance with desired origin and normal
        plane = vtkPlane::New();
        plane->SetNormal(norm_vec);
        plane->SetOrigin(origin_vec);
        // create a cutter instance and pass the above plane
        cutter = vtkPlaneCutter::New();
        cutter->SetPlane(plane);
        cutter->SetInputData(mbds);
        cutter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
        // composing the data set (necessary filter for multi-block data set)
        comp = vtkCompositeDataGeometryFilter::New();
        comp->SetInputConnection(cutter->GetOutputPort());
        comp->UpdatePiece(myid, numProcs, 0);

        if ( !isIOrank ){
            // this rank is not responsible for I/O and sends its data to an I/O rank
            const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
        else {
            // determine the number of delivering ranks+
            const int numRecv = std::min( procsPerIOrank, numProcs-myid );

            // retrieve and collect all the requested data
            for (int i = 1; i < numRecv; i++){
                vtkPolyData* pd = vtkPolyData::New();
                controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
                app->AddInputData(pd);
                pd->Delete();
            }

            // append the local contribution of this rank
            vtkPolyData* outputCopy = vtkPolyData::New();
            outputCopy->ShallowCopy(comp->GetOutput());
            app->AddInputData(outputCopy);
            outputCopy->Delete();
        }

        cutter->Delete();
        plane->Delete();
        comp->Delete();
    }

    if ( isIOrank ){
        // final update of the appended data
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
    }

    // clean up objects in all processes
    for ( int iter = 1; iter <= numGrids; iter++ ){ mbds->GetBlock(iter)->Delete(); }
    app->Delete();
    mbds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
// This function is to plot a field value interpolated on a vertical x-z slice
// In contrast to other routines, cell values are set here
void scalarVerSliceRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkPlaneCutter* cutter;
    vtkPlane* plane;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);
    const int nscal = *(args->nscal);


    // determine file name
    const char* scalarName = "t";
    char nmbr[9]; char fileName[21];  // 4+3 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/tV_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // ID of the scalar
    const int iTscal = *(args->itscal);
    assert( iTscal > 0 && "ERROR: Argument iTscal starts with 1 (Fortran / MGLET style)" );

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );

                // computing the scalar pointer
                int ip3t = nscal * ip3 - (nscal-1);

                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop, 2);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac+1; i++){
                    float val = 0.5 * (float) px[i-1] + 0.5 * (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft+1; j++){
                    float val = 0.5 * (float) py[j-1] + 0.5 * (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop+1; k++){
                    float val = 0.5 * (float) pz[k-1] + 0.5 * (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac+1, jj-rrgt-rlft+1, kk-rbot-rtop+1 );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);

                auto* scalarInput = args->array[4];     // 3 = P
                auto* scalarAvgInput = args->array[14];  // 13 = scalar P_AVG

                vtkNew<vtkFloatArray> scalar1;
                scalar1->SetName("scalar_1");
                scalar1->SetNumberOfComponents(1);
                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj) + (kk*jj*ii) * (1-1);
                            float val = (float) ( scalarInput[ip3t-1+add] - 0.0 * scalarAvgInput[ip3t-1+add] );
                            scalar1->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }
                assert( sanityCounter == grid->GetNumberOfCells() && "ERROR: Grid disagrees with data content");
                grid->GetCellData()->AddArray(scalar1);

                vtkNew<vtkFloatArray> scalar2;
                scalar2->SetName("scalar_2");
                scalar2->SetNumberOfComponents(1);
                sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj) + (kk*jj*ii) * (2-1);
                            float val = (float) ( scalarInput[ip3t-1+add] - 0.0 * scalarAvgInput[ip3t-1+add] );
                            scalar2->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }
                assert( sanityCounter == grid->GetNumberOfCells() && "ERROR: Grid disagrees with data content");
                grid->GetCellData()->AddArray(scalar2);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

     // this ranks is an I/O ranks and gathers the data to write a file
    vtkAppendPolyData* app = vtkAppendPolyData::New();

    for ( int iPlane = 0; iPlane < 2; iPlane++){

        // setting the parameters
        double origin_vec[3] = {0.0, 0.01 + 10.0 * float(iPlane), 0.0};
        double norm_vec[3] = {0.0, 1.0, 0.0};
        // create a plane instance with desired origin and normal
        plane = vtkPlane::New();
        plane->SetNormal(norm_vec);
        plane->SetOrigin(origin_vec);
        // create a cutter instance and pass the above plane
        cutter = vtkPlaneCutter::New();
        cutter->SetPlane(plane);
        cutter->SetInputData(pds);
        // cutter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
        // composing the data set (necessary filter for multi-block data set)
        comp = vtkCompositeDataGeometryFilter::New();
        comp->SetInputConnection(cutter->GetOutputPort());
        comp->UpdatePiece(myid, numProcs, 0);

        if ( !isIOrank ){
            // this rank is not responsible for I/O and sends its data to an I/O rank
            const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
        else {
            // determine the number of delivering ranks+
            const int numRecv = std::min( procsPerIOrank, numProcs-myid );
            // retrieve and collect all the requested data
            for (int i = 1; i < numRecv; i++){
                vtkPolyData* pd = vtkPolyData::New();
                controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
                app->AddInputData(pd);
                pd->Delete();
            }
            // append the local contribution of this rank
            vtkPolyData* outputCopy = vtkPolyData::New();
            outputCopy->ShallowCopy(comp->GetOutput());
            app->AddInputData(outputCopy);
            outputCopy->Delete();
        }
        cutter->Delete();
        plane->Delete();
        comp->Delete();
    }

    if ( isIOrank ){
        // final update of the appended data
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
    }

    // clean up objects in all processes
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    app->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}


// This will be called by all processes
// This function is to plot a field value interpolated on a vertical x-z slice
void pVerSliceRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkPlaneCutter* cutter;
    vtkPlane* plane;
    vtkRectilinearGrid* grid;
    vtkPartitionedDataSet* pds;
    vtkCompositeDataGeometryFilter* comp;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    int igrid;
    int prtCount = 0;
    int maxFineCell = 0;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = *(args->numgrids);
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );
    const int lmin = *(args->lvlmin);
    const int lmax = *(args->lvlmax);

    // determine file name
    const char* scalarName = "p";
    char nmbr[9]; char fileName[21];  // 4+3 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/pV_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    // Inserting grids into partitioned data set
    pds = vtkPartitionedDataSet::New();
    pds->SetNumberOfPartitions( numGrids );

    for ( int ilvl = lmin; ilvl <= lmax; ilvl++ ){

        // results for extended pressure field are stored in field "hilf" (args->array[20])
        // args->cp_compute_p_extended( &ilvl );

        const int ngrid_lvl = get_ngrids_lvl( args, ilvl );
        for ( int iter = 1; iter <= ngrid_lvl; iter++ ){

            // getting the value of igrid (possibly also -1)
            args->cp_iterate_grids_lvl( &ilvl, &iter, &igrid );
            assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");
            args->cp_compute_max_finecell( &igrid, &maxFineCell );

            if ( maxFineCell == 1 || ilvl == lmax ) {
                int ii; int jj; int kk; int ip1; int ip3;
                int rfro; int rbac; int rrgt; int rlft; int rbot; int rtop;

                // getting the grid dimension and the pointers
                args->cp_mgdims( &kk, &jj, &ii, &igrid );
                args->cp_get_ip1( &ip1, &igrid );
                args->cp_get_ip3( &ip3, &igrid );
                // determine reduction from actual grid size
                get_grid_reduction( args, igrid, rfro, rbac, rrgt, rlft, rbot, rtop);

                vtkNew<vtkFloatArray> xCoords;
                auto* px = args->xi[0];
                int istart = ip1 - 1;
                int iend = ip1 + ii;
                for (int i = istart+rfro; i < iend-rbac; i++){
                    float val = (float) px[i];
                    xCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> yCoords;
                auto* py = args->xi[1];
                int jstart = ip1 - 1;
                int jend = ip1 + jj;
                for (int j = jstart+rrgt; j < jend-rlft; j++){
                    float val = (float) py[j];
                    yCoords->InsertNextValue( val );
                }

                vtkNew<vtkFloatArray> zCoords;
                auto* pz = args->xi[2];
                int kstart = ip1 - 1;
                int kend = ip1 + kk;
                for (int k = kstart+rbot; k < kend-rtop; k++){
                    float val = (float) pz[k];
                    zCoords->InsertNextValue( val );
                }

                auto* scalarInput = args->array[3];     // 3 = P
                auto* scalarAvgInput = args->array[13];  // 13 = scalar P_AVG

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName(scalarName);

                int sanityCounter = 0;
                for (int k = rbot; k < kk-rtop; k++){
                    for (int j = rrgt; j < jj-rlft; j++){
                        for (int i = rfro; i < ii-rbac; i++){
                            int add = k + j*(kk) + i*(kk*jj);
                            float val = (float) ( scalarInput[ip3-1+add] - scalarAvgInput[ip3-1+add] );
                            scalars->InsertNextValue(val);
                            sanityCounter++;
                        }
                    }
                }

                // configuring the rectilinear grids
                grid = vtkRectilinearGrid::New();
                grid->SetDimensions( ii-rfro-rbac, jj-rrgt-rlft, kk-rbot-rtop );
                grid->SetXCoordinates(xCoords);
                grid->SetYCoordinates(yCoords);
                grid->SetZCoordinates(zCoords);
                assert( sanityCounter == grid->GetNumberOfPoints() && "ERROR: Grid disagrees with data content");
                grid->GetPointData()->SetScalars(scalars);

                pds->SetPartition( prtCount, grid );
                prtCount++;
            }
        }
    }

     // this ranks is an I/O ranks and gathers the data to write a file
    vtkAppendPolyData* app = vtkAppendPolyData::New();

    for ( int iPlane = 0; iPlane < 3; iPlane++){

        // setting the parameters
        double origin_vec[3] = {0.0, 4.0 + ((float) iPlane) * 20.0, 0.0};
        double norm_vec[3] = {0.0, 1.0, 0.0};
        // create a plane instance with desired origin and normal
        plane = vtkPlane::New();
        plane->SetNormal(norm_vec);
        plane->SetOrigin(origin_vec);
        // create a cutter instance and pass the above plane
        cutter = vtkPlaneCutter::New();
        cutter->SetPlane(plane);
        cutter->SetInputData(pds);
        // cutter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);
        // composing the data set (necessary filter for multi-block data set)
        comp = vtkCompositeDataGeometryFilter::New();
        comp->SetInputConnection(cutter->GetOutputPort());
        comp->UpdatePiece(myid, numProcs, 0);

        if ( !isIOrank ){
            // this rank is not responsible for I/O and sends its data to an I/O rank
            const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
            controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
        }
        else {
            // determine the number of delivering ranks+
            const int numRecv = std::min( procsPerIOrank, numProcs-myid );
            // retrieve and collect all the requested data
            for (int i = 1; i < numRecv; i++){
                vtkPolyData* pd = vtkPolyData::New();
                controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
                app->AddInputData(pd);
                pd->Delete();
            }
            // append the local contribution of this rank
            vtkPolyData* outputCopy = vtkPolyData::New();
            outputCopy->ShallowCopy(comp->GetOutput());
            app->AddInputData(outputCopy);
            outputCopy->Delete();
        }
        cutter->Delete();
        plane->Delete();
        comp->Delete();
    }

    if ( isIOrank ){
        // final update of the appended data
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
    }

    // clean up objects in all processes
    for ( int iter = 0; iter < prtCount; iter++ ){ pds->GetPartition(iter)->Delete(); }
    app->Delete();
    pds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// This will be called by all processes
// This function is to plot a field value interpolated on a slice
void velocityGlyphRoutine( vtkMultiProcessController* controller, void* arg )
{
    // pointers to VTK objects
    vtkGlyph3D* glyph;
    vtkArrowSource* arrow;
    vtkRectilinearGrid* grid;
    vtkMultiBlockDataSet* mbds;
    vtkCompositeDataGeometryFilter* comp;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    // starting the timer
    begin = std::chrono::steady_clock::now();

    // casting back the pointer
    ParallelIsoArgs_tmp* args = reinterpret_cast<ParallelIsoArgs_tmp*>(arg);

    // obtain all constant parameters form the argument
    const int nbl = 2;

    const int xStep = 2;
    const int yStep = 2;
    const int zStep = 2;

    const int myid = *(args->myid);
    const int numProcs = *(args->numprocs);
    const int numGrids = get_ngrids_lvl( args, *(args->ilevel) );
    const int procsPerIOrank = *(args->procPerIO);
    const bool isIOrank = ( (myid % procsPerIOrank) == 0 );

    const char* vectorName = "vel";

    // determine file name
    char nmbr[9]; char fileName[22];  // 4+4 for quantity, 8 for number, 5 for ending +1
    std::snprintf(nmbr, sizeof(nmbr), "%08d", *(args->istep) );
    strcpy(fileName, "VTK/vel_");
    strcat(fileName, nmbr);
    strcat(fileName, ".pvtp");

    int igrid;

    // Inserting grids into multi-block data set
    mbds = vtkMultiBlockDataSet::New();
    mbds->SetNumberOfBlocks( numGrids );

    for ( int iter = 1; iter <= numGrids; iter++ ){

        // getting the value of igrid (possibly also -1)
        args->cp_iterate_grids_lvl( args->ilevel, &iter, &igrid );
        assert( igrid > 0 && "ERROR: Invalid value of igrid detected.");

        if ( igrid > 0 ) {
            int ii; int jj; int kk; int ip1; int ip3;
            // getting the grid dimension and the pointers
            args->cp_mgdims( &kk, &jj, &ii, &igrid );
            args->cp_get_ip1( &ip1, &igrid );
            args->cp_get_ip3( &ip3, &igrid );

            vtkNew<vtkFloatArray> xCoords;
            auto* px = args->xi[0];
            int istart = ip1 - 1;
            int iend = ip1 + ii;
            for (int i = istart+nbl; i < iend-nbl; i+=xStep){
                float val = (float) px[i];
                xCoords->InsertNextValue( val );
            }

            vtkNew<vtkFloatArray> yCoords;
            auto* py = args->xi[1];
            int jstart = ip1 - 1;
            int jend = ip1 + jj;
            for (int j = jstart+nbl; j < jend-nbl; j+=yStep){
                float val = (float) py[j];
                yCoords->InsertNextValue( val );
            }

            vtkNew<vtkFloatArray> zCoords;
            auto* pz = args->xi[2];
            int kstart = ip1 - 1;
            int kend = ip1 + kk;
            for (int k = kstart+nbl; k < kend-nbl; k+=zStep){
                float val = (float) pz[k];
                zCoords->InsertNextValue( val );
            }

            auto* xVelocity = args->array[0];  // u
            auto* yVelocity = args->array[1];  // v
            auto* zVelocity = args->array[2];  // w

            vtkNew<vtkFloatArray> velocity;
            velocity->SetName(vectorName);
            velocity->SetNumberOfComponents(3);

            int nt = (ii-2*nbl)/xStep * (jj-2*nbl)/yStep * (kk-2*nbl)/zStep;
            velocity->SetNumberOfTuples( nt );

            // setting the velocity vector field
            // the loops realize the conversion from Fortran to C array
            int counter = 0;
            for (int k = nbl; k < kk-nbl; k+=zStep){
                for (int j = nbl; j < jj-nbl; j+=yStep){
                    for (int i = nbl; i < ii-nbl; i+=xStep){
                        int add = k + j*(kk) + i*(kk*jj);
                        float u = (float) xVelocity[ip3-1+add];
                        float v = (float) yVelocity[ip3-1+add];
                        float w = (float) zVelocity[ip3-1+add];
                        float abs = sqrt( u*u + v*v + w*w );
                        if ( abs > 0.00000001 ){
                            velocity->SetTuple3( counter, u/abs, v/abs, w/abs );
                        } else {
                            velocity->SetTuple3( counter, 0.0, 0.0, 0.0 );
                        }
                        counter++;
                    }
                }
            }

            assert( counter == nt && "ERROR: Invalid value of igrid detected.");

            // configuring the rectilinear grids
            grid = vtkRectilinearGrid::New();
            grid->SetDimensions( (ii-2*nbl)/xStep, (jj-2*nbl)/yStep, (kk-2*nbl)/zStep );
            grid->SetXCoordinates(xCoords);
            grid->SetYCoordinates(yCoords);
            grid->SetZCoordinates(zCoords);
            grid->GetPointData()->SetVectors(velocity);

            mbds->SetBlock( iter, grid );
        }
    }

    // create an instance of the glyph 3D filter
    arrow = vtkArrowSource::New();
    arrow->SetTipResolution(5);
    arrow->SetTipLength(0.2);
    arrow->SetTipRadius(0.1);
    arrow->Update();

    glyph = vtkGlyph3D::New();
    glyph->SetInputData(mbds);
    glyph->SetSourceData(arrow->GetOutput());
    glyph->SetVectorModeToUseVector();
    glyph->SetScaleModeToScaleByVector();
    glyph->SetScaleFactor(0.15);

    // composing the data set (necessary filter for multi-block data set)
    comp = vtkCompositeDataGeometryFilter::New();
    comp->SetInputConnection(glyph->GetOutputPort());
    comp->UpdatePiece(myid, numProcs, 0);

    if ( !isIOrank )
    {
        // this rank is not responsible for I/O and sends its data to an I/O rank
        const int destRank = ( myid / procsPerIOrank ) * procsPerIOrank;
        controller->Send(comp->GetOutput(), destRank, ISO_OUTPUT_TAG);
    }
    else
    {
        // this ranks is an I/O ranks and gathers the data to write a file
        vtkAppendPolyData* app = vtkAppendPolyData::New();

        // determine the number of delivering ranks+
        const int numRecv = std::min( procsPerIOrank, numProcs-myid );

        // retrieve and collect all the requested data
        for (int i = 1; i < numRecv; i++){
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, (myid+i), ISO_OUTPUT_TAG);
            app->AddInputData(pd);
            pd->Delete();
        }

        // append the local contribution of this rank
        vtkPolyData* outputCopy = vtkPolyData::New();
        outputCopy->ShallowCopy(comp->GetOutput());
        app->AddInputData(outputCopy);

        // final update of the appended data
        app->Update();

        vtkXMLPPolyDataWriter* parallelWriter = vtkXMLPPolyDataWriter::New();
        parallelWriter->SetInputData(app->GetOutput());
        parallelWriter->SetController(args->subcontroller);
        parallelWriter->SetFileName(fileName);
        parallelWriter->SetNumberOfPieces(args->subcontroller->GetNumberOfProcesses());
        parallelWriter->SetStartPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetEndPiece(args->subcontroller->GetLocalProcessId());
        parallelWriter->SetDataModeToBinary();
        parallelWriter->Update();
        parallelWriter->Write();

        // clean up objects present on I/O ranks
        parallelWriter->Delete();
        outputCopy->Delete();
        app->Delete();
    }

    // clean up objects in all processes
    glyph->Delete();
    arrow->Delete();
    comp->Delete();
    for ( int iter = 1; iter <= numGrids; iter++ ){ mbds->GetBlock(iter)->Delete(); }
    mbds->Delete();

    // starting the timer
    end = std::chrono::steady_clock::now();
    if ( myid == 0 ){
        std::cout << " VTK file " << fileName << " written in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()
                  << " ms" << std::endl;
    }
    return;
}



// Main function that is called from MGLET
// (known restriction: Fields are assumed to be in double precision)

extern "C" void render_vtk_function_c( void (*cp_mgdims)(int*,int*,int*,const int*), void (*cp_get_ip1)(int*,const int*),
                                       void (*cp_get_ip3)(int*,const int*), void (*cp_iterate_grids_lvl)(const int*,int*,int*),
                                       void (*cp_compute_q)(const int*), void (*cp_compute_lambda2)(const int*),
                                       void (*cp_compute_p_extended)(const int*), void (*cp_compute_max_finecell)(const int*,int*),
                                       void (*cp_mgbasb)(int*,int*,int*,int*,int*,int*,const int*), void (*cp_compute_event)(const int*),
#ifdef _DOUBLE_PRECISION_
                                       double* u, double* v, double* w, double* p, double* t, double* au, double* av, double* aw, double* ap, double* at,
                                       double* auu, double* avv, double* aww, double* auv, double* auw, double* avw,
                                       double* x, double* y, double* z, double* dx, double* dy, double* dz, double* ddx, double* ddy, double* ddz, double* hilf,
#else
                                       float* u, float* v, float* w, float* p, float* t, float* au, float* av, float* aw, float* ap, float* at,
                                       float* auu, float* avv, float* aww, float* auv, float* auw, float* avw,
                                       float* x, float* y, float* z, float* dx, float* dy, float* dz, float* ddx, float* ddy, float* ddz, float* hilf,
#endif
                                       int* myid, int* numprocs, int* istep, int* idim3d, int* idim2d, int* idim1d,
                                       int* nscal, int* lvlmin, int* lvlmax ) {

    // determine if active in this time step
    // - 20 for x96
    // - 10 for x48
    const int tstep_interval = 20;
    if ( (*istep) % tstep_interval != 0 ){ return; }

    // set the shear velocity
    // - for 04_DNS2: sqrt( 0.00005903 * 5.0 )
    const float utau = sqrt( 0.00005903 * 5.0 );

    // setting the parameters (tune for performance)
    int procPerIO = 2*28;

    // Note that this will create a vtkMPIController if MPI
    // is configured, vtkThreadedController otherwise.
    vtkMPIController* controller = vtkMPIController::New();
    int a; char** b;
    controller->Initialize( &a, &b, 1 );
    assert( controller->GetNumberOfProcesses() == *numprocs );
    assert( controller->GetLocalProcessId() == *myid );

    // creating a seperate controller for I/O ranks only
    vtkProcessGroup* subgroup = vtkProcessGroup::New();
    subgroup->Initialize(controller);
    subgroup->RemoveAllProcessIds();
    for (int iproc = 0; iproc < *numprocs; iproc+=procPerIO){ subgroup->AddProcessId(iproc); }
    vtkMPIController* subcontroller = controller->CreateSubController(subgroup);

    if ( (*myid) % procPerIO == 0 ){
        // this process is I/O rank and size of the subcontroller is checked
        assert( subcontroller->GetNumberOfProcesses() == ((*numprocs)+procPerIO-1)/procPerIO );
    } else{
        // this process is not an I/O rank and subcontroller is not defined (nullptr)
        assert( !subcontroller );
    }

    // Setting the arguments to be handed to the main routines
    // all of them are packaged in the struct "args"
    // ----------------------------------------------
    ParallelIsoArgs_tmp args;
    double isosurface = 0.0;
    int itscal = -1;
    int ilevel = -1;
    const int mglet_lvlmin = *lvlmin;
    const int mglet_lvlmax = *lvlmax;
    int mylvlmin = mglet_lvlmin;
    int mylvlmax = mglet_lvlmax;
    double compression = 0.95;

    args.subcontroller = subcontroller;
    args.procPerIO = &procPerIO;

    int retVal = 1; args.retVal = &retVal;
    args.isosurface = &isosurface;
    args.itscal = &itscal;
    args.myid = myid; args.numprocs = numprocs; args.istep = istep;
    args.ilevel = &ilevel;
    args.nscal = nscal;
    args.lvlmin = &mylvlmin; args.lvlmax = &mylvlmax;
    args.compression = &compression;

    // // wrapping the function pointers
    args.cp_mgdims = cp_mgdims; args.cp_get_ip1 = cp_get_ip1;
    args.cp_get_ip3 = cp_get_ip3; args.cp_iterate_grids_lvl = cp_iterate_grids_lvl;
    args.cp_mgbasb = cp_mgbasb;
    args.cp_compute_q = cp_compute_q;
    args.cp_compute_lambda2 = cp_compute_lambda2;
    args.cp_compute_p_extended = cp_compute_p_extended;
    args.cp_compute_max_finecell = cp_compute_max_finecell;
    args.cp_compute_event = cp_compute_event;

    // wrapping the coordinates / cell measures
    args.xi[0] = x; args.xi[1] = y; args.xi[2] = z;
    args.dxi[0] = dx; args.dxi[1] = dy; args.dxi[2] = dz;
    args.ddxi[0] = ddx; args.ddxi[1] = ddy; args.ddxi[2] = ddz;

    // wrapping the 3D fields (ordering to be improved)
    args.array[0] = u; args.array[1] = v; args.array[2] = w;
    args.array[3] = p; args.array[4] = t;
    args.array[10] = au; args.array[11] = av; args.array[12] = aw;
    args.array[13] = ap; args.array[14] = at;
    args.array[20] = hilf;

    int counter = get_ngrids_tot( &args );
    args.numgrids = &counter;

    // ----------------------------------------------

    // isosurface = 0.005;  // empirical value
    // controller->SetSingleMethod( lambda2Routine, &args );
    // controller->SingleMethodExecute();

    // ----------------------------------------------

    // ilevel = 2;
    // controller->SetSingleMethod( pHorSliceRoutine, &args );
    // controller->SingleMethodExecute();

    // ----------------------------------------------

    // settings to be made inside the function
    // e.g. for 48 cells/D: only sample every 8th cell
    // ilevel = 3;
    // controller->SetSingleMethod( velocityGlyphRoutine, &args );
    // controller->SingleMethodExecute();

    // ----------------------------------------------

    mylvlmin = 1; mylvlmax = 1;

    // isosurface = 1.0 * utau;  // setting u_tau
    // controller->SetSingleMethod( uRoutine, &args );
    // controller->SingleMethodExecute();

    // controller->SetSingleMethod( wRoutine, &args );
    // controller->SingleMethodExecute();

    // isosurface = 2.0 * utau * utau;  // (u_tau)^2
    // controller->SetSingleMethod( eventRoutine, &args );
    // controller->SingleMethodExecute();

    // isosurface = 0.05; itscal = 1;
    // controller->SetSingleMethod( tiRoutine, &args );
    // controller->SingleMethodExecute();

    // isosurface = 0.3333; itscal = 2;
    // controller->SetSingleMethod( tiRoutine, &args );
    // controller->SingleMethodExecute();

    // settings to be made inside the function
    itscal = 1;
    controller->SetSingleMethod( scalarVerSliceRoutine, &args );
    controller->SingleMethodExecute();

    // ----------------------------------------------

    subgroup->Delete();
    if ( (*myid) % procPerIO == 0 ){ subcontroller->Delete(); }
    controller->Finalize(1);  // argument "1" tells function not call MPI_Finalize internally
    controller->Delete();

    return;
}
