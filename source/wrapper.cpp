/**
 * Minimal 2D Jacobi example. Code which is commented out demos how to
 * add a PPMWriter for output.
 */
/*
#include <libgeodecomp/io/simpleinitializer.h>
#include <libgeodecomp/io/ppmwriter.h>
#include <libgeodecomp/io/simplecellplotter.h>
#include <libgeodecomp/io/tracingwriter.h>
#include <libgeodecomp/parallelization/serialsimulator.h>
#include <libgeodecomp/parallelization/hiparsimulator.h>
*/
#include "../libgeodecomp/src/libgeodecomp.h"
#include <libflatarray/short_vec.hpp>
#include <cmath>

using namespace LibGeoDecomp;

class Cell
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasUpdateLineX,
        //public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >, 
	public APITraits::HasStencil<Stencils::Moore<3, 1> >, //contains all spatial neighbours
        //public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<3>,	// periodische randbedingung - cube ist konstant
        //public APITraits::HasCubeTopology<3>,
        public APITraits::HasPredefinedMPIDataType<double>       
    {};

    inline explicit Cell(double v = 0) :
        temp(v)
    {}

    template<typename NEIGHBORHOOD>
    static void updateLineX(Cell *target, long *x, long endX, const NEIGHBORHOOD& hood, const int nanoStep )
    {
        LIBFLATARRAY_LOOP_PEELER(double, 8, long, x, endX, updateLineImplmentation, target, hood, nanoStep);
    }

    template<typename DOUBLE, typename NEIGHBORHOOD>
    static void updateLineImplmentation(
    		long *x,
		long endX,
		Cell *target,
		const NEIGHBORHOOD& hood,
		const int /* nanoStep */)
    {

//	vectorized SoA
//	int num_vars = 1;
//	for(int i = 0; i < num_vars; ++i){
            DOUBLE factor = 1.0 / 5.0;
//            DOUBLE factor = 1.0 / 9.0;
//            DOUBLE factor = 1.0 / 7.0;
//            DOUBLE factor = 1.0 / 27.0;

            for (; *x < endX; ++*x) 
            {	
//		*** 2D5PT ***
                DOUBLE buf = hood[FixedCoord<-1,  0,  0>()].temp;
                	buf += hood[FixedCoord< 0, -1,  0>()].temp;
                        buf += hood[FixedCoord< 0,  0,  0>()].temp;
                        buf += hood[FixedCoord< 0,  1,  0>()].temp;
                        buf += hood[FixedCoord< 1,  0,  0>()].temp;
			buf *= factor;
			/*
// 		*** 2D9PT ***
                DOUBLE buf = hood[FixedCoord<-1, -1,  0>()].temp;
                	buf += hood[FixedCoord<-1,  0,  0>()].temp;
                	buf += hood[FixedCoord<-1,  1,  0>()].temp;
			
			buf += hood[FixedCoord< 0, -1,  0>()].temp;
			buf += hood[FixedCoord< 0,  0,  0>()].temp;
			buf += hood[FixedCoord< 0,  1,  0>()].temp;
			
			buf += hood[FixedCoord< 1, -1,  0>()].temp;
			buf += hood[FixedCoord< 1,  0,  0>()].temp;
			buf += hood[FixedCoord< 1,  1,  0>()].temp;
			buf *= factor;			
// 		*** 3D7PT ***
		DOUBLE buf = hood[FixedCoord<-1,  0,  0>()].temp;		
			buf += hood[FixedCoord< 0, -1,  0>()].temp;
			buf += hood[FixedCoord< 0,  0,  0>()].temp;
			buf += hood[FixedCoord< 0,  1,  0>()].temp;	
			buf += hood[FixedCoord< 1,  0,  0>()].temp;
			buf += hood[FixedCoord< 0,  0, -1>()].temp;
			buf += hood[FixedCoord< 0,  0,  1>()].temp;
			buf *= factor;
// 		*** 3D27PT ***
		DOUBLE buf = hood[FixedCoord<-1, -1, -1>()].temp;
			buf += hood[FixedCoord<-1,  0, -1>()].temp;
			buf += hood[FixedCoord<-1,  1, -1>()].temp;
			
			buf += hood[FixedCoord< 0, -1, -1>()].temp;
			buf += hood[FixedCoord< 0,  0, -1>()].temp;
			buf += hood[FixedCoord< 0,  1, -1>()].temp;
			
			buf += hood[FixedCoord< 1, -1, -1>()].temp;
			buf += hood[FixedCoord< 1,  0, -1>()].temp;
			buf += hood[FixedCoord< 1,  1, -1>()].temp;

			// *** MIDDLE ***
			buf += hood[FixedCoord<-1, -1,  0>()].temp;
			buf += hood[FixedCoord<-1,  0,  0>()].temp;
			buf += hood[FixedCoord<-1,  1,  0>()].temp;
			
			buf += hood[FixedCoord< 0, -1,  0>()].temp;
			buf += hood[FixedCoord< 0,  0,  0>()].temp;
			buf += hood[FixedCoord< 0,  1,  0>()].temp;
			
			buf += hood[FixedCoord< 1, -1,  0>()].temp;
			buf += hood[FixedCoord< 1,  0,  0>()].temp;
			buf += hood[FixedCoord< 1,  1,  0>()].temp;

			// *** FRONT ***
			buf += hood[FixedCoord<-1, -1,  1>()].temp;
			buf += hood[FixedCoord<-1,  0,  1>()].temp;
			buf += hood[FixedCoord<-1,  1,  1>()].temp;
			
			buf += hood[FixedCoord< 0, -1,  1>()].temp;
			buf += hood[FixedCoord< 0,  0,  1>()].temp;
			buf += hood[FixedCoord< 0,  1,  1>()].temp;
			
			buf += hood[FixedCoord< 1, -1,  1>()].temp;
			buf += hood[FixedCoord< 1,  0,  1>()].temp;	
			buf += hood[FixedCoord< 1,  1,  1>()].temp;
			buf *= factor;
*/ 			
		&target[*x].temp << buf;
	
            }
            
    }

    double temp;
};



class CellInitializer : public SimpleInitializer<Cell>
{
public:
    CellInitializer(const unsigned dimX, const unsigned dimY, const unsigned dimZ, const unsigned num_timesteps, int numVars, double *sourceTotal_, double *spikes_ , int *spikeLoc) :
    SimpleInitializer<Cell>(Coord<3>(dimX, dimY, dimZ), num_timesteps),
    dimX(dimX),
    dimY(dimY),
    dimZ(dimZ),
    numberOfVars(numVars)
    {
        sourceTotal = sourceTotal_;
        spikes = spikes_;
        spikeLocation = spikeLoc;
    }

    virtual void grid(GridBase<Cell, 3> *ret)
    {
		
        CoordBox<3> rect = ret->boundingBox();        
        for (CoordBox<3>::Iterator i = rect.begin(); i != rect.end(); ++i)
        {
	    //seede den random nr generator 
	    // RANDOM_NUMBER max?
	    ret->set(*i, Cell(Random::gen_d()));
	
	    // DEBUG_GRID == 1
	    //ret->set(*i, Cell(1.0));
	}
	
	/// Todo: Set multiple first spikes
        Coord<3> c( (unsigned) spikeLocation[1],
	            (unsigned) spikeLocation[2],
		    (unsigned) spikeLocation[3]);
		
        if (rect.inBounds(c)) 
        {
	    std::cout << "WARNING---------------------------------------------------\n"
	              << "WARNING: We're at Location " <<  spikeLocation[1] << ", " << spikeLocation[2] << ", " << spikeLocation[3] << "\n"
	              << "WARNING: and we 're setting the initial spike " << std::setprecision (15) << spikes[0] << " into the grid\n"
		      << "WARNING---------------------------------------------------\n";
	
            ret->set(c, Cell(spikes[0]));
        }
		    
        // update sourceTotal on every node
	for( int currentVar =0; currentVar< numberOfVars; currentVar++ )
        {
	    sourceTotal[currentVar] = sourceTotal[currentVar] + spikes[currentVar];
        }
    }

private:
    unsigned dimX;
    unsigned dimY;
    unsigned dimZ;
    int numberOfVars;
    double *sourceTotal;
    double *spikes;
    int *spikeLocation;
};


class ToleranceChecker : public Clonable<ParallelWriter<Cell>, ToleranceChecker>
{
public:
    typedef ParallelWriter<Cell>::GridType GridType;
	
    ToleranceChecker(const unsigned outputPeriod = 1, int numberOfVars = 1, double errorTol = 2500.0,  double *sourceTotal = NULL) :
    Clonable<ParallelWriter<Cell>, ToleranceChecker>("", outputPeriod),
    localSum(0),
    globalSum(0),
    num_vars(numberOfVars),
    err_tol(errorTol)
    {
        src_total = sourceTotal;
    }
	
    void stepFinished(
		const GridType& grid,
		const RegionType& validRegion,
		const CoordType& globalDimensions,
		unsigned step,
		WriterEvent event,
		std::size_t rank,
		bool lastCall)
    {		
	
	//for( int j = 0 ; j < num_vars ; ++j )
	//{
	    for (RegionType::Iterator i = validRegion.begin(); i != validRegion.end(); ++i) { 
	        localSum += grid.get(*i).temp;
	        //gridSum[j] += grid.get(*i).temp[j]; 
	    }
	//} 
			
	// so we're done summing the local grid
        if(lastCall)
        {
            //tag fuer einzigartigkeit
            MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            if( rank == 0 )
            	std::cout << "globalSum(" << step << ") = " << std::setprecision (15) << globalSum << "\n";
            
            // we're still missing the feature how many of our grids (variables in this case) should be summed
            //(see MG_BUFINIT.F)
            
            int j = 0;
            // checking error tolerance
            //for( int j = 0 ; j < num_vars ; ++j )
            //{
            	if( ( std::abs(src_total[j] - globalSum) / src_total[j] ) > err_tol )
            		std::cout << "error_tol(" << step << ") has not been met\n";
            //}
            
            // reset local and global sum since we're done with this step
            localSum = 0;
            globalSum = 0;
       }
    }
        
private: 
    double localSum;
    double globalSum; 
    int num_vars;
    double err_tol;
    double *src_total;
};

class SpikeJab : public Steerer<Cell>
{
public:
    using Steerer<Cell>::CoordType;
    using Steerer<Cell>::GridType;
    using Steerer<Cell>::Topology;

    SpikeJab (const unsigned ioPeriod, int numVars = 1, int numSpikes = 1, double *sourceTotal_ = NULL, double *spikes_ = NULL, int *spikeLocation_ = NULL) :
        Steerer<Cell>(ioPeriod),
        numberOfVars(numVars),
        numberOfSpikes(numSpikes)
    {
	sourceTotal = sourceTotal_;
	spikes = spikes_;
	spikeLocation = spikeLocation_;
    }

    void nextStep(
        GridType *grid,
        const Region<Topology::DIM>& validRegion,
        const CoordType& globalDimensions,
        unsigned step,
        SteererEvent event,
        std::size_t rank,
        bool lastCall,
        SteererFeedback *feedback)
    {
        // gibt es einen cooleren weg dem steerer zu sagen dass er den allerletzten step nichtmehr machen braucht?
        if ( step == numberOfSpikes * getPeriod() ) 
	{
	    return;
	}
		
	currentSpike = ( (int)step / getPeriod() );
	
	const Coord<3> currentCoord( spikeLocation[(currentSpike * 4) + 1], 
				     spikeLocation[(currentSpike * 4) + 2],
				     spikeLocation[(currentSpike * 4) + 3]);
	
	// setting spikes even in ghostzones							
	if(validRegion.count(currentCoord))
	{
	    Cell cell = grid->get(currentCoord);
	    cell.temp = spikes[(currentSpike * numberOfVars) + 0];
			
	    std::cout << "WARNING---------------------------------------------------\n"
		      << "WARNING: We're at Location " <<  spikeLocation[(currentSpike * 4) + 1] << ", " << spikeLocation[(currentSpike * 4) + 2] << ", " << spikeLocation[(currentSpike * 4) + 3] << "\n"
		      << "WARNING: We're rank: " << rank << "\n"
		      << "WARNING: and we 're setting spike " << std::setprecision (15) << cell.temp << " into grid at timestep " << step << "\n"
		      << "WARNING---------------------------------------------------\n";
		
			grid->set(currentCoord, cell);	
			
	}
	// everybody refresh sourceTotal - whether or not it belongs to this local grid
	if( lastCall )
	{
	    for( int currentVar = 0; currentVar < numberOfVars; currentVar++ )
	    {
	    	sourceTotal[currentVar] = sourceTotal[currentVar] + spikes[(currentSpike * numberOfVars) + currentVar];
	    }
	}
		  
    }

private:
    int numberOfVars;
    int numberOfSpikes;
    double *sourceTotal;
    double *spikes;
    int *spikeLocation;
    int currentSpike;
};


//LINE DOMINANT VS ROW DOMINANT - WATCH OUT !
//double spikes[*num_spikes][*num_vars], double spike_loc[*num_spikes][4]
extern "C" void simulate_(int *nx, int *ny, int *nz, int *num_vars, int *num_spikes, int *num_tsteps, double *err_tol, double *source_total, double *spikes, int *spike_loc)
{
    Typemaps::initializeMaps();
    {
	//SerialSimulator<Cell> sim(new CellInitializer(dimX, dimY, num_timesteps));
	CellInitializer *init =  new CellInitializer(*nx, *ny, *nz, (*num_tsteps * *num_spikes), *num_vars, source_total, spikes, spike_loc);
	
	//CheckerBoarding		
	HiParSimulator::HiParSimulator<Cell, RecursiveBisectionPartition<3> > sim(
			init,
			//MPILayer().rank() ? 0 : new TracingBalancer(new NoOpBalancer()),
			0,
			*num_tsteps,
			1);
	
	ToleranceChecker *toleranceChecker = new ToleranceChecker(1, *num_vars, *err_tol, source_total);
	sim.addWriter(toleranceChecker);

	SpikeJab *spikeJab = new SpikeJab(*num_tsteps, *num_vars, *num_spikes, source_total, spikes, spike_loc);
	sim.addSteerer(spikeJab);
	
	sim.run();
    }    
}
