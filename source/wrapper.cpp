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
#include <cmath>

using namespace LibGeoDecomp;

class Cell
{
public:
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasStencil<Stencils::VonNeumann<2, 1> >,
        //public APITraits::HasOpaqueMPIDataType<Cell>,
        public APITraits::HasTorusTopology<2>,	// periodische randbedingung - cube ist konstant
        public APITraits::HasPredefinedMPIDataType<double>       
    {};

    inline explicit Cell(double v = 0) :
        temp(v)
    {}

    template<typename COORD_MAP>
    void update(const COORD_MAP& neighborhood, const unsigned& nanoStep)
    {
		//vectorized SoA
//	int num_vars = 1;
//	for(int i = 0; i < num_vars; ++i){
    
/*	        temp = (neighborhood[Coord<2>( 0, -1)].temp +
        	        neighborhood[Coord<2>(-1,  0)].temp +
               		neighborhood[Coord<2>( 0,  0)].temp +
                	neighborhood[Coord<2>( 1,  0)].temp +
                	neighborhood[Coord<2>( 0,  1)].temp) * (1.0 / 5.0);

*/
// 		*** 2D5PT ***

		temp = (neighborhood[FixedCoord<-1,  0>()].temp +
			neighborhood[FixedCoord< 0, -1>()].temp +
			neighborhood[FixedCoord< 0,  0>()].temp +
			neighborhood[FixedCoord< 0,  1>()].temp +
			neighborhood[FixedCoord< 1,  0>()].temp) * (1.0 / 5.0);
			
/*		temp = (neighborhood[Coord<2>(-1,  0)].temp +
			neighborhood[Coord<2>( 0, -1)].temp +
			neighborhood[Coord<2>( 0,  0)].temp +
			neighborhood[Coord<2>( 0,  1)].temp +
			neighborhood[Coord<2>( 1,  0)].temp) * (1.0 / 5.0);
*/
/*		temp[i] = (neighborhood[Coord<2>(-1,  0)].temp[i] +
			   neighborhood[Coord<2>( 0, -1)].temp[i] +
			   neighborhood[Coord<2>( 0,  0)].temp[i] +
			   neighborhood[Coord<2>( 0,  1)].temp[i] +
			   neighborhood[Coord<2>( 1,  0)].temp[i]) * (1.0 / 5.0);

*/
// 		*** 2D9PT ***
/*		temp = (neighborhood[Coord<2>(-1, -1)].temp +
			neighborhood[Coord<2>(-1,  0)].temp +
			neighborhood[Coord<2>(-1,  1)].temp +
			
			neighborhood[Coord<2>( 0, -1)].temp +
			neighborhood[Coord<2>( 0,  0)].temp +
			neighborhood[Coord<2>( 0,  1)].temp +
			
			neighborhood[Coord<2>( 1, -1)].temp +
			neighborhood[Coord<2>( 1,  0)].temp +
			neighborhood[Coord<2>( 1,  1)].temp) * (1.0 / 9.0);
*/			
//	}
   }

    double temp;
};



class CellInitializer : public SimpleInitializer<Cell>
{
public:
    CellInitializer(const unsigned dimX, const unsigned dimY, const unsigned num_timesteps, int numVars, double *sourceTotal_, double *spikes_ , int *spikeLoc) :
    SimpleInitializer<Cell>(Coord<2>(dimX, dimY), num_timesteps),
    dimX(dimX),
    dimY(dimY),
    numberOfVars(numVars)
    {
		sourceTotal = sourceTotal_;
		spikes = spikes_;
		spikeLocation = spikeLoc;
	}

    virtual void grid(GridBase<Cell, 2> *ret)
    {
        CoordBox<2> rect = ret->boundingBox();
        
        std::cout << "---------------------CELLINIT------------------------------\n"
				  << " X : " << dimX << " Y : " << dimY << "\n"
				  << " BoundingBox : " << rect.dimensions << "\n";

        for (unsigned int y = 0; y < dimY; ++y) 
        {
            for (unsigned int x = 0; x < dimX; ++x) 
            {
                Coord<2> c(x, y);
                if (rect.inBounds(c)) 
                {
					///TODO: array
					// initial spike
					if( (x == (unsigned) spikeLocation[1]) &&
						(y == (unsigned) spikeLocation[2]) )
					{
						std::cout << "WARNING---------------------------------------------------\n"
								  << "WARNING: We're at Location " <<  spikeLocation[1] << ", " << spikeLocation[2] << ", " << spikeLocation[3] << "\n"
								  << "WARNING: and we 're setting the initial spike " << std::setprecision (15) << spikes[0] << " into the grid\n"
								  << "WARNING---------------------------------------------------\n";
			
						ret->set(c, Cell(spikes[0]));
						
						for( int currentVar =0; currentVar< numberOfVars; currentVar++ )
						{
							sourceTotal[currentVar] = sourceTotal[currentVar] + spikes[currentVar];
						}
						
						
					}
					else
					{						
						// RANDOM_NUMBER max?
						//ret->set(c, Cell(Random::gen_d()));
						
						// DEBUG_GRID == 1
						ret->set(c, Cell(0.0)); 
					}
				}
            }
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
		if( step == 1 )
		{
			CoordBox<2> box = validRegion.boundingBox();      
			std::cout << "rank " << rank << "\'s dim: " << box.dimensions.toString() << "\n";
		}

		// valid region  = 48 || boundingbox = 52 ??? aber 50 waere eigentlich ganz nett 
		// valid region zieht die inneren ghost zones ab! 
		// haengt zusaetzlich noch bin der Topology ab

		// FRAGE: SIND HIER DIE GHOSTZONES DABEI??
		CoordBox<2> box = validRegion.boundingBox(); 
		
		//for( int j = 0 ; j < num_vars ; ++j )
		//{
			for (CoordBox<2>::Iterator i = box.begin(); i != box.end(); ++i) { 
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

			// checking error tolerance
			for( int j = 0 ; j < num_vars ; ++j )
			{
				if( ( std::abs(src_total[j] - globalSum) / src_total[j] ) > err_tol )
					std::cout << "error_tol(" << step << ") has not been met\n";
			}
			//grids_to_sum kram aus mg_bufinit.f    
				//      ERROR_ITER = ABS ( SOURCE_TOTAL(IVAR) - GSUM ) / SOURCE_TOTAL(IVAR)
			//      IF ( ERROR_ITER > ERROR_TOL ) THEN
			
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
		
		const Coord<2> currentCoord( spikeLocation[(currentSpike * 4) + 1], 
									 spikeLocation[(currentSpike * 4) + 2]);
		
		// setting spikes even in ghostzones							
		if(validRegion.count(currentCoord))
		{
			Cell cell = grid->get(currentCoord);
			cell.temp = spikes[(currentSpike * numberOfVars) + 0];
//				for( int currentVar =0; i< numberOfVars; currentVar++ )
//				{
//					cell.temp[currentVar] = spikes[(currentSpike * numberOfVars) + currentVar];//<< "\n";
//				}
				
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
/*
    if (MPILayer().rank() == 0) {
		printf("Dimensions: %i * %i * %i\n", *nx, *ny, *nz);
		printf("Number Of Variables: %i\n", *num_vars);
		printf("Number Of Spikes: %i\n", *num_spikes);
		printf("Number of Timesteps: %i\n", *num_tsteps);
		printf("Error Tolerance: %f\n", *err_tol);

		//Warum ist der erste Spike der gesetzt wird, nur fuer Variable 1 =  x*y*z???

		for(int k =0; k<(*num_spikes); k++){
			for(int l =0; l< (*num_vars); l++){
				std::cout << "spikes[" << k << "][" << l <<"]: " << std::setprecision (15) << spikes[(k * *num_vars) + l] << "\n";
			}
		}
		
		for(int k =0; k<(*num_spikes); k++){
			for(int l =0; l< 4; l++){
				std::cout << "loc[" << k << "][" << l <<"]: " << std::setprecision (15) << spike_loc[(k * 4) + l] << "\n";
			}
		}   
	}
*/	
/*    for( int j = 0 ; j < *num_vars ; ++j )
    {	
        std::cout << "source_total[" << j << "] = " << source_total[j] << "\n";
    }*/

    Typemaps::initializeMaps();
    {
		//SerialSimulator<Cell> sim(new CellInitializer(dimX, dimY, num_timesteps));
		CellInitializer *init =  new CellInitializer(*nx, *ny, (*num_tsteps * *num_spikes), *num_vars, source_total, spikes, spike_loc);
		
		//CheckerBoarding		
		HiParSimulator::HiParSimulator<Cell, RecursiveBisectionPartition<2> > sim(
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
