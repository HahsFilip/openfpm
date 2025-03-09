#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "sph.cuh"
template <int dim>
class ComputationalSpace{

    private:
            // static constexpr int dim = 2;
        using ParticleVector =  vector_dist_gpu<dim,double,aggregate<double,  double, double, double,Point<dim, double>, double , Point<dim, double>, double, Point<dim, double>,double, bool>>;
        // enum FIELDS  {

        // };

        double smoothing_distance;
        double extent;
        const double rho_zero = 6.5;
        const double coeff_sound = 2.0;
        const double gamma_ = 1.0;
        const double visco =0.9;
        const double p_zero = 10000;
        double Eta2;
        const double CFL_number = 0.2;
        const double MassFluid = 0.000614125;
        double B;
        double dt = 10e-4;
        double max_visc;
        const bool write_particles;
        size_t n;
        Box<dim, double> box;
        Ghost<dim, double> ghost;
        std::string average_vel_file_name;
        size_t bc[dim];

        void uniform_fill(){
            size_t dom[dim];
            double dp_tmp[dim]; 
            Box<dim, double> box_2 = box;

            for(int i = 0; i < dim; ++i){
                dom[i] = n;
                dp_tmp[i] =(box.getHigh(i) -box.getLow(i)) /(double)(n);
                std::cout << dp_tmp[i] << std::endl;
                box_2.setLow(i, box_2.getLow(i) + dp_tmp[i]/2.0);
                box_2.setHigh(i, box_2.getHigh(i) + dp_tmp[i]/2.0);

            }
            auto new_it = DrawParticles::DrawBox(particle_vec, dom, box, box_2);
            // auto it =  particle_vec.getGridIterator(dom);
            printf("\n");

            while (new_it.isNext())
            {
                particle_vec.add();
         
                auto key = new_it.get();
   

                for(int i = 0; i < dim; ++i){
                    particle_vec.getLastPos()[i] =new_it.get().get(i);//+dp_tmp[i]*0.05*(float)rand() / RAND_MAX;
                    // printf("%lf ", particle_vec.getLastPos()[i]);

                }   

                // next point
                ++new_it;
            }

            auto it2 = particle_vec.getDomainIterator();
            double i =0;
            while (it2.isNext()){
                auto p = it2.get();
                // printf("random: %f",0.1*(float)rand() / RAND_MAX);
                particle_vec.template getProp<SPH::FLUID_MASS>(p) =MassFluid;   
                particle_vec.template getProp<SPH::DENSITY>(p) = 0.0;
                for(int j =0; j< dim; ++j){
                    particle_vec.template getProp<SPH::VELOCITY>(p)[j] = 0.1*((float)rand()-0.5) / RAND_MAX;
                }
                particle_vec.template getProp<SPH::VELOCITY_OLD>(p) =Point<dim, double>(0.0);

                for(int j = 0; j < dim; j++){
                particle_vec.template getProp<SPH::D_V>(p)[j] = 0.0;
                particle_vec.template getProp<SPH::VELOCITY_OLD>(p)[j] = 0.0;


                }
                particle_vec.template getProp<SPH::DENSITY_OLD>(p) = 0.0;
                particle_vec.template getProp<SPH::D_RHO>(p) = 0.0;
                particle_vec.template getProp<SPH::REAL_ID>(p) = i;
                i += 1;

                ++it2;
            }
            particle_vec.map();
            particle_vec.template ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();

        }

 
    public: 
    
        ParticleVector particle_vec;
        ComputationalSpace(double extent, int n_val, double visc = 0.9, bool write_particles = true ):extent(extent), n(n_val), visco(visc), write_particles(write_particles){
            double zeros[dim];
            double extents[dim];
            printf("dp: %lf\n",(extent/(double)n_val) );
            smoothing_distance = sqrt(3.0*(extent/(double)n_val)*(extent/(double)n_val));
            double ghost_thickness = 2.0*smoothing_distance;
            printf("Smoothing distance: %lf\n", smoothing_distance);
            for(int i = 0; i < dim; ++i){
                zeros[i] = 0.0;
                extents[i] = extent;
                bc[i] = PERIODIC;
            }
            std::stringstream s_stream; 
            s_stream << "average_velocity_"<<visco<<".txt";
            average_vel_file_name = s_stream.str();
            B = coeff_sound*coeff_sound*rho_zero/gamma_;
            box = Box<dim, double>(zeros,extent); 
            std::cout << "box: " << box.getHigh(0) << " " << box.getLow(0) <<std::endl;
            ghost = Ghost<dim, double>(ghost_thickness);
            particle_vec = ParticleVector(0,box, bc, ghost);
            Eta2 = 0.01*smoothing_distance*smoothing_distance;
            uniform_fill();
	        openfpm::vector<std::string> names({"ID","MASS","SPH::DENSITY","SPH::PRESSURE","SPH::VELOCITY", "SPH::DENSITY_OLD", "SPH::VELOCITY_OLD", "SPH::D_RHO", "SPH::D_V", "SPH::REAL_ID"});
            particle_vec.setPropNames(names);
        	particle_vec.hostToDevicePos();
            particle_vec.template hostToDeviceProp<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();
            CalcDensity();

             
        }


        void CalcDensity(){
            // particle_vec.deleteGhost();
            particle_vec.map();
            particle_vec.template ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>(RUN_ON_DEVICE);

            auto NN = particle_vec.getCellListGPU(2.0*smoothing_distance);

            // particle_vec.ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V>();
            particle_vec.updateCellListGPU(NN);

            auto it = particle_vec.getDomainIteratorGPU();
      

            CUDA_LAUNCH(SPH::CalcDensityGPU, it, particle_vec.toKernel(), NN.toKernel(), smoothing_distance);
            particle_vec.deviceToHostPos();
            particle_vec.template deviceToHostProp<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();

        }

 

        

        void WriteParticles(int i){
            // particle_vec.template ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();
            // particle_vec.template ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();
            // particle_vec.template ghost_get<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();
            particle_vec.deviceToHostPos();
            particle_vec.template hostToDeviceProp<SPH::PARTICLE_ID,SPH::FLUID_MASS,SPH::DENSITY,SPH::PRESSURE,SPH::VELOCITY, SPH::DENSITY_OLD, SPH::VELOCITY_OLD, SPH::D_RHO, SPH::D_V, SPH::REAL_ID>();
            particle_vec.deleteGhost();
            // particle_vec.map();
            particle_vec.write_frame("particles", i);

        }

};

int main(int argc, char *argv[]){

    openfpm_init(&argc,&argv);
    const double dp = 0.0085;
    std::vector<double> viscocities = {0.1, 0.2, 0.4, 0.8, 0.9, 0.95, 0.99};
    ComputationalSpace<2> test(0.0085 * 250, 250, viscocities[0], false);
    test.WriteParticles(0);
    openfpm_finalize();
    return 0;
}
