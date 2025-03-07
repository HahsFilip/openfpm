#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

template <int dim>
class ComputationalSpace{

    private:
            // static constexpr int dim = 2;
        using ParticleVector =  vector_dist<dim,double,aggregate<double,  double, double, double,Point<dim, double>, double , Point<dim, double>, double, Point<dim, double>,double, bool>>;
        // enum FIELDS  {
            static constexpr int PARTICLE_ID = 0;
            static constexpr int FLUID_MASS = 1;
            static constexpr int DENSITY = 2;
            static constexpr int PRESSURE = 3;
            static constexpr int VELOCITY = 4;
            static constexpr int DENSITY_OLD = 5;
            static constexpr int VELOCITY_OLD = 6;
            static constexpr int D_RHO = 7; 
            static constexpr int D_V = 8;
            static constexpr int REAL_ID = 9;
            static constexpr int WALL = 10;
        // };
        constexpr double normalization() const {
            if constexpr(dim == 1)
                return 2.0/3.0;
            else if constexpr(dim == 2)
                return 10.0 / (7.0 * M_PI);
            else if constexpr(dim == 3)
                return 1.0 / (1.0*M_PI);
            else
                return 1.0;  // Fallback (this may not be correct for higher dimensions)
        }
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
         
                // The position of the particle is given by the point-index (0,0,0) ; (0,0,1) ... (127,127,127)
                // multiplied by the spacing

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
                particle_vec.template getProp<FLUID_MASS>(p) =MassFluid;   
                particle_vec.template getProp<DENSITY>(p) = 0.0;
                for(int i =0; i < dim; ++i){
                    particle_vec.template getProp<VELOCITY>(p)[i] = 0.1*((float)rand()-0.5) / RAND_MAX;
                }
                particle_vec.template getProp<VELOCITY_OLD>(p) =Point<dim, double>(0.0);

                for(int j = 0; j < dim; j++){
                    // particle_vec.template getProp<VELOCITY>(p)[j] +=0.001*(float)rand() / RAND_MAX;
                particle_vec.template getProp<D_V>(p)[j] = 0.0;
                particle_vec.template getProp<VELOCITY_OLD>(p)[j] = 0.0;


                }
                particle_vec.template getProp<DENSITY_OLD>(p) = 0.0;
                particle_vec.template getProp<D_RHO>(p) = 0.0;
                particle_vec.template getProp<REAL_ID>(p) = i;
                i += 1;

                ++it2;
            }
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

        }
        __device__ __host__ double spline_kernel(Point<dim, double> a, Point<dim, double> b){
            Point<dim, double> r_diff = b-a;
            double r = r_diff.norm();
            double x = r/smoothing_distance;
            // double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
            double multiplier = normalization() / std::pow(smoothing_distance, dim);
            double kernel_val;

            if( x <1){
                kernel_val = 1.0-1.5*x*x + 0.75*x*x*x;
            }else{
                if(x >=1 && x <2){
                    kernel_val = (1.0/4.0)*(2.0-x)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }
  
            return multiplier*kernel_val;
        }

        __device__ __host__ double spline_kernel(double r){
            double x = r/smoothing_distance;
            double multiplier = normalization() / std::pow(smoothing_distance, dim);
            double kernel_val;
            if( x <1){
                kernel_val = 1.0-1.5*x*x + 0.75*x*x*x;
            }else{
                if(x >=1 && x <2){
                    kernel_val = (1.0/4.0)*(2.0-x)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }
            return multiplier*kernel_val;
        }

        __device__ __host__ Point<dim, double> grad_spline_kernel(Point<dim, double> a, Point<dim, double> b){
            Point<dim, double> tmp = b-a;
            double r = tmp.norm();
            Point<dim, double> direction = tmp/tmp.norm();
            double x = r/smoothing_distance;
            // double multipler = 1.0/(M_PI*smothing_distance*smoothing_distance*smoothing_distance*smoothing_distance);
            double multiplier = normalization() / std::pow(smoothing_distance, dim + 1);
            double kernel_val;
            if(x>=0 && x <=1){
                kernel_val = (9.0/4.0)*x*x -3.0*x;
            }else{
                if(x >=1 && x <=2){
                    kernel_val = -(3.0/4.0)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }

            return direction*multiplier*kernel_val;
        }
       
        __device__ __host__ double dot(Point<dim, double> a, Point<dim, double> b){
            double result = 0;
            for(int i = 0; i < dim; ++i){
                result += a[i]*b[i];
            }
            return result;
        }

        __device__ __host__ double StateEquation(double rho){
            return B*(pow(rho/rho_zero, gamma_)-1.0);
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
	        openfpm::vector<std::string> names({"ID","MASS","DENSITY","PRESSURE","VELOCITY", "DENSITY_OLD", "VELOCITY_OLD", "D_RHO", "D_V", "REAL_ID"});
            particle_vec.setPropNames(names);
            CalcDensity();
            
        }


        void CalcDensity(){
            // particle_vec.deleteGhost();
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

            auto NN = particle_vec.getCellList(2.0*smoothing_distance);

            // particle_vec.ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V>();
            particle_vec.updateCellList(NN);

            auto it = particle_vec.getDomainIterator();

            while (it.isNext()){
                auto p = it.get();
                auto Np = NN.getNNIteratorBox(NN.getCell(particle_vec.getPos(p)));
                particle_vec.template getProp<DENSITY>(p) = 0.0;
                Point<dim, double> pos_1 = particle_vec.getPos(p);
        
                while(Np.isNext()==true){
                    auto np = Np.get();

                    if(p.getKey() == np ){++Np; continue;}
        
                    Point<dim, double> pos_2 = particle_vec.getPos(np);
                    Point<dim, double> r = pos_2-pos_1;
           
                    double local_mass = particle_vec.template getProp<FLUID_MASS>(np);
                    particle_vec.template getProp<DENSITY>(p) +=  local_mass*spline_kernel(pos_1, pos_2);
                    ++Np;
                }
        
                ++it;
            }
            particle_vec.map();
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
        
        }

 

        

        void WriteParticles(int i){
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            particle_vec.deleteGhost();
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

            // particle_vec.map();
            particle_vec.write_frame("particles", i);

        }

};

int main(int argc, char *argv[]){

    openfpm_init(&argc,&argv);
    const double dp = 0.0085;
    std::vector<double> viscocities = {0.1,0.2,0.4,0.8,0.9,0.95,0.99};
    ComputationalSpace<2> test(0.0085*250, 250, viscocities[0], false); 
    test.WriteParticles(0);
        // ComputationalSpace<1> test(1, 10, 2); 
    for(int i = 0; i < viscocities.size(); ++i){
        ComputationalSpace<2> test_2(0.0085*150, 150, viscocities[i], false); 
        // test_2.Run(100000, 10.0);

    }    


    openfpm_finalize();

    return 0;
}