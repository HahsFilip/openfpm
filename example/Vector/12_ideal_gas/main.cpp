#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#define DIM 1
#undef MASS
#undef DENSITY
#undef PRESSURE
#undef VELOCITY
#undef VELOCITY_OLD
#undef D_V


template <int dim>
class ComputationalSpace{

    private:
        using ParticleVector =  vector_dist<dim,double,aggregate<size_t,  double, double, double,Point<dim, double>, double , Point<dim, double>, double, Point<dim, double>>>;
        enum FIELDS  {
            PARTICLE_ID = 0,
            FLUID_MASS = 1,
            DENSITY = 2,
            PRESSURE = 3,
            VELOCITY = 4,
            DENSITY_OLD = 5,
            VELOCITY_OLD = 6,
            D_RHO = 7, 
            D_V = 8
        };
   
        double smoothing_distance;
        double extent;
        const double rho_zero = 1000.0;
        const double coeff_sound = 20.0;
        const double gamma_ = 7.0;
        const double visco = 0.1;
        double Eta2;
        const double MassFluid = 0.000614125;
        double B;


        size_t n;
        Box<dim, double> box;
        Ghost<dim, double> ghost;

        size_t bc[dim];

        void uniform_fill(){
            size_t dom[dim];
            for(int i = 0; i < dim; ++i){
                dom[i] = n;
            }
            auto new_it = DrawParticles::DrawBox(particle_vec, dom, box, box);
            // auto it =  particle_vec.getGridIterator(dom);
            printf("\n");

            while (new_it.isNext())
            {
                particle_vec.add();
         
                auto key = new_it.get();
         
                // The position of the particle is given by the point-index (0,0,0) ; (0,0,1) ... (127,127,127)
                // multiplied by the spacing

                for(int i = 0; i < dim; ++i){
                    particle_vec.getLastPos()[i] =new_it.get().get(i);//+it.getSpacing(i)*0.05*(float)rand() / RAND_MAX;
                    printf("%lf ", particle_vec.getLastPos()[i]);

                }   
                printf("\n");

                // next point
                ++new_it;
            }
            auto it2 = particle_vec.getDomainIterator();
            while (it2.isNext()){
                auto p = it2.get();
                // printf("random: %f",0.1*(float)rand() / RAND_MAX);
                particle_vec.template getProp<FIELDS::FLUID_MASS>(p) =MassFluid;   
                particle_vec.template getProp<FIELDS::DENSITY>(p) = 0.0;
                particle_vec.template getProp<FIELDS::VELOCITY>(p) = Point<dim, double>(0.0);//+ 0.01*(float)rand() / RAND_MAX;
                particle_vec.template getProp<FIELDS::DENSITY_OLD>(p) = 0.0;
                particle_vec.template getProp<FIELDS::VELOCITY_OLD>(p) =Point<dim, double>(0.0);
                particle_vec.template getProp<FIELDS::D_RHO>(p) = 0.0;
                particle_vec.template getProp<FIELDS::D_V>(p) = Point<dim, double>(0.0);
                ++it2;
            }
        }
        double spline_kernel(Point<dim, double> a, Point<dim, double> b){
            Point<dim, double> r_diff = b-a;
            double r = r_diff.norm();
            double x = r/smoothing_distance;
            double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
            double kernel_val;
            if( x <1){
                kernel_val = 1-1.5*x*x + 0.75*x*x*x;
            }else{
                if(x >=1 && x <=2){
                    kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }
            return multipler*kernel_val;
        }
        Point<dim, double> grad_spline_kernel(Point<dim, double> a, Point<dim, double> b){
            double r = (a-b).norm();
            Point<dim, double> direction = (b-a).normalized();
            double x = r/smoothing_distance;
            double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance*smoothing_distance);
            double kernel_val;
            if(x>=0 && x <=1){
                kernel_val = (9.0/4.0)*x*x -3.0*x;
            }else{
                if(x >=1 && x <=2){
                    kernel_val = -0.75*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }

            return direction*multipler*kernel_val;
        }
        double StateEquation(double rho){
            return B*(pow(rho/rho_zero, gamma_)-1.0);
        }
        void CalcPressure(){
            particle_vec.map();
            particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();
            auto it = particle_vec.getDomainIterator();

            while(it.isNext()){
                auto a = it.get();
                double tmp_rho = particle_vec.template getProp<FIELDS::DENSITY>(a);

                particle_vec.template getProp<FIELDS::PRESSURE>(a) = StateEquation(tmp_rho);
                ++it;
            }
            particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();


        }
    
    public: 
    
        ParticleVector particle_vec;
        ComputationalSpace(double extent, int n_val, double ghost_thickness):extent(extent), n(n_val){
            double zeros[dim];
            double extents[dim];
            smoothing_distance = sqrt(3.0*(extent/n));
            printf("Smoothing distance: %lf", smoothing_distance);
            for(int i = 0; i < dim; ++i){
                zeros[i] = 0.0;
                extents[i] = extent;
                bc[i] = PERIODIC;
            }

            B = coeff_sound*coeff_sound*rho_zero/gamma_;
            box = Box<dim, double>(zeros,extent); 
            ghost = Ghost<dim, double>(ghost_thickness);
            particle_vec = ParticleVector(0,box, bc, ghost);
            uniform_fill();
	        openfpm::vector<std::string> names({"ID","MASS","DENSITY","PRESSURE","VELOCITY", "DENSITY_OLD", "VELOCITY_OLD", "D_RHO", "D_V"});
            particle_vec.setPropNames(names);
            CalcDensity();
            CalcPressure();
            
        }

        void CalcDensity(){
            // particle_vec.ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();
            particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();

            auto NN = particle_vec.getCellList(10.0*smoothing_distance);
            auto it = particle_vec.getDomainIterator();
            while (it.isNext()){
                auto p = it.get();
                auto Np = NN.getNNIteratorBox(NN.getCell(particle_vec.getPos(p)));
                particle_vec.template getProp<FIELDS::DENSITY>(p) = 0.0;
        
                while(Np.isNext()){
                    auto np = Np.get();
                    if(p.getKey() == np){++Np; continue;}
        
                    Point<dim, double> pos_1 = particle_vec.getPos(p);
                    Point<dim, double> pos_2 = particle_vec.getPos(np);
                    Point<dim, double> r = pos_1-pos_2;
                    // printf("r norm: %lf\n",r.norm());
                    if(r.norm()< 10e-10){
                        ++Np;
                        continue;
                    }
                    double local_mass = particle_vec.template getProp<FIELDS::FLUID_MASS>(np);
                    particle_vec.template getProp<FIELDS::DENSITY>(p) +=  local_mass*spline_kernel(pos_1, pos_2);
                    ++Np;
                }
                // printf("density end: %lf\n",vec.getProp<FIELDS::DENSITY>(p));
        
                ++it;
            }
            particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();
        
        }




        

        void WriteParticles(int i){
        
            // particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();
            // particle_vec.template ghost_get<FIELDS::PARTICLE_ID,FIELDS::FLUID_MASS,FIELDS::DENSITY,FIELDS::PRESSURE,FIELDS::VELOCITY, FIELDS::DENSITY_OLD, FIELDS::VELOCITY_OLD, FIELDS::D_RHO, FIELDS::D_V>();

            // particle_vec.map();k
            printf("Smoothing distance: %lf", smoothing_distance);
            particle_vec.template write("particles", BINARY);
        }

};

int main(int argc, char *argv[]){

    openfpm_init(&argc,&argv);
    const double dp = 0.0085;

    // ComputationalSpace<1> test(1, 10, 2); 
    // ComputationalSpace<2> test_2(0.085*0.5, 0.085*0.5/dp, 2*dp); 
    ComputationalSpace<1> test_2(1, 1000, 0.2); 
    
    // test_2.CalcDensity();
    test_2.WriteParticles(0);
	openfpm_finalize();

    return 0;
}