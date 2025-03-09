
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
namespace SPH{
    constexpr int PARTICLE_ID = 0;
    constexpr int FLUID_MASS = 1;
    constexpr int DENSITY = 2;
    constexpr int PRESSURE = 3;
    constexpr int VELOCITY = 4;
    constexpr int DENSITY_OLD = 5;
    constexpr int VELOCITY_OLD = 6;
    constexpr int D_RHO = 7; 
    constexpr int D_V = 8;
    constexpr int REAL_ID = 9;
    constexpr int WALL = 10;
    template<int dim>
    constexpr double normalization()  {
        if constexpr(dim == 1)
            return 2.0/3.0;
        else if constexpr(dim == 2)
            return 10.0 / (7.0 * M_PI);
        else if constexpr(dim == 3)
            return 1.0 / (1.0*M_PI);
        else
            return 1.0;  // Fallback (this may not be correct for higher dimensions)
    }

    template<int dim>
    __device__ double spline_kernel(Point<dim, double> a, Point<dim, double> b, double smoothing_distance){
        Point<dim, double> r_diff = b-a;
        double r = r_diff.norm();
        double x = r/smoothing_distance;
        // double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
        double multiplier = SPH::normalization<dim>() / std::pow(smoothing_distance, dim);
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

    template<int dim>
    __device__ double spline_kernel(double r, double smoothing_distance){
        double x = r/smoothing_distance;
        double multiplier = normalization<dim>() / std::pow(smoothing_distance, dim);
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
    template<int dim>
    __device__ Point<dim, double> grad_spline_kernel(Point<dim, double> a, Point<dim, double> b, double smoothing_distance){
        Point<dim, double> tmp = b-a;
        double r = tmp.norm();
        Point<dim, double> direction = tmp/tmp.norm();
        double x = r/smoothing_distance;
        // double multipler = 1.0/(M_PI*smothing_distance*smoothing_distance*smoothing_distance*smoothing_distance);
        double multiplier = SPH::normalization<dim>() / std::pow(smoothing_distance, dim + 1);
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
    template<int dim>
    __device__ double dot(Point<dim, double> a, Point<dim, double> b){
        double result = 0;
        for(int i = 0; i < dim; ++i){
            result += a[i]*b[i];
        }
        return result;
    }
   
    template<int dim>
    __device__ double StateEquation(double rho, double rho_zero, double gamma_, double B){
        return B*(pow(rho/rho_zero, gamma_)-1.0);
    }
    
    template< typename vec_type, typename NN_type>
    __global__ void CalcDensityGPU(vec_type particles, NN_type NN, double smoothing_distance){
        constexpr int dim = decltype(particles)::dims; // deduce dimension from the particle vector

        auto p = GET_PARTICLE(particles);
        auto Np = NN.getNNIteratorBox(NN.getCell(particles.getPos(p)));
        particles.template getProp<SPH::DENSITY>(p) = 0.0;
        Point<dim, double> pos_1 = particles.getPos(p);

        while(Np.isNext()){
            auto np = Np.get();


            Point<dim, double> pos_2 = particles.getPos(np);
            Point<dim, double> r = pos_2-pos_1;
            // printf("spline_kernel: %lf\n", spline_kernel<dim>(pos_1, pos_2, smoothing_distance));

            double local_mass = particles.template getProp<SPH::FLUID_MASS>(np);
            particles.template getProp<SPH::DENSITY>(p) +=  local_mass*spline_kernel<dim>(pos_1, pos_2, smoothing_distance);
            ++Np;
        }

    }
    template< typename vec_type, typename NN_type>
    __global__ void CalcPressureGPU(vec_type particles, NN_type NN, double smoothing_distance, double rho_zero, double gamma_, double B){
        auto a = GET_PARTICLE(particles);

            double tmp_rho = particles.template getProp<DENSITY>(a);
            double pressure =  StateEquation(tmp_rho, rho_zero, gamma_, B);
            if(pressure < 0){
                std::cout << 'PROBLEM, PRESSURE < 0 ' << std::endl;
            }
            particles.template getProp<PRESSURE>(a) = pressure;

        particles.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
    }

}