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
            // static constexpr int dim = 2;
        using ParticleVector =  vector_dist<dim,double,aggregate<double,  double, double, double,Point<dim, double>, double , Point<dim, double>, double, Point<dim, double>,double>>;
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
        const double rho_zero = 9;
        const double coeff_sound = 2.0;
        const double gamma_ = 1.0;
        const double visco = 0.2;
        const double p_zero = 10000;
        double Eta2;
        const double CFL_number = 0.2;
        const double MassFluid = 0.000614125;
        double B;
        double dt = 10e-4;
        double max_visc;
        size_t n;
        Box<dim, double> box;
        Ghost<dim, double> ghost;

        size_t bc[dim];

        void uniform_fill(){
            size_t dom[dim];
            double dp_tmp[dim]; 
            for(int i = 0; i < dim; ++i){
                dom[i] = n;
                dp_tmp[i] =(box.getHigh(i) -box.getLow(i)) /n;
            }
            Box<dim, double> box_2 = box;
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
                    particle_vec.getLastPos()[i] =new_it.get().get(i)+dp_tmp[i]*0.05*(float)rand() / RAND_MAX;
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
                particle_vec.template getProp<VELOCITY>(p) = Point<dim, double>(0.0);//+ 0.01*(float)rand() / RAND_MAX;
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
        }
        // double spline_kernel(Point<dim, double> a, Point<dim, double> b){
        //     Point<dim, double> r_diff = b-a;
        //     double r = r_diff.norm();
        //     double x = r/smoothing_distance;
        //     // double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
        //     double multiplier = normalization() / std::pow(smoothing_distance, dim);
        //     double kernel_val;

        //     if( x <1){
        //         kernel_val = 1.0-1.5*x*x + 0.75*x*x*x;
        //     }else{
        //         if(x >=1 && x <2){
        //             kernel_val = (1.0/4.0)*(2.0-x)*(2.0-x)*(2.0-x);
        //         }else{
        //             kernel_val = 0;
        //         }
        //     }
        //     // if( x <1){
        //     //     kernel_val = 1-1.5*x*x + 0.75*x*x*x;
        //     // }else{
        //     //     if(x >=1 && x <2){
        //     //         kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
        //     //     }else{
        //     //         kernel_val = 0;
        //     //     }
        //     // }          
        //     // if( x <1){
        //     //     kernel_val = 1-1.5*x*x + 0.75*x*x*x;
        //     // }else{
        //     //     if(x >=1 && x <2){
        //     //         kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
        //     //     }else{
        //     //         kernel_val = 0;
        //     //     }
        //     // }
        //     return multiplier*kernel_val;
        // }

        // double spline_kernel(double r){
        //     double x = r/smoothing_distance;
        //     double multiplier = normalization() / std::pow(smoothing_distance, dim);
        //     // double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
        //     double kernel_val;
        //     // if( x <1){
        //     //     kernel_val = 1-1.5*x*x + 0.75*x*x*x;
        //     // }else{
        //     //     if(x >=1 && x <2){
        //     //         kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
        //     //     }else{
        //     //         kernel_val = 0;
        //     //     }
        //     // }
        //     if( x <1){
        //         kernel_val = 1.0-1.5*x*x + 0.75*x*x*x;
        //     }else{
        //         if(x >=1 && x <2){
        //             kernel_val = (1.0/4.0)*(2.0-x)*(2.0-x)*(2.0-x);
        //         }else{
        //             kernel_val = 0;
        //         }
        //     }
        //     return multiplier*kernel_val;
        // }

        // Point<dim, double> grad_spline_kernel(Point<dim, double> a, Point<dim, double> b){
        //     Point<dim, double> tmp = b-a;
        //     double r = tmp.norm();
        //     Point<dim, double> direction = tmp/tmp.norm();
        //     double x = r/smoothing_distance;
        //     // double multipler = 1.0/(M_PI*smothing_distance*smoothing_distance*smoothing_distance*smoothing_distance);
        //     double multiplier = normalization() / std::pow(smoothing_distance, dim + 1);
        //     double kernel_val;
        //     if(x>=0 && x <=1){
        //         kernel_val = (9.0/4.0)*x*x -3.0*x;
        //     }else{
        //         if(x >=1 && x <=2){
        //             kernel_val = -(3.0/4.0)*(2.0-x)*(2.0-x);
        //         }else{
        //             kernel_val = 0;
        //         }
        //     }

        //     return direction*multiplier*kernel_val;
        // }
        double spline_kernel(Point<dim, double> a, Point<dim, double> b){
            Point<dim, double> r_diff = b-a;
            double r = r_diff.norm();
            double x = r/smoothing_distance;
            double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
            double kernel_val;
            if( x <1){
                kernel_val = 1-1.5*x*x + 0.75*x*x*x;
            }else{
                if(x >=1 && x <2){
                    kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }
            return multipler*kernel_val;
        }

        double spline_kernel(double r){
            double x = r/smoothing_distance;
            double multipler = 1.0/(M_PI*smoothing_distance*smoothing_distance*smoothing_distance);
            double kernel_val;
            if( x <1){
                kernel_val = 1-1.5*x*x + 0.75*x*x*x;
            }else{
                if(x >=1 && x <2){
                    kernel_val = 0.25*(2.0-x)*(2.0-x)*(2.0-x);
                }else{
                    kernel_val = 0;
                }
            }
            return multipler*kernel_val;
        }

        Point<dim, double> grad_spline_kernel(Point<dim, double> a, Point<dim, double> b){
            Point<dim, double> tmp = b-a;
            double r = tmp.norm();
            Point<dim, double> direction = tmp/tmp.norm();
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
        double dot(Point<dim, double> a, Point<dim, double> b){
            double result = 0;
            for(int i = 0; i < dim; ++i){
                result += a[i]*b[i];
            }
            return result;
        }

        double StateEquation(double rho){
            return B*(pow(rho/rho_zero, gamma_)-1.0);
        }
        void CalcPressure(){
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            auto it = particle_vec.getDomainIterator();

            while(it.isNext()){ 
                auto a = it.get();
                double tmp_rho = particle_vec.template getProp<DENSITY>(a);

                particle_vec.template getProp<PRESSURE>(a) = StateEquation(tmp_rho);
                ++it;
            }
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
        }
        double Pi(const Point<dim, double> dr, const Point<dim, double> dv, double rho_a, double rho_b, double mass_b ){
            double dot_result = dot(dr, dv);
            double dot_rr2 = dot_result/(dr.norm()+Eta2);
            max_visc = std::max(dot_rr2, max_visc);
            if(dot_result < 0){
                const float amubar=smoothing_distance*dot_rr2;
                const float robar=(rho_a+rho_b)*0.5f;
                const double pi_visc=(-visco*coeff_sound*amubar/robar);
        		return pi_visc;
            }else{
                return 0.0;
            }
        }

 
    public: 
    
        ParticleVector particle_vec;
        ComputationalSpace(double extent, int n_val):extent(extent), n(n_val){
            double zeros[dim];
            double extents[dim];
            printf("dp: %lf",(extent/(double)n_val) );
            smoothing_distance = sqrt(3.0*(extent/(double)n_val)*(extent/(double)n_val));
            double ghost_thickness = 2.0*smoothing_distance;
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
            Eta2 = 0.01*smoothing_distance*smoothing_distance;
            uniform_fill();
	        openfpm::vector<std::string> names({"ID","MASS","DENSITY","PRESSURE","VELOCITY", "DENSITY_OLD", "VELOCITY_OLD", "D_RHO", "D_V", "REAL_ID"});
            particle_vec.setPropNames(names);
            CalcDensity();
            CalcPressure();
            initializeOldFields();
            
        }

        void calculate_timestep(){
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            auto it = particle_vec.getDomainIterator();
            double v_max = 0;
            double dv_max = 0;

            while(it.isNext()){ 
                auto a = it.get();
                Point<dim, double>  tmp_v = particle_vec.template getProp<VELOCITY>(a);
                Point<dim, double>  tmp_a = particle_vec.template getProp<D_V>(a);
                if(tmp_v.norm() > v_max){
                    v_max = tmp_v.norm();
                }
                if(tmp_a.norm() > dv_max){
                    dv_max = tmp_a.norm();
                }
                ++it;
            }
            Vcluster<> & v_cl = create_vcluster();
            v_cl.max(v_max);
            v_cl.max(dv_max);
            v_cl.max(max_visc);
            v_cl.execute();
        	double dt_f = (dv_max)?sqrt(smoothing_distance/dv_max):std::numeric_limits<int>::max();
        	const double dt_cv = smoothing_distance/(std::max(coeff_sound,dv_max*10.) + smoothing_distance*max_visc);

            std::cout << "V_MAX: " << v_max  << std::endl;
            std::cout << "max_visc: " << max_visc  << std::endl;
            std::cout << "dt_f: " << dt_f  << std::endl;
            std::cout << "dt_cv: " << dt_cv  << std::endl;
            dt=(CFL_number)*std::min(dt_f,dt_cv);
            std::cout << "dt: " << dt  << std::endl;


        }
        double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
        {
            const double qq=r/smoothing_distance;
            //-Cubic Spline kernel
            double wab;
            const double W_dap = 1.0/spline_kernel(smoothing_distance/1.5);
            const double a2 = 1.0/M_PI/smoothing_distance/smoothing_distance/smoothing_distance;
            const double a2_4 = 0.25*a2;


            if(r>smoothing_distance)
            {
                double wqq1=2.0f-qq;
                double wqq2=wqq1*wqq1;
        
                wab=a2_4*(wqq2*wqq1);
            }
            else
            {
                double wqq2=qq*qq;
                double wqq3=wqq2*qq;
        
                wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
            }
        
            //-Tensile correction.
            double fab=wab*W_dap;
            fab*=fab; fab*=fab; //fab=fab^4
            const double tensilp1=(prs1/(rhoa*rhoa))*(prs1>0? 0.01: -0.2);
            const double tensilp2=(prs2/(rhob*rhob))*(prs2>0? 0.01: -0.2);
        
            return (fab*(tensilp1+tensilp2));
        }
        
        void UpdateOld(){
                particle_vec.map();
                particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
                auto it = particle_vec.getDomainIterator();
    
                while(it.isNext()){ 
                    auto a = it.get();
                    double tmp_rho = particle_vec.template getProp<DENSITY>(a);
    
                    particle_vec.template getProp<VELOCITY_OLD>(a) =  particle_vec.template getProp<VELOCITY>(a);
                    particle_vec.template getProp<DENSITY_OLD>(a) =  particle_vec.template getProp<DENSITY>(a);

                    ++it;
                }
        }

        void CalcForces(){
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

            auto NN = particle_vec.getCellList(2.0*smoothing_distance);

            // particle_vec.ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V>();
            particle_vec.updateCellList(NN);

            auto it = particle_vec.getDomainIterator();

            while (it.isNext()){
                auto p = it.get();
                auto Np = NN.getNNIteratorBox(NN.getCell(particle_vec.getPos(p)));
                particle_vec.template getProp<D_V>(p) = particle_vec.template getProp<D_V>(p)*0.0;
                Point<dim, double> pos_1 = particle_vec.getPos(p);

                Point<dim, double> vel_1 = particle_vec.template getProp<VELOCITY>(p);
                double rho_1 = particle_vec.template getProp<DENSITY>(p);
                double p_1 = particle_vec.template getProp<PRESSURE>(p);
                Point<dim, double> tmp_dv = 0.0;
                while(Np.isNext()==true){
                    auto np = Np.get();

                    if(p.getKey() == np){++Np; continue;}
        
                    Point<dim, double> pos_2 = particle_vec.getPos(np);
                    Point<dim, double> vel_2 = particle_vec.template getProp<VELOCITY>(np);
                    double rho_2 = particle_vec.template getProp<DENSITY>(np);
                    double p_2 = particle_vec.template getProp<PRESSURE>(np);

                    Point<dim, double> r = pos_1-pos_2;
                    Point<dim, double> dvel = vel_1-vel_2;
        
                    double local_mass = particle_vec.template getProp<FLUID_MASS>(np);
                    double pi_val = Pi(r,dvel, rho_1, rho_2, local_mass);
                    //  tmp_dv =  particle_vec.template getProp<D_V>(p);
                    double tensile_val = Tensile(r.norm(),rho_1,rho_2,p_1,p_2);
                    Point<dim, double> tmp_diff= local_mass*((p_1+p_2)/(rho_1*rho_2) + pi_val+tensile_val)*grad_spline_kernel(pos_1, pos_2);
                    // if(tmp_diff.norm() >100){
                    //     std::cout <<"tmp_diff: " <<  tmp_diff.norm() <<std::endl;
                    // }

                    for(int i = 0; i<dim;++i){
                        particle_vec.template getProp<D_V>(p)[i] +=  tmp_diff[i];
                    } 
                    // particle_vec.template getProp<D_V>(p) -=  local_mass*((p_1+p_2)/(rho_1*rho_2))*grad_spline_kernel(pos_1, pos_2);
                    ++Np;

                }
                // std::cout <<  tmp_dv_outer <<std::endl;
                // std::cout <<"tmp_dv: " <<  tmp_dv.norm() <<std::endl;
        
                ++it;
            }
            particle_vec.map();
        }

// Add this function inside your ComputationalSpace class.
    void initializeOldFields() {
        // Ensure the initial density derivative and forces have been computed.
        CalcDRho();
        CalcForces();
        
        auto it = particle_vec.getDomainIterator();
        while (it.isNext()) {
            auto p = it.get();
            // "Back-step" the velocity and density
            particle_vec.template getProp<VELOCITY_OLD>(p) =
                particle_vec.template getProp<VELOCITY>(p) - dt * particle_vec.template getProp<D_V>(p);
            particle_vec.template getProp<DENSITY_OLD>(p) =
                particle_vec.template getProp<DENSITY>(p) - dt * particle_vec.template getProp<D_RHO>(p);
            ++it;
        }
    }

        void CalcDensity(){

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

                    if(p.getKey() == np){++Np; continue;}
        
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

        void CalcDRho(){
            particle_vec.map();
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

            auto NN = particle_vec.getCellList(2.0*smoothing_distance);

            // particle_vec.ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V>();
            particle_vec.updateCellList(NN);

            auto it = particle_vec.getDomainIterator();

            while (it.isNext()){
                auto p = it.get();
                auto Np = NN.getNNIteratorBox(NN.getCell(particle_vec.getPos(p)));
                particle_vec.template getProp<D_RHO>(p) = particle_vec.template getProp<D_RHO>(p)*0.0;

                Point<dim, double> pos_1 = particle_vec.getPos(p);
                Point<dim, double> vel_1 = particle_vec.template getProp<VELOCITY>(p);
                
                while(Np.isNext()==true){
                    auto np = Np.get();

                    if(p.getKey() == np){++Np; continue;}
        
                    Point<dim, double> pos_2 = particle_vec.getPos(np);
                    Point<dim, double> vel_2 = particle_vec.template getProp<VELOCITY>(np);
                    Point<dim, double> rel_vel = vel_1-vel_2;
                    double local_mass = particle_vec.template getProp<FLUID_MASS>(np);

                    particle_vec.template getProp<D_RHO>(p) +=  local_mass*dot(grad_spline_kernel(pos_2, pos_1), rel_vel);
                    ++Np;
                }
        
                ++it;
            }
        }
        void VerletTime(int n){

            for(int i = 0; i<n; ++i){
                particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

                particle_vec.map();
                std::cout << "TIMESTEP: " << i << std::endl;
                max_visc = 0.0;\
                if(i%1== 0){
                    WriteParticles(i);
                }
                CalcDRho();
                CalcPressure();
                CalcForces();
                calculate_timestep();
                auto it = particle_vec.getDomainIterator();
    
                while(it.isNext()){ 
                    auto a = it.get();
                    Point<dim, double> v_plus_one = particle_vec.template getProp<VELOCITY_OLD>(a) + 2.0*dt*particle_vec.template getProp<D_V>(a); 
                    Point<dim, double> r_plus_1 = particle_vec.getPos(a);
                    r_plus_1 += dt*particle_vec.template getProp<VELOCITY>(a) + 0.5*dt*dt*particle_vec.template getProp<D_V>(a);

                    for(int j = 0; j< dim; ++j){
                        particle_vec.getPos(a)[j] = r_plus_1[j]; 
                    }
                    double rho_plus_one = particle_vec.template getProp<DENSITY_OLD>(a) + 2.0*dt*particle_vec.template getProp<D_RHO>(a);
                    particle_vec.template getProp<DENSITY_OLD>(a) = particle_vec.template getProp<DENSITY>(a);
                    particle_vec.template getProp<VELOCITY_OLD>(a) = particle_vec.template getProp<VELOCITY>(a);
                    particle_vec.template getProp<VELOCITY>(a) = v_plus_one;
                    particle_vec.template getProp<DENSITY>(a) = rho_plus_one;

                    // double tmp_rho = particle_vec.template getProp<DENSITY>(a);
    
                    // particle_vec.template getProp<VELOCITY_OLD>(a) =  particle_vec.template getProp<VELOCITY>(a);
                    // particle_vec.template getProp<DENSITY_OLD>(a) =  particle_vec.template getProp<DENSITY>(a);

                    ++it;
                }
                particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
                particle_vec.map();


            }
        
        }


        

        void WriteParticles(int i){
            particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            particle_vec.map();
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();
            // particle_vec.template ghost_get<PARTICLE_ID,FLUID_MASS,DENSITY,PRESSURE,VELOCITY, DENSITY_OLD, VELOCITY_OLD, D_RHO, D_V, REAL_ID>();

            // particle_vec.map();
            particle_vec.write_frame("particles", i);

        }

};

int main(int argc, char *argv[]){

    openfpm_init(&argc,&argv);
    const double dp = 0.0085;

        // ComputationalSpace<1> test(1, 10, 2); 
    ComputationalSpace<2> test_2(0.0085*50, 50); 
    // ComputationalSpace<2> test_2(1, 100); 
    // test_2.CalcDRho();
    // test_2.CalcForces();
    // test_2.CalcDensity();
    test_2.WriteParticles(0);
    test_2.VerletTime(50000);
	openfpm_finalize();

    return 0;
}