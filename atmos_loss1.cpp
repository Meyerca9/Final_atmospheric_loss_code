#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;
int main(){
    double m1,m2,m_imp,v_imp,rho,h,R,impactor_bulk_density,r,pi,M1,r_cap,M1cap,
    new_atmos,cap_factor,m_min,m_max,r_min,m_eject,theta,M2,M_solid,M_atmos,
    x,X_loss,v_esc,Mt,mu,k,T,G,g,p_base,A,V,Rd,N,r_gi,sum,q,m_cap;


    impactor_bulk_density =  2000; // kg/m^3
    pi =  3.1415926535897;
    mu = 4.8e-23; // g mean molecular weight of atmosphere
    k = 1.38e-23; // J/K boltzmann constant
    G = 6.674e-11; // Gravitational constant
    Rd = 287.058; // Gas constant in J/(kg*K)
    // M1 = solid mass of target (kg)
    // m_imp = solid mass of projectile
    // m1 = Mass of the targets Atmosphere (kg)
    // m2 = Mass of the projectiles atmosphere (kg)
    // v_imp = impact velocity (m/s)
    // theta = impact angle (enter radians with decimal)
    // M_lr = largest remnant mass
    // M_slr = second largest remnant mass
    //M_frag = cumulative mass of fragmented objects
    //v_sc = super-catastrophic critical velocity
    
    sum = 0;
    for (int i = 0; i <=10000; ++i){

    std::cout <<"Enter the target mass";
    std::cin >> M1;

    std::cout <<"Enter the projectiles mass";
    std::cin >> m_imp;

    std::cout <<"Enter the Mass of the targets atmosphere";
    std::cin >> m1;

    std::cout <<"Enter the impact velocity";
    std::cin >> v_imp;
  
    std::cout <<"Enter the impact angle";
    std::cin >> theta;

    std::cout <<"Enter the temperature at the base of the target atmosphere";
    std::cin >> T;

    std::cout <<"Enter the type of atmosphere: \n Adiabatic: 1\n Isothermal: 2\n";
    char answer;
    std::cin >> answer;
    
    

    // M1cap = targets cap mass
    // M2cap = projectiles cap mass 
    // v_esc = escape velocity (m/s)
    // R = target radius (km)
    // r = impactor radius (km)
    // m_eject = mass that is ejected above the tangent line
    // X_loss = mass loss due to velocity of impactor
    // Mt = Total impactor mass
    // g = surface gravity 
    // V = volume of earth
    // A = surface area of target



    R = std::cbrt(M1/(1.33333*pi*5510)); // radius of target in km (5510 is impactor bulk density of earth) in meters
    r = std::cbrt(m_imp/(1.33333*pi*impactor_bulk_density)); // radius of projectile in meters
    g = (G*M1)/(pow(R,2)); // surface gravity 
    A = 4*pi*pow(R,2); // surface area of target m^2 ------
    p_base = (m1*g)/A; //pressure at bottom of an atmosphere in Pa
    rho = p_base/(Rd*T); // atmospheric density at the surface
    h = k*T/(mu*g)*1e3; // atmospheric scale height in km
    v_esc = sqrt((2*G*M1)/R)/1000; // escape velocity in km/s
    cap_factor = sqrt((pi*R)/(2*h))*M1cap; // Impactor mass needed to eject all the atmospheric mass for planetesimals
    M1cap = ((2*pi)*(rho)*pow(h,2)*(R)); // maximum atmospheric mass that a single planetesimal impact can eject above the tangent plane
    m_min = 4*pi*rho*pow(h,3); // minimum mass needed in order to eject any atmosphere
    m_cap = pow(pi*h/(8*R), 1.0 / 2.0)*m1; // masses greater than this are able to eject all mass above the tangent line which is h/2r of the atmosphere
    r_min = pow(((3*rho)/impactor_bulk_density), 1.0 / 3.0)*h; // if impactor is less than r_min no mass is ejected
    r_cap = pow(((3*sqrt(2*pi)*rho)/(4*impactor_bulk_density)), 1.0/3.0)*sqrt(h*R); // impactor radius needed in order to eject entire cap
    Mt = ((2*r)/r_min)*pow(1-pow((r_min/r),2),-1)*m1; // total impactor mass needed to eject the atmosphere as a function of impact radius r
    r_gi = pow((2*h*pow(R,2)), 1/3); //radius required for giant impact
    x = (v_imp/v_esc)*(m_imp/M1);
    

    if ((pi/2) - theta > sqrt(h/R))
    {
        m_eject = ((pow(sin(theta),2)*cos(theta))/2)*m_imp;
    }
        else if (theta < sqrt(h/R)) 
        {
            m_eject = (r_min/(2*r))*(1-pow((r_min/r),2))*m_imp;      
        }
    
    if (r>r_cap)
    {
        m_eject = M1cap; // maximum mass ejected 
    }
   
    if (r < r_min){ // if impactor radius is smaller than minimum radius required to eject any atmosphere
        m_eject =0;     
    }

    if (m_imp <= m_min){ // if mass of impactor is less than minimum mass required to eject any atmosphere
        m_eject = 0;    
    }
    if (m_imp >= m_cap){
        m_eject = M1cap;
    }
    
    if (answer == '1'){ // Adiabatic atmosphere 
        X_loss = .4*x + (1.8*pow(x,2)) -1.2*pow(x,3); // global mass loss for an adiabatic atmosphere
           
    }
    if (answer == '2'){ // Isothermal atmosphere
        X_loss = .4*x + (1.4*pow(x,2)) -0.8*pow(x,3); // global mass loss for an isothermal atmosphere 
    }

    if (m_eject > M1cap){
        m_eject = M1cap;
    }

    if (r_min<r<r_cap){ // N represents the number of impactors needed to eject the atmosphere
        N = 6*pow((R/r),2)*(rho*h)/(impactor_bulk_density*r_min);
    }

        else if (r_cap < r < r_gi){
            N = ((2*R)/h);  
        }



new_atmos = (m_eject/m1);   
std::cout<<"m_eject: "<<m_eject<<endl;
std::cout<<"m_eject/M1cap: "<<m_eject/M1cap<<endl;
std::cout<<"The fraction of the atmosphere that is ejected is: " << new_atmos<<" of the entire atmosphere"<<endl;
std::cout<<"r: "<<r<<endl;


        //sum = new_atmos;
        //sum +=i;
    }
    //sum = new_atmos;
    //std::cout<<sum<<endl;


}