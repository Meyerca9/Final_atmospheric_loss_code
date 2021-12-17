#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;
int main() {
    double m1,m2,m_imp,v_imp,rho,h,R,impactor_bulk_density,r,pi,M1,r_cap,M1cap,M2cap,
    new_atmos,cap_factor,mcap_ratio,m_min,m_max,r_min,m_eject,theta,M2,M_solid,M_atmos,
    x,X_loss,v_esc,Mt,mu,k,T,G,g,p_base,A,V,Rd,N,r_gi,sum,v_sesc,M_lr,M_slr,M_frag,v_sc,u,
    M_tot,Q_sc,Q_RD_prime,u_alpha,Q_RD_star,c,u_bar,qg,b,b_crit,gamma,Q,Q_ero,v_ero;


    impactor_bulk_density =  2000; // kg/m^3
    pi =  3.1415926535897;
    mu = 4.8e-23; // g mean molecular weight of atmosphere
    k = 1.38e-23; // J/K boltzmann constant
    G = 6.674e-11; // Gravitational constant
    Rd = 287.058; // Gas constant in J/(kg*K)
    c = 5; // dissipation of energy in small target bodies
    u_bar = .37; // material parameter

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
    

    std::cout <<"Enter the target mass";
    std::cin >> M1;

    std::cout <<"Enter the projectiles mass";
    std::cin >> m_imp;

    std::cout <<"Enter the Mass of the targets atmosphere";
    std::cin >> m1;

    std::cout <<"Enter the Mass of the projectiles atmosphere";
    std::cin >> m2;

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


    R = std::cbrt(M1/(1.33333*pi*5510)); // radius of target in km (5510 is impactor bulk density of earth) in m
    r = std::cbrt(m_imp/(1.33333*pi*impactor_bulk_density)); // radius of projectile in m
    g = (G*M1)/(pow(R,2)); // surface gravity 
    A = 4*pi*pow(R,2); // surface area of target m^2 ------
    p_base = (m1*g)/A; //pressure at bottom of an atmosphere in Pa
    rho = p_base/(Rd*T); // atmospheric density at the surface
    h = k*T/(mu*g)*1e3; // atmospheric scale height in km
    v_esc = sqrt((2*G*M1)/R); // escape velocity in km/s
    cap_factor = sqrt((pi*R)/(2*h)); // Impactor mass needed to eject all the atmospheric mass for planetesimals
    M1cap = ((2*pi)*(rho)*pow(h,2)*(R)); // maximum atmospheric mass that a single planetesimal impact can eject above the tangent plane
    mcap_ratio = sqrt((pi*h)/(8*R)); // impactor mass needed to eject all the atmospheric mass above tangent plane if impact veloity is comparable to escape velocity 
    m_min = 4*pi*rho*pow(h,3); // minimum mass needed in order to eject any atmosphere
    m_max = sqrt(2)*rho*pow((pi*h*R),3.0/2.0); // masses greater than this are able to eject all mass above the tangent line which is h/2r of the atmosphere
    r_min = pow(((3*rho)/impactor_bulk_density), 1.0 / 3.0)*h; // if impactor is less than r_min no mass is ejected
    r_cap = pow(((3*sqrt(2*pi)*rho)/(4*impactor_bulk_density)), 1.0/3.0)*sqrt(h*R); // impactor radius needed in order to eject entire cap
    Mt = ((2*r)/r_min)*pow(1-pow((r_min/r),2),-1)*m1; // total impactor mass needed to eject the atmosphere as a function of impact radius r
    r_gi = pow((2*h*pow(R,2)), 1/3); //radius required for giant impact
     
    gamma = (M1/m_imp); // target to impactor ratio
    qg = pow(1/8*(32*pi*c)/5,((3*u_bar)/2)); //coefficient for catastrophic distribution threshold
    x = (v_imp/v_esc)*(m_imp/M1);
    v_sesc = sqrt((2*G*m_imp)/(R+r)); // surface escape velocity of the interacting mass
    M_tot = M1+m_imp; // total mass
    u = (M1*m_imp)/M_tot; // reduced mass
    Q_RD_star = pow(qg*(impactor_bulk_density*G),(3*u_bar)/2)*pow((R+r),3*u_bar); // catastrophic disruption threshold
    Q_RD_prime = 1.9e2; // critical impact energy was assumed
    Q_sc = 1.8*Q_RD_prime; // Q_sc = super-catastrophic disruption energy
    v_sc = sqrt((2*Q_sc*M_tot)/u); // super-catastrophic critical velocity
    Q = .5*M1*pow(v_imp,2); // kinetic energy of the target mass
    Q_ero = Q_RD_prime*((2*m_imp)/M_tot); // critical erosive energy
    v_ero = sqrt((2*Q_ero*M_tot)/u);
    b = sin(theta); // impact parameter
    b_crit = (R/(R+r));


    // solid mass loss

    if (v_imp < v_sesc){ // perfect merge regime 
        M_lr = M1+m_imp;
        M_slr = 0;
        M_frag = 0;
    }
    if (v_imp >= v_sc){
        M_lr = .1*M_tot*pow((Q/Q_sc),-3/2);
    }
    if (v_imp >= v_ero){
        M_lr = M_tot*(1-(Q/(2*Q_RD_prime)));
    }
    if (b < b_crit && v_imp < v_ero){
        M_lr = M_tot*(1-(Q/(2*Q_RD_prime)));
    }

std::cout<<M_lr<<endl;
std::cout<<M_slr<<endl;
std::cout<<M_frag<<endl;
std::cout<<M_tot<<endl;



}
    

    
