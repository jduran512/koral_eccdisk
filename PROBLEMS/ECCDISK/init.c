int init_dsandvels_katotorus(FTYPE r, FTYPE th, FTYPE phi, FTYPE *rhoout, FTYPE *uuout);

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

#ifdef CONSISTENTGAMMA
ldouble gamma=calc_gammaintfromTei(1.e6,1.e6); //good faith estimate
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;
ldouble phi=geomBL.zz;

init_dsandvels_katotorus(r, th, phi, &rho, &uint); 
uintorg=uint;

if(rho<0.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    pgas = GAMMAM1 * uint;
    ell*=-1.;

    ldouble ult,ulr,ulph,ucov[4],ucon[4],drdt,dphidt;
    drdt=E_DISK*sin(phi)/sqrt(r*(1.+E_DISK*cos(phi)));
    dphidt=sqrt((1.+E_DISK*cos(phi))/r/r/r);
    ult = sqrt(-1./(geomBL.gg[0][0] + drdt*drdt*geomBL.gg[1][1] + 2.*dphidt*geomBL.gg[0][3] + dphidt*dphidt*geomBL.gg[3][3]));
    ulr = ult*drdt;
    ulph = ult*dphidt;

    ucon[0]=ult;
    ucon[1]=ulr;
    ucon[2]=0.;
    ucon[3]=ulph;
    
    //indices_12(ucov,ucon,geomBL.GG);
    indices_21(ucon,ucov,geomBL.gg);

    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
  
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2]; //derived in BL?
    pp[4]=ucon[3];

    //printf("ix iy iz uint uintbackup = %i %i %i %e %e \n",ix,iy,iz,uint,ppback[1]);

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

    


#ifdef RADIATION
    //distributing pressure
    ldouble P,aaa,bbb;
    P=GAMMAM1*uint;
    //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
    aaa=4.*SIGMA_RAD/3.;
    bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
    ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
    ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

    E=calc_LTE_EfromT(T4);
    Fx=Fy=Fz=0.;
    uint=calc_PEQ_ufromTrho(T4,rho,ix,iy,iz);

    pp[UU]=uint;//my_max(uint,ppback[1]);
    pp[EE0]=E;//my_max(E,ppback[EE0]);



    pp[FX0]=Fx;
    pp[FY0]=Fy;
    pp[FZ0]=Fz;

#ifdef EVOLVEPHOTONNUMBER
    pp[NF0]=calc_NFfromE(pp[EE0]);
#endif

    //transforming from BL lab radiative primitives to code non-ortonormal primitives
    prad_ff2lab(pp,pp,&geomBL);

#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
#ifdef MAGNFIELD 
    //MYCOORDS vector potential to calculate B's
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;

#if(NTORUS==1)
    //standard single poloidal loop
    Acov[3]=my_max(pp[RHO]-5.e-1*RHO0,0.);
#endif


    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];
#endif

   }

// for global field models instead of disk only
#ifdef MAGNFIELD 
// vector potential components (Acov[mu]) used to generate magnetic field
ldouble Acov[4];
Acov[0] = Acov[1] = Acov[2] = 0.;   // set t, r, theta components to zero

#if(NTORUS==2)
// B-field from Teixeira 2018, modified to follow eccentric disk orbits
ldouble rterm, thterm, thmin, thmax, thmid, dth;

// compute semi-major axis from (r, phi, eccentricity = E_DISK)
ldouble a = r * (1. + E_DISK * cos(phi)) / (1. - E_DISK*E_DISK);

// define polar angle region where loops are present
thmin = 0.2 * M_PI;                       // lower theta bound
thmax = (1. - 0.2) * M_PI;                // upper theta bound
thmid = fabs(M_PI/2. - th);               // distance from midplane
dth   = M_PI / fabs(thmax - thmin);       // normalization factor

// case 1: inside transition semi-major axis
if(a >= A_MIN && a < Atr2) 
{
  rterm = my_max(pow(a - A_MIN, V_B)*pow(10.,40.) - 0.02, 0.);   // radial profile
  if(th <= thmin) thterm = 0.;                                  // cut below theta min
  if(th >= thmax) thterm = 0.;                                  // cut above theta max
  if(th > thmin && th < thmax) thterm = pow(cos(thmid*dth), 1 + H_B); // angular taper
}
// case 2: between transition and cutoff semi-major axes
else if(a >= Atr2 && a < Acut)
{
  rterm = my_max(pow(Atr2 - A_MIN, V_B)*pow(10.,40.) - 0.02, 0.) // taper with cosine factor
        * cos(M_PI*(a - Atr2)/(Acut - Atr2)/2.);
  if(th <= thmin) thterm = 0.;
  if(th >= thmax) thterm = 0.;
  if(th > thmin && th < thmax) thterm = pow(cos(thmid*dth), 1 + H_B*(Atr2/a)); // theta taper
}
// case 3: outside cutoffs means no magnetic field
else if(a < A_MIN || a > Acut)
{
  rterm = 0.;
}

// assign phi component of vector potential
Acov[3] = rterm * thterm;

// compute magnetic field components from vector potential
pp[B1] = Acov[1];   // Br
pp[B2] = Acov[2];   // Btheta
pp[B3] = Acov[3];   // Bphi
#endif
#endif


//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],ix,iy,iz);

//thermodynamics and relativistic electrons
#ifdef EVOLVEELECTRONS
#ifdef  RELELECTRONS
int ie;
for (ie=0; ie < NRELBIN; ie++)
{
  if(relel_gammas[ie]<RELEL_INJ_MIN || relel_gammas[ie]>RELEL_INJ_MAX)
  pp[NEREL(ie)] = 0.0 ; //Always initialize with zero nonthermal
  else
  pp[NEREL(ie)] = pow(relel_gammas[ie],-RELEL_HEAT_INDEX);
}

ldouble unth_tmp = calc_relel_uint(pp);
if (unth_tmp > 0)
{
    ldouble scalefac = 1.e-5*pp[UU]/unth_tmp;
    for(ie=0;ie<NRELBIN;ie++)
    {
      pp[NEREL(ie)] *= scalefac;
    }
}
#endif //RELELECtRONS
ldouble rhogas=pp[RHO];
ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);

//slightly colder electrons initially: Should be same T for super-Critical problem?
ldouble ue=pp[UU];//1./100.*pp[UU];
ldouble ui=pp[UU];//(1.-1./100.)*pp[UU];
ldouble Te,Ti;
pp[ENTRE]=calc_Sefromrhou(calc_thermal_ne(pp)*MU_E*M_PROTON,ue,ELECTRONS);
pp[ENTRI]=calc_Sefromrhou(rhogas,ui,IONS);

#ifdef CONSISTENTGAMMA
Ti=solve_Teifromnmu(pp[RHO]/MU_I/M_PROTON, M_PROTON, ui, IONS); //solves in parallel for gamma and temperature
Te=solve_Teifromnmu(pp[RHO]/MU_E/M_PROTON, M_ELECTR, ue, ELECTRONS); //solves in parallel for gamma and temperature
    
#ifdef RELELECTRONS
gamma = calc_gammaint_relel(pp, Te, Ti);
#else
gamma=calc_gammaintfromTei(Te,Ti); //good faith estimate
#endif 

set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif //CONSISTENTGAMMA
pp[ENTRE]=calc_SefromrhoT(rhogas,Te,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(rhogas,Ti,IONS);
#endif //EVOLVEELECTRONS

//to conserved
p2u(pp,uu,&geom);



/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
