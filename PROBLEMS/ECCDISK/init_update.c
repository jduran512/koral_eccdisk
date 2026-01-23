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

    //MYCOORDS vector potential
    ldouble Acov[4];
    Acov[0] = Acov[1] = Acov[2] = 0.;   //set t,r,th=0


    // =====================================================
    // NTORUS == 1 (classic SANE loop)
    // =====================================================
#if (NTORUS==1)
    Acov[3] = my_max(pp[RHO] - 0.5*RHO0, 0.);
    pp[B1] = Acov[1];
    pp[B2] = Acov[2];
    pp[B3] = Acov[3];
#endif



    // =====================================================
    // NTORUS == 2 OR NTORUS == 3
    // MAD geometry (Teixeira+2018), but NTORUS==3 uses SANE MAXBETA
    // =====================================================
#if (NTORUS==2 || NTORUS==3)
    {
        ldouble rBL  = geomBL.xx;
        ldouble thBL = geomBL.yy;
        ldouble phiBL = geomBL.zz;

        //eccentric semi-major axis a(r,phi)
        ldouble a = rBL * (1. + E_DISK*cos(phiBL)) / (1. - E_DISK*E_DISK);

        ldouble rterm=0., thterm=0.;


        // ------------------------------
        // Region 1: A_MIN < a < Atr2
        // ------------------------------
        if(a >= A_MIN && a < Atr2)
        {
            rterm = my_max(pow(a - A_MIN, V_B)*pow(10.,40.) - 0.02, 0.);

            //theta taper for both NTORUS==2 and NTORUS==3
            ldouble thlo = TH_CUT * M_PI;          //for NTORUS==3 this is new
            ldouble thhi = (1. - TH_CUT) * M_PI;

            if(thBL > thlo && thBL < thhi)
            {
                ldouble thmid = fabs(M_PI/2. - thBL);
                ldouble dth   = M_PI / fabs(thhi - thlo);
                thterm = pow(cos(thmid*dth), 1 + H_B);
            }
        }


        // ------------------------------
        // Region 2: Atr2 < a < Acut
        // ------------------------------
        else if(a >= Atr2 && a < Acut)
        {
            rterm = my_max(pow(Atr2 - A_MIN, V_B)*pow(10.,40.) - 0.02, 0.)
                  * cos(M_PI*(a - Atr2)/(Acut - Atr2)/2.);

            ldouble thlo = TH_CUT * M_PI;
            ldouble thhi = (1. - TH_CUT) * M_PI;

            if(thBL > thlo && thBL < thhi)
            {
                ldouble thmid = fabs(M_PI/2. - thBL);
                ldouble dth   = M_PI / fabs(thhi - thlo);
                thterm = pow(cos(thmid*dth), 1 + H_B*(Atr2/a));
            }
        }

        else
            rterm = 0.;   //outside region


        //A_phi
        Acov[3] = rterm * thterm;

        pp[B1] = Acov[1];
        pp[B2] = Acov[2];
        pp[B3] = Acov[3];
    }
#endif

#endif //MAGNFIELD



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
