
int init_dsandvels_katotorus(FTYPE r, FTYPE th, FTYPE phi, FTYPE *rhoout, FTYPE *uuout)
{
  ldouble rho,uint,Td,Ad,A0,Rd,Zd,Hd,Td0,Hfac; //Ad is semi-major axis
  
  Rd=r*sin(th); //cylindrical radius
  Zd=r*cos(th); //cylindrical height
  Hd=HR_DISK*Rd; 
  Ad=Rd*(1.+E_DISK*cos(phi))/(1.-E_DISK*E_DISK);
  Hfac=(1+E_DISK*cos(phi))/(1+E_DISK); // ratio H(phi=0)/H(phi)
  
  if(Ad<A_MIN || Ad>A_MAX) {*rhoout=-1.;return 0;}
  if(abs(Zd) >= 1.5*Hd) {*rhoout=-1.;return 0;}
  
  if(Ad<A_0)
  {
    rho = RHO0*exp(Ad-A_0)*exp(-Zd*Zd/2./Hd/Hd)*pow(Hfac,1.);
  }
  else
  {
    rho = RHO0*(A_0/Ad)*exp(-Zd*Zd/2./Hd/Hd)*pow(Hfac,1.);
  }
  Td = T_DISK0*(A_0/Ad);
  uint = uint=calc_PEQ_ufromTrho(Td,rho,0,0,0);

  //printf("r th phi const rho T0 Td uint = %e %f %f %e %e %e %e %e \n",r,th,phi,mugas_mp_over_kB,rho,T_DISK0,Td,uint);

  *rhoout = rho;
  *uuout = uint;

  return 0;

}

