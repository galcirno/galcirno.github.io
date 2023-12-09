! ADGZPRESCOLD.  PRESCOLD, CALCULATION OF THE HALF-LIFE FOR ALPHA DECAY,  ADGZ0000
! 1   CLUSTER RADIOACTIVITY AND COLD FISSION PROCESSES.  M. GONCALVES,    ADGZ0000
! 2   S.B. DUARTE, F. GARCIA, O. RODRIGUEZ.                               ADGZ0000
! REF. IN COMP. PHYS. COMMUN. 107 (1997) 246                              ADGZ0000
! CPC FREE FORMAT DATA 4468 CARDS                                                 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! %%%                     PRESCOLD.FOR  --- next 918 lines                    %%% 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
C       %    This program calculates the half lives, barrier  penetrability   % 
C       %    and yields for the alpha decay, cluster radioactivity and cold   % 
C       %    fission process, by taking into account two different inertial   % 
C       %    coefficients according to the shape parametrization choosen to   % 
C       %    to describe the dynamical evolution of the nuclear system, ie,   % 
C       %    cluster like or fission like shape.                              % 
C       %                                                                     % 
C       %                    Present Version: 01/08/97                        % 
C       %                                                                     % 
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                                                                
      Implicit None                                                             
      Real*8 Ak, Ak1, Ak2, Deltm, Deltm1, Deltm2, E, Excess,                    
     1       Excess_daugh, Excess_emit, Excess_parent, Mas, Mas1, Mas2,         
     2       Uma, Zk, Zk1, Zk2                                                  
      Real*8 Aux1, Aux2                                                         
      Integer*4 I, Idaughter, Iemitted, Key_like, L, Mass_kind,                 
     1          N, NAatomic, NAdaughter, NAemitted, NAparent, NAtemp,           
     2          Number_of_Lines, NZatomic, NZdaughter, NZemitted,               
     3          NZParent, NZtemp                                                
      Integer*4 Ieb, Idb                                                        
                                                                                
      Dimension NAatomic(3000), NZatomic(3000), Excess(3000)                    
                                                                                
c------------------------------------------------------------------------------ 
      Open(1,File='Mass_exc.dat')                                               
      Open(2,File='Out.dat')                                                    
      Open(3,File='In.dat')                                                     
c------------------------------------------------------------------------------ 
c     Key_like is the key to choose between the Cluster like (1) or             
c     Fission like (2) shapes. Mass_kind is the key to choose between           
c     The Werner-Wheeler (1) of the Effective (2) inertia coefficients.         
      Read(3,*) NAparent, NZparent, Key_like, Mass_kind, L                      
      Uma = 931.501d0     !    The nuclear mass Unit                            
      N = 100             !    Number of steps for numerical integration        
      Number_of_Lines = 2930  !    Number of lines in the input file            
      Aux1 = 1.d100                                                             
      Aux2 = 40.0                                                               
c------------------------------------------------------------------------------ 
      Do i = 1, Number_of_Lines                      !   Loop for reading the   
         Read(1,*)NAatomic(i),NZatomic(i),Excess(i)  !   input data file.       
         If ( (NAatomic(i) .eq. NAparent) .and.                                 
     1        (NZatomic(i) .eq. NZparent) ) Excess_parent = Excess(i)           
      Enddo                                                                     
c------------------------------------------------------------------------------ 
      Ieb = 1                                                                   
c-----------------------------------------------------------------------------  
      Do Iemitted = Ieb, Number_of_Lines                                        
         Write(*,*)'Iemitted:',Iemitted                                         
         NAemitted = NAatomic(Iemitted)                                         
         NZemitted = NZatomic(Iemitted)                                         
         Excess_emit = Excess(Iemitted)                                         
         Idb = Iemitted                                                         
         Do Idaughter = Idb, Number_of_Lines                                    
            NAdaughter = NAatomic(Idaughter)                                    
            NZdaughter = NZatomic(Idaughter)                                    
            Excess_daugh = Excess(Idaughter)                                    
            NAtemp = NAemitted + NAdaughter                                     
            NZtemp = NZemitted + NZdaughter                                     
            If ((NAtemp.eq.NAparent).and.(NZtemp.eq.NZparent))Then              
               Mas  = Uma*NAparent  + Excess_parent                             
               Mas1 = Uma*NAemitted + Excess_emit                               
               Mas2 = Uma*NAdaughter + Excess_daugh                             
               Deltm = Excess_parent                                            
               Deltm1 = Excess_emit                                             
               Deltm2 = Excess_daugh                                            
               Ak = NAparent                                                    
               Ak1 = NAemitted                                                  
               Ak2 = NAdaughter                                                 
               Zk = NZparent                                                    
               ZK1 = NZemitted                                                  
               ZK2 = NZdaughter                                                 
               E =Deltm-(Deltm1+Deltm2)                                         
               If (E .Gt. 0.d0) Then                                            
                  If (Key_like .eq. 1) Then                                     
                     Call Penetrability_Cluster_Like(Mas1, Mas2,                
     1               Zk, Zk1, Zk2, Ak, Ak1, Ak2, E, N, L, Mass_kind)            
                  Endif                                                         
                  If (Key_like .eq. 2) Then                                     
                     Call Penetrability_Fission_Like(Mas, Mas1, Mas2,           
     1               Zk, Zk1, Zk2, Ak, Ak1, Ak2, E, N, L, Mass_kind)            
                  Endif                                                         
               Else                                                             
                  Write(2,5)Ak, Zk, Ak1, Zk1, Aux1, Aux2   !  Decay Forbiden    
               Endif                                                            
            Endif       ! End of Iparent If                                     
         Enddo       ! End of Idaugther loop                                    
      Enddo       ! End of Iemitted loop                                        
                                                                                
 5    Format(4F5.0, E12.4, F12.4)                                               
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Gamow penetrability factor calculation with the effective potential.     C
C     This potential is obtained as a sum of the Coulomb potential and the     C
C     surface potential for two spherical segments in separation. The          C
C     Coulomb part is taken from reference [Gaudin, J. de Physique ATTENTION]. C
C     The surface part is obtained by considering that the two fission         C
C     products have null total excitation energy and there is no neutron       C
C     emission. Also the frequency of assaults on the potential barrier is     C
C     set constant for all parent nuclei.                                      C
C     Here we are considering the Cluster-like shape to describe the dynamical C
C     evolution of the dimolecular phase.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Penetrability_Cluster_Like(Mas1, Mas2, Zk, Zk1,                
     +           Zk2, Ak, Ak1, Ak2, E, N, L, Mass_kind)                         
      Implicit None                                                             
      Real*8 Ak, Ak1, Ak2, Alpha, Bconst1, Bconst2, Cc12, Ce,                   
     1       Coeff_reduced, Coul0, Coul01, Coul02, Cs, C1, Dz, E, Ec012,        
     2       Es01, Es012, Es02, EmV, EmVm, EmV0, EmVmax, Eta, Factor,           
     3       Gamow, Hbc, Lambda, Lambda0, Log_of_Tau, Log_of_2, Mas1,           
     4       Mas2, Mu, Mu_reduced, One_third, Penetr, Rbar, Rhomas1,            
     5       Rhomas2, Rn0, Rp, Rtouch, Pi, R01, R02, R2, Sigma, Simpar,         
     6       Sn, Sn1, Sn2, Spar, Tau, V, Volp0, Volp02, Vol1, Vol10,            
     7       Vol2, Vol20, Z, Zinner, Zk, Zk1, Zk2                               
      Integer*4 K, L, Loop_Exceed, Mass_kind, My_error, N                       
                                                                                
      Data Pi/3.141592653589793d0/                                              
c-------------------------------------------------------------------------------
      C1 = 8.d0/9.d0*Pi     !   Setting some constants useful to the calculation
      Ce = 1.4399784d0                                                          
      Hbc= 197.327d0  !197.47d0                                                             
      One_third = 1.d0/3.d0                                                     
c-------------------------------------------------------------------------------
      If (Mass_kind .eq. 1) Rn0 = 1.37d0            !   Nuclear radius parameter
      If (Mass_kind .eq. 2) Rn0 = 1.20d0            !   Nuclear radius parameter
c-------------------------------------------------------------------------------
      Rp = Rn0*Ak**One_third           !    Parent nucleus radius calculation   
      Volp0 = 4.d0*One_third*Pi*Rp**3                                           
      Vol10 = Zk1*Volp0/Zk                                                      
      Vol20 = Zk2*Volp0/Zk                                                      
      R01 = (0.75d0*Vol10/Pi)**One_third !   Emitted nucleus radius calculation 
      R02 = (0.75d0*Vol20/Pi)**One_third !  Daughter nucleus radius calculation 
c------------------------------------------------------------------------------ 
      Rbar = Rp - R01                                  !   Inner  turning point 
      Rtouch = R01 + R02                               !   Outter turning point 
c------------------------------------------------------------------------------ 
      Mu_reduced = Mas1*Mas2 / (Mas1 + Mas2)     ! Reduced mass and coefficient 
      Coeff_reduced = Dsqrt(2.d0*Mu_reduced)/Hbc ! for the asymptotic region.   
      Rhomas1 = Mas1/Vol10      !   Asymptotic mass density of the two emitted  
      Rhomas2 = Mas2/Vol20      !   fragments. They will be useful for the      
      Bconst1 = Rhomas1         !   matching in the Touching Point.             
      Bconst2 = Rhomas2                                                         
c------------------------------------------------------------------------------ 
      Volp02 = Volp0*Volp0      !   Eletrical constants useful to the Coulomb   
      Factor = Ce*Zk*Zk/Volp02  !   energy calculation.                         
      Cc12 = Zk1*Zk2*Ce                                                         
c-----------------------------------------------------------------------------  
      Coul0  = 0.6d0 * Ce * Zk**2/Rp     !   Calculation of the Coulomb and     
      Coul01 = 0.6d0 * Ce * Zk1**2/R01   !   Surface self-energies of the       
      Coul02 = 0.6d0 * Ce * Zk2**2/R02   !   parent nuclei and of the final     
      Ec012 = Coul01 + Coul02            !   fragments.                         
      Sigma = (E + Ec012 - Coul0)/(4.d0*Pi*(Rp**2 - R01**2 - R02**2))           
      Cs = 2.d0 * Pi * Sigma                                                    
      Es01 = 2.d0 * Cs * R01*R01                                                
      Es02 = 2.d0 * Cs * R02*R02                                                
      Es012= Es01 + Es02                                                        
c------------------------------------------------------------------------------ 
      Zinner = Rbar                                                             
      Dz = (Rtouch - Zinner)/(2.*N)   !    Setting some variable to start the   
      Z = Zinner                      !    integration loop for the prescission 
      Spar = 0.d0                     !    region. We will integrate by using   
      Simpar = 0.d0                   !    the Simpsom Method.                  
      EmV0 = 0.d0                                                               
      EmVmax = 0.d0                                                             
c------------------------------------------------------------------------------ 
      Mu = 0.d0                                                                 
c------------------------------------------------------------------------------ 
      Do K = 1, N-1           ! Starting the loop for Simpson's integration     
c------------------------------------------------------------------------------ 
         Z = Z + Dz            ! Beggining of the odd part of the Simpsom's Sum 
         Call Pnos_Cl(Rp, R01, Ec012, Es012, Cs, C1, Factor,                    
     +             Z, R2, V, Vol1, Vol2, Alpha, Loop_Exceed)                    
         If (Loop_Exceed .Ge. 10000) Then                                       
            Log_of_Tau = -100.d0                                                
            Return                                                              
         Endif                                                                  
         EMVm = V-E                                                             
         If (EMVm .gt. EMVmax) EMVmax = EMVm                                    
         If (Mass_kind .eq. 1) Call MuWWcl (Bconst1, Bconst2, R01,              
     +      R2, Z, Vol1, Vol2, Mu)                                              
         If (Mass_kind .eq. 2) Call MuEffcl(E, Mas1, Rhomas1,                   
     +      Rhomas2, Volp0, Vol10, Vol1, Vol2, Alpha, Mu)                       
         EmV = Dsqrt(Mu*Dabs(V-E))                                              
         Simpar = Simpar + EmV       ! End of the odd part of the Simpsom's Sum 
c------------------------------------------------------------------------------ 
         Z = Z + Dz           ! Beggining of the even part of the Simpsom's Sum 
         Call Pnos_Cl(Rp, R01, Ec012, Es012, Cs, C1, Factor,                    
     +             Z, R2, V, Vol1, Vol2, Alpha, Loop_Exceed)                    
         If (Loop_Exceed .Ge. 10000) Then                                       
            Log_of_Tau = -100.d0                                                
            Return                                                              
         Endif                                                                  
         EmVm = V - E                                                           
         If (EmVm .gt. EmVmax) EmVmax = EmVm                                    
         If (Mass_kind .eq. 1) Call MuWWcl (Bconst1, Bconst2, R01,              
     +      R2, Z, Vol1, Vol2, Mu)                                              
         If (Mass_kind .eq. 2) Call MuEffcl(E, Mas1, Rhomas1,                   
     +      Rhomas2, Volp0, Vol10, Vol1, Vol2, Alpha, Mu)                       
         EmV = Dsqrt(Mu*Dabs(V-E))                                              
         Spar = Spar + EmV          ! End of the even part of the Simpsom's Sum 
      Enddo       ! End of the loop for integration in the post scission region 
c------------------------------------------------------------------------------ 
      Z = Z + Dz                                             !  Calculating the 
      Call Pnos_Cl(Rp, R01, Ec012, Es012, Cs, C1, Factor,       !  potential in the                                                                             
     +          Z, R2, V, Vol1, Vol2, Alpha, Loop_Exceed)    !  touching point  
      If (Loop_Exceed .Ge. 10000) Then                       !  to complete the 
         Log_of_Tau = -100.d0                                !  the integral    
         Return                                              !  calculation.    
      Endif                                                                     
      EmVm = V - E                                                              
      If (EmVm .gt. EmVmax) EmVmax = EmVm                                       
      If (Mass_kind .eq. 1) Call MuWWcl (Bconst1, Bconst2, R01, R2,             
     +   Z, Vol1, Vol2, Mu)                                                     
      If (Mass_kind .eq. 2) Call MuEffcl(E, Mas1, Rhomas1, Rhomas2,             
     +   Volp0, Vol10, Vol1, Vol2, Alpha, Mu)                                   
      EmV = Dsqrt(Mu*Dabs(V-E))                                                 
      Simpar = Simpar + EmV         ! End of the odd part of the Simpsom's Sum  
                                                                                
      Z = Z + Dz                                                                
      V = Cc12/Z                                                                
      EmV = Coeff_reduced*Dsqrt(Dabs(V - E)) ! Sn1 below is the total integral  
      Spar = Spar + EmV                      ! for the prescission region.      
                                                                                
      Sn1 = 1.41421356d0 * (EmV0 +4.d0*Simpar + 2.d0*Spar-EmV)*                 
     +                        Dz/(3.d0*Hbc)                                     
c------------------------------------------------------------------------------ 
c--------  Integral calculation for the Post Sscission region    -------------- 
      CALL Only_Coulomb(Cc12, E, Rtouch, Coeff_reduced,                         
     +                  Mu_reduced, L, Sn2, My_Error)                           
c------------------------------------------------------------------------------ 
      If (My_Error .Eq. 1) Then                                                 
         Log_of_Tau = -100.d0                                                   
         Return                                                                 
      Endif                                                                     
c------------------------------------------------------------------------------ 
      Sn = Sn1 + Sn2                            !    Total Gamow integral value 
c------------------------------------------------------------------------------ 
      Lambda0 = 1.d22                  !   Frequency of assaults on the barrier 
      Gamow = -2.d0*Sn                                                          
      If (Gamow .lt. -400.d0) Gamow = -400.d0                                   
      Penetr = Dexp(Gamow)                                                      
      Lambda = Lambda0 * Penetr                                                 
      Log_of_2 = 0.30102999d0                                                   
      Tau = Log_of_2 / Lambda                                                   
      Log_of_Tau = Dlog10(Tau)                                                  
      Eta = Abs(Ak2-Ak1)/Ak         !   Eta is the mass asymetry of the process 
      If (Log_of_Tau .ge. -40.d0 .and. Log_of_Tau .le. 40.d0) Then              
         Write(2,5)Ak, Zk, Ak1, Zk1, Penetr, Log_of_Tau                         
      Else                                                                      
         Write(2,5)Ak, Zk, Ak1, Zk1, 1.d100, 40.d0   !  Decay Forbiden          
      Endif                                                                     
c------------------------------------------------------------------------------ 
 5    Format(4F5.0, E12.4, F12.4)                                               
c------------------------------------------------------------------------------ 
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Gamow penetrability factor calculation with the effective potential.     C
C     This potential is obtained as a sum of the Coulomb potential and the     C
C     surface potential for two spherical segments in separation. The          C
C     Coulomb part is taken from reference [Gaudin, J. de Physique ATTENTION]. C
C     The surface part is obtained by considering that the two fission         C
C     products have null total excitation energy and there is no neutron       C
C     emission. Also the frequency of assaults on the potential barrier is     C
C     set constant for all parent nuclei.                                      C
C     Here we are considering the Fission-like shape to describe the dynamical C
C     evolution of the dimolecular phase.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine Penetrability_Fission_Like(Mas, Mas1, Mas2, Zk, Zk1,           
     +           Zk2, Ak, Ak1, Ak2, E, N, L, Mass_kind)                         
      Implicit None                                                             
      Real*8 Ak, Ak1, Ak2, Bconst1, Bconst2, Bt, Cc12, Ce,                      
     1       Coeff_reduced, Coul0, Coul01, Coul02, Cs, C1, Dz, E, Ec012,        
     2       Es01, Es012, Es02, EmV, EmVm, EmV0, EmVmax, Eta, Factor,           
     3       Gamow, Hbc, Lambda, Lambda0, Log_of_Tau, Log_of_2, Mas,            
     4       Mas1, Mas2, Mu, Mu_reduced, One_third, Penetr, Rbar,               
     5       Rhomas1, Rhomas2, Rn0, Rp, RT, Rtouch, R1t, R2t, Pi, R01,          
     6       R02, R1, R2, Sigma, Sum, Sn, Sn1, Sn2, Tau, V,                     
     7       Volp0, Volp02, Vol10, Vol20, Z, Zinner, Zk,                        
     8       Zk1, Zk2, Z0                                                       
      Integer*4 Dimen, K, L, Mass_kind, My_error, N                             
      Parameter (Dimen=200)                                                     
                                                                                
      Dimension Bt(Dimen), Rt(Dimen), R1t(Dimen), R2t(Dimen)                    
                                                                                
      Data Pi/3.141592653589793d0/                                              
c-------------------------------------------------------------------------------
      N = Dimen                                                                 
c-------------------------------------------------------------------------------
      C1 = 8.d0/9.d0*Pi     !   Setting some constants useful to the calculation
      Ce = 1.4399784d0                                                          
      Hbc= 197.47d0                                                             
      One_third = 1.d0/3.d0                                                     
c-------------------------------------------------------------------------------
      If (Mass_kind .eq. 1) Rn0 = 1.31d0            !   Nuclear radius parameter
      If (Mass_kind .eq. 2) Rn0 = 1.17d0            !   Nuclear radius parameter
c-------------------------------------------------------------------------------
      Rp = Rn0*Ak**One_third           !    Parent nucleus radius calculation   
      Volp0 = 4.d0*One_third*Pi*Rp**3                                           
      Vol10 = Zk1*Volp0/Zk                                                      
      Vol20 = Zk2*Volp0/Zk                                                      
      R01 = (0.75d0*Vol10/Pi)**One_third !   Emitted nucleus radius calculation 
      R02 = (0.75d0*Vol20/Pi)**One_third !  Daughter nucleus radius calculation 
c------------------------------------------------------------------------------ 
      Rbar = Rp - R01                                  !   Inner  turning point 
      Rtouch = R01 + R02                               !   Outter turning point 
c------------------------------------------------------------------------------ 
      Mu_reduced = Mas1*Mas2 / (Mas1 + Mas2)     ! Reduced mass and coefficient 
      Coeff_reduced = Dsqrt(2.d0*Mu_reduced)/Hbc ! for the asymptotic region.   
      Rhomas1 = Mas1/Vol10      !   Asymptotic mass density of the two emitted  
      Rhomas2 = Mas2/Vol20      !   fragments. They will be useful for the      
      Bconst1 = Rhomas1         !   matching in the Touching Point.             
      Bconst2 = Rhomas2                                                         
c------------------------------------------------------------------------------ 
      Volp02 = Volp0*Volp0      !   Eletrical constants useful to the Coulomb   
      Factor = Ce*Zk*Zk/Volp02  !   energy calculation.                         
      Cc12 = Zk1*Zk2*Ce                                                         
c-----------------------------------------------------------------------------  
      Coul0  = 0.6d0 * Ce * Zk**2/Rp     !   Calculation of the Coulomb and     
      Coul01 = 0.6d0 * Ce * Zk1**2/R01   !   Surface self-energies of the       
      Coul02 = 0.6d0 * Ce * Zk2**2/R02   !   parent nuclei and of the final     
      Ec012 = Coul01 + Coul02            !   fragments.                         
      Sigma = (E + Ec012 - Coul0)/(4.d0*Pi*(Rp**2 - R01**2 - R02**2))           
      Cs = 2.d0 * Pi * Sigma                                                    
      Es01 = 2.d0 * Cs * R01*R01                                                
      Es02 = 2.d0 * Cs * R02*R02                                                
      Es012= Es01 + Es02                                                        
c------------------------------------------------------------------------------ 
      Zinner = 1.d-8              !     Setting some variable to start the      
      Z = Zinner                  ! integration loop for the prescission region 
      EmV0 = 0.d0                                                               
      EmVmax = 0.d0                                                             
c------------------------------------------------------------------------------ 
      Mu = 0.d0                                                                 
      Sum = 0.d0                                                                
      Call Configuration (N, Mass_kind, Rp, R01, R02, Bconst1, Bconst2,         
     +     Mas, Mas1, Mu_reduced, Ak, Ak2, Rt, R1t, R2t, Bt)                    
c------------------------------------------------------------------------------ 
      Do K = 1, N     ! Starting the loop for integration in prescission region 
c------------------------------------------------------------------------------ 
         Z0 = Z                                                                 
         Z = Rt(k)                                                              
         dZ = Dabs(Z - Z0)                                                      
         R1 = R2t(k)                                                            
         R2 = R1t(k)                                                            
         Call Pnos_Fl(Ec012, Es012, Cs, C1, Factor, Z, R1, R2, V)               
         EMVm = V-E                                                             
         IF (EMVm .gt. EMVmax) EMVmax = EMVm                                    
         If (Mass_kind .eq. 1) Call MuWWfl (N, Rt, Bt, Z, Mu)                   
         If (Mass_kind .eq. 2) Call MuEfffl(N, Rt, Bt, Z, Mu)                   
         EmV = Dsqrt(Mu*Dabs(V-E))                                              
         Sum = Sum + EmV*dZ                                                     
c------------------------------------------------------------------------------ 
      Enddo       !   End of the loop for integration in the prescission region 
c------------------------------------------------------------------------------ 
      Sn1 = Dsqrt(2.d0) * Sum / Hbc                                             
c------------------------------------------------------------------------------ 
c--------  Integral calculation for the Post Scission region    --------------  
      CALL Only_Coulomb(Cc12, E, Rtouch, Coeff_reduced,                         
     +                  Mu_reduced, L, Sn2, My_Error)                           
c------------------------------------------------------------------------------ 
      If (My_Error .Eq. 1) Then                                                 
         Log_of_Tau = -100.d0                                                   
         Return                                                                 
      Endif                                                                     
c------------------------------------------------------------------------------ 
      Sn = Sn1 + Sn2                            !    Total Gamow integral value 
c------------------------------------------------------------------------------ 
      Lambda0 = 1.d22                  !   Frequency of assaults on the barrier 
      Gamow = -2.d0*Sn                                                          
      If (Gamow .lt. -400.d0) Gamow = -400.d0                                   
      Penetr = Dexp(Gamow)                                                      
      Lambda = Lambda0 * Penetr                                                 
      Log_of_2 = Dlog10(2.d0)                                                   
      Tau = Log_of_2 / Lambda                                                   
      Log_of_Tau = Dlog10(Tau)                                                  
      Eta = Abs(Ak2-Ak1)/Ak         !   Eta is the mass asymetry of the process 
      If (Log_of_Tau .ge. -40.d0 .and. Log_of_Tau .le. 40.d0) Then              
         Write(2,5)Ak, Zk, Ak1, Zk1, Penetr, Log_of_Tau                         
      Else                                                                      
         Write(2,5)Ak, Zk, Ak1, Zk1, 1.d100, 40.d0   !  Decay Forbiden          
      Endif                                                                     
c------------------------------------------------------------------------------ 
 5    Format(4F5.0, E12.4, F12.4)                                               
c------------------------------------------------------------------------------ 
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C.......CALCULO DE NOSSO POTENCIAL PARA REGIAO DE PRE-CISAO.............        
c---- This subroutine enable to calculate the total potencial of our work for   
c---- two spherical fragments in separation. Z is the variable distance between 
c---- the geometrical centers of the two spheres. The present potential consists
c---- in a sum of the Coulomb energy and the surface energy of the fragments    
c---- during the overlaping, and only the Coulomb energy for the region post-   
c---- scission. In this subroutine Alpha comes out as the useful factor to      
c---- obtain the effective inertial coefficient.                                
      Subroutine Pnos_Cl(Rp, R01, Ec012, Es012, Cs, C1, Factor,                 
     +                Z, R2, V, Vol1, Vol2, Alpha,Loop_Exceed)                  
      Implicit None                                                             
      Real*8 Aa, Alpha, Cs, C1, Ec, Ecn, Ec012, Es, Esn, Es012, Fat,            
     +       Factor, Pi, Rp, R01, R2, Tau1, Tau2, Tet1, Tet2, V, Vol1,          
     +       Vol2, Z, Z1, Z1bar, Z2, Z2bar                                      
      Integer*4 Loop_Exceed                                                     
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      Call NewtZ(Rp, R01, Z, R2, Loop_Exceed)  !   Numerical solution of the third                                                                           
      Z2 = (R2**2 - R01**2 + Z**2)/(2.d0*Z)    !   degree equation for R2 as a  
      Z1 = Z - Z2                              !   function of Z.               
      Aa = Dsqrt(R2**2 - Z2**2)                !   Aa is the radius of the circle of                                                                       
      Tet1 = Dacos(Z1/R01)                     !   base.                        
      Tet2 = 2.d0*Pi - Dasin(Aa/R2)                                             
      Tau1 = Tet1 - Pi                         !   Tau1 and Tau2 are the angles of                                                                              
      Tau2 = Tet2 - Pi                         !   the Gaudin's paper           
      Call Ecoul(C1, Factor, Tau1, Tau2, Aa, Ec)                                
      Call Sup(Cs, Tet1, Tet2, Aa, Es)         !   Ec and Es are the configuration                                                                   
      Ecn = Ec - Ec012                         !   dependent Coulomb and Surface
      Esn = Es - Es012                         !   energies.                    
      V = Ecn + Esn                                                             
c----------------------------------------------!   Vol1 and Vol2 are the configuration                                                                   
      Vol1 = Pi*(R01+Z1)**2*(2.d0*R01-Z1)/3.d0 !   dependent volumes of the two sphe-                                                                           
      Vol2 = Pi*(R2 +Z2)**2*(2.d0*R2 -Z2)/3.d0 !   rical segments.              
c----------------------------------------------!   Alpha calculation            
      Z1bar = Pi*(R01**2 - Z1**2)**2/(4.d0*Vol1)                                
      Z2bar = Pi*(R2 **2 - Z2**2)**2/(4.d0*Vol2)                                
      Fat = 2.d0/(Z*(R2 - Z2))                                                  
      Alpha = 1.d0 - Fat*(Z1*(Z1bar + Z2bar) + Z1bar**2 - Z2bar**2)             
c------------------------------------------------------------------------------ 
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C.......CALCULO DE NOSSO POTENCIAL PARA REGIAO DE PRE-CISAO.............        
c---- This subroutine enable to calculate the total potencial of our work for   
c---- two spherical fragments in separation. Z is the variable distance between 
c---- the geometrical centers of the two spheres. The present potential consists
c---- in a sum of the Coulomb energy and the surface energy of the fragments    
c---- during the overlaping, and only the Coulomb energy for the region post-   
c---- scission. In this subroutine Alpha comes out as the useful factor to      
c---- obtain the effective inertial coefficient.                                
      Subroutine Pnos_Fl(Ec012, Es012, Cs, C1, Factor, Z, R1, R2, V)            
      Implicit None                                                             
      Real*8 Aa, Alpha, Cs, C1, Ec, Ecn, Ec012, Es, Esn, Es012, Fat,            
     +       Factor, Pi, R1, R2, Tau1, Tau2, Tet1, Tet2, V, Vol1,               
     +       Vol2, Z, Z1, Z1bar, Z2, Z2bar                                      
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      Z2 = (R2**2 - R1**2 + Z**2)/(2.d0*Z)                                      
      Z1 = Z - Z2                                                               
      Aa = Dsqrt(R2**2 - Z2**2)                !   Aa is the radius of the circle of                                                                       
      Tet1 = Dacos(Z1/R1)                      !   base.                        
      Tet2 = 2.d0*Pi - Dasin(Aa/R2)                                             
      Tau1 = Tet1 - Pi                         !   Tau1 and Tau2 are the angles of                                                                              
      Tau2 = Tet2 - Pi                         !   the Gaudin's paper           
      Call Ecoul(C1, Factor, Tau1, Tau2, Aa, Ec)                                
      Call Sup(Cs, Tet1, Tet2, Aa, Es)         !   Ec and Es are the configuration                                                                   
      Ecn = Ec - Ec012                         !   dependent Coulomb and Surface
      Esn = Es - Es012                         !   energies.                    
      V = Ecn + Esn                                                             
c----------------------------------------------!   Vol1 and Vol2 are the configuration                                                                   
      Vol1 = Pi*(R1+Z1)**2*(2.d0*R1-Z1)/3.d0 !   dependent volumes of the two sphe-                                                                           
      Vol2 = Pi*(R2+Z2)**2*(2.d0*R2-Z2)/3.d0 !   rical segments.                
c----------------------------------------------!   Alpha calculation            
      Z1bar = Pi*(R1**2 - Z1**2)**2/(4.d0*Vol1)                                 
      Z2bar = Pi*(R2**2 - Z2**2)**2/(4.d0*Vol2)                                 
      Fat = 2.d0/(Z*(R2 - Z2))                                                  
      Alpha = 1.d0 - Fat*(Z1*(Z1bar + Z2bar) + Z1bar**2 - Z2bar**2)             
c------------------------------------------------------------------------------ 
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Calculation of the Gamow Exponential for the post-scission region by      
c---- using the de Carvalho analitical expression. (Internal communication)     
      Subroutine Only_Coulomb(Cc12, E, Rtouch, CoeffWW, MuWW, L,                
     +                        Sn, My_Error)                                     
      Implicit None                                                             
      Real*8 Aux, C, Cc12, CoeffWW, CoverD, CoverD2, D, DoverC, E, EL,          
     +       ELoverE, G, MuWW, Rc, Raiz1, Raiz2, Rtouch, Sn,                    
     +       Term1, Term2, Term3                                                
      Integer*4 L, My_Error                                                     
                                                                                
      Rc = Cc12/E                                                               
      C = Rtouch                                                                
      D = Rc                                                                    
      My_Error = 0                                                              
      If (C .Ge. D) Then                                                        
         My_Error = 1                                                           
         Return                                                                 
      Endif                                                                     
      CoverD = C/D                                                              
      DoverC = D/C                                                              
      Call Centrif(MuWW, C, L, EL)                                              
      ELoverE = EL/E                                                            
      CoverD2 = CoverD*CoverD                                                   
      Aux = ELoverE*(ELoverE + DoverC - 1.d0)                                   
      Raiz1 = Dsqrt( Aux )                                                      
      Aux = 1.d0 + 4.d0*CoverD2*ELoverE                                         
      Raiz2 = Dsqrt( Aux )                                                      
                                                                                
      Term1 = CoverD*Dsqrt(ELoverE)*Dlog((Raiz1 + ELoverE +0.5d0*DoverC)        
     .        /(2.d0*CoverD*ELoverE/(1.d0 + Raiz2) + 0.5d0*DoverC) )            
      Aux = 0.5d0*(1.d0 - (1.d0 - 2.d0*CoverD)/Raiz2)                           
      Term2 = Dacos(Dsqrt(Aux))                                                 
      Aux = CoverD2*ELoverE + CoverD - CoverD2                                  
      Term3 = Dsqrt( Aux )                                                      
      G = Term1 + Term2 - Term3                                                 
      Sn= CoeffWW*Cc12*G/Dsqrt(E)                                               
                                                                                
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
c---- Centrifugal barrier calculation                                           
      Subroutine Centrif(MuWW, Z, L, VL)                                        
      Implicit None                                                             
      Real*8 Hbc, MuWW, VL, Z                                                   
      Integer*4 L                                                               
      Data Hbc/197.47d0/                                                        
                                                                                
      VL = Hbc*Hbc*L*(L+1)/(2.d0*MuWW*Z*Z)                                      
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Auxiliar function for the Coulomb energy of Gaudin.                       
      Double Precision Function F(X)                                            
      Implicit None                                                             
      Real*8 CosX, CosXh, Pi, SinX, SinXh, X                                    
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      CosX = Dcos(X)                                                            
      SinX = Dsin(X)                                                            
      CosXh= Dcos(X/2.d0)                                                       
      SinXh= Dsin(X/2.d0)                                                       
      F = 1.d0 - X*CosX/SinX - Pi*SinXh/(2.*CosXh)                              
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Auxiliar function for the Coulomb energy of Gaudin.                       
      Double Precision Function FL(X)                                           
      Implicit None                                                             
      Real*8 Pi, Tem1, Tem2, X, Xd, Xh                                          
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      Xd = 2.d0*X                                                               
      Xh = X/2.d0                                                               
      Tem1 = (Xd-Dsin(Xd))/(2.d0*Dsin(X)**2)                                    
      Tem2 = (Pi*Dtan(Xh)**2)/4.d0                                              
      FL = Tem1 - Tem2                                                          
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Auxiliar function for the Coulomb energy of Gaudin.                       
      Double Precision Function G(X)                                            
      Implicit None                                                             
      Real*8 Sin3, Tanh, Tanh2, X, Xh                                           
                                                                                
      Sin3 = Dsin(x)**3                                                         
      Xh = X/2.d0                                                               
      Tanh = Dsin(Xh)/Dcos(Xh)                                                  
      Tanh2= Tanh*Tanh                                                          
      G = 0.1d0*Tanh*(15.d0 + 10.d0*Tanh2 + 3.d0*Tanh2*Tanh2)+2.d0/Sin3         
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Calculation of the factor Epson of the Gaudin's Coulomb energy            
      Subroutine Epson(X1, X2, Eps)                                             
      Implicit None                                                             
      Real*8 Css1, Css2, Ctg1, Ctg2, Eps, F, FL, G, Sin1, Sin2, Pi,             
     +       Ter1, Ter2, Ter3, Ter4, X1, X12, X2                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      X12 = X1 + X2                                                             
      Sin1 = Dsin(X1)                                                           
      Sin2 = Dsin(X2)                                                           
      Ctg1 = Dcos(X1) / Sin1                                                    
      Ctg2 = Dcos(X2) / Sin2                                                    
      Css1 = 1.d0/(Sin1*Sin1)                                                   
      Css2 = 1.d0/(Sin2*Sin2)                                                   
                                                                                
      Ter1 = (Css2 - Css1)*( F(X2)*Css2 - F(X1)*Css1 )                          
      Ter2 = (Ctg2 + Ctg1)*(FL(X2)*Css2 +FL(X1)*Css1 )                          
      Ter3 = Css2*Css1*( F(X12) + Dsin(X12)**2/3.d0 )                           
      Ter4 = Pi*( G(X2) + G(X1) )/8.d0                                          
      Eps = Ter1 - Ter2 + Ter3 + Ter4                                           
                                                                                
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Coulomb energy according to the Gaudin's prescription                     
      Subroutine Ecoul(C1, Factor, X1, X2, Aa, Ec)                              
      Implicit None                                                             
      Real*8 Aa, C1, Ec, Eps, Factor, X1, X2                                    
                                                                                
      Call Epson(-X1, X2, Eps)           !   Epson calculates the Coulomb       
      Ec = Factor * C1 *Aa**5 * Eps       !   energy of two spherical fragments 
      Return                             !   dependent on the configuration.    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Surface energy calculation regarding the effective surface tension.       
      Subroutine Sup(Cs, X1, X2, Aa, Es)                                        
      Implicit None                                                             
      Real*8 Aa, Cs, Es, S1, S2, X1, X2                                         
                                                                                
      S1 = (1.d0 + Dcos(X1))/Dsin(X1)**2 !   The surface relations were taken   
      S2 = (1.d0 + Dcos(X2))/Dsin(X2)**2 !   from Gaudin's paper!               
      Es = Cs * Aa**2 * (S1 + S2)        !   The constant Cs keeps all nuclear  
      Return                             !   characteristics through the Q-     
      End                                !   value of the decay process.        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- The search for the root of the R2 equation.                               
      Subroutine NewtZ(Rp, R01, Z, R2, Loop_Exceed)                             
      Implicit None                                                             
                                                                                
      Real*8 Const, Epss, F, FL, Rp, R01, R2, R2T, Z                            
      Integer*4 Loop_Exceed                                                     
                                                                                
      Epss = 1.d-4                                                              
      R2t = Rp                                                                  
      Const = 16.d0*Rp**3 + Z**3 - R01**2*(8.d0*R01 + 6.d0*Z +                  
     .        3.d0*R01**2/Z)                                                    
      Do Loop_Exceed = 1, 10001                                                 
         F = 8.d0*R2t**3 + 3.d0*R2t**4/Z + 6.d0*R2t**2*                         
     .      (Z - R01**2/Z) - Const                                              
         FL=12.d0*( R2t**2*(2.d0*Z + R2t) + R2t*(Z**2 - R01**2) )               
         R2 = R2t - F/FL                                                        
         If (Dabs(R2-R2t) .lt. Epss) Goto 2                                     
         R2t = R2                                                               
      Enddo                                                                     
                                                                                
 2    Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Werner-Wheeler inertial coefficient taken from Poenaru et. al., Z. Phys.  
c---- A - Atomic Nuclei 333, 291 (1989).                                        
      Subroutine MuWWcl(Bconst1, Bconst2, R01, R2, Z, V1til, V2til, B)          
      Implicit None                                                             
                                                                                
      Real*8 B, Bconst1, Bconst2, B1, B11, B12, B2, Dens1, Dens2, M,            
     +       D1, H1, M1, M2, Pi, R, R01, R1, R1L, R2, V1,                       
     +       V1til, V2, V2til, Z, Z1L, Z2L                                      
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      R1 = R2          !  Exchanging our notation (R1 is the radius of the      
      R2 = R01         !  emitted cluster) with the Poenaru's notation (R1      
      Dens1 = Bconst2  !  is the radius of the larger fragment). At this        
      Dens2 = Bconst1  !  point we introduced two different mass densities      
      V1 = V2til       !  to allow the matching in the touching point           
      V2 = V1til       !  between the Werner-Wheeler inertial coefficient       
      M1 = Dens1 * V1  !  and the asymptotic free reduced mass.                 
      M2 = Dens2 * V2                                                           
      M = M1 + M2                                                               
                                                                                
      R = Z                                                                     
      D1 = 0.5d0 * (R + (R1**2 - R2**2)/R)                                      
      H1 = R1 - D1                                                              
      R1L = -0.5d0 * H1/R1                                                      
      Z1L = ( -M2 + Pi*( Dens1*R1*R1L*(R1+D1)**2 ) )/M                          
      Z2L = Z1L + 1.d0                                                          
                                                                                
      B11 = M1*Z1L**2 - 2.d0*Dens1*Pi*Z1L*R1*R1L*(R1+D1)**2                     
      B12 = Pi*Dens1*(R1*R1L)**2*( 2.d0*R1**2/H1 - 4.5d0*R1                     
     +      - 3.5d0*D1 + 6.d0*R1*Dlog(2.d0*R1/H1) )                             
      B1 = B11 + B12                                                            
      B2 = M2*Z2L**2                                                            
      B  = B1 + B2                                                              
                                                                                
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      Subroutine MuWWfl(N, Rt, Bt, Z, Mufl)                                     
      Real*8 Bt, Btmax, Btmin, Mufl, Rt, Rtmax, Rtmin, Z                        
      Integer*4 j, k, N                                                         
      Dimension Bt(N), Rt(N)                                                    
                                                                                
      k = 2                                                                     
      j = 1                                                                     
      If (Z .le. RT(1)) Then                                                    
         Mufl = Bt(1)                                                           
         j = 0                                                                  
      Endif                                                                     
      If (Z .ge. RT(N)) Then                                                    
         Mufl = BT(N)                                                           
         j = 0                                                                  
      Endif                                                                     
      Do While (j .ne. 0)                                                       
         If (Z .lt. RT(k)) Then                                                 
            Rtmax = RT(k)                                                       
            Rtmin = RT(k-1)                                                     
            BTmax = BT(k)                                                       
            BTmin = BT(k-1)                                                     
            Mufl = (Z-RTmin)/(RTmax-RTmin)*(BTmax-BTmin) + BTmin                
            j = 0                                                               
         Endif                                                                  
         k = k + 1                                                              
      Enddo                                                                     
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c---- Effective inertial coefficient calculation. Details are presented in the  
c---- Ms. Thesis of M. Goncalves.                                               
      Subroutine MuEffcl(E, Mas1, Rhomas1, Rhomas2, Volp0, Vol10, V1til,        
     +                   V2til, Alpha, Mucl)                                    
      Implicit none                                                             
      Real*8 Alpha, E, Fdpros, Mas1, Mucl, Mtil, Mtil1, Mtil2, Qtil,            
     +       Rhomas1, Rhomas2, Volp0, Vol10, V1til, V2til                       
                                                                                
      Fdpros = 1.d0 - V1til/Vol10                                               
      Qtil = Fdpros * E                                                         
      Mtil1 = ( Rhomas1 + Qtil/Volp0 )*V1til                                    
      Mtil2 = ( Rhomas2 + Qtil/Volp0 )*V2til + Fdpros * Mas1                    
      Mtil = Mtil1*Mtil2 / (Mtil1 + Mtil2)                                      
      Mucl = Alpha*Alpha * Mtil                                                 
                                                                                
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      Subroutine MuEfffl(N, Rt, Bt, Z, Mufl)                                    
      Real*8 Bt, Btmax, Btmin, Mufl, Rt, Rtmax, Rtmin, Z                        
      Integer*4 j, k, N                                                         
      Dimension Bt(N), Rt(N)                                                    
                                                                                
      k = 2                                                                     
      j = 1                                                                     
      If (Z .le. RT(1)) Then                                                    
         Mufl = Bt(1)                                                           
         j = 0                                                                  
      Endif                                                                     
      If (Z .ge. RT(N)) Then                                                    
         Mufl = BT(N)                                                           
         j = 0                                                                  
      Endif                                                                     
      Do While (j .ne. 0)                                                       
         If (Z .lt. RT(k)) Then                                                 
            Rtmax = RT(k)                                                       
            Rtmin = RT(k-1)                                                     
            BTmax = BT(k)                                                       
            BTmin = BT(k-1)                                                     
            Mufl = (Z-RTmin)/(RTmax-RTmin)*(BTmax-BTmin) + BTmin                
            j = 0                                                               
         Endif                                                                  
         k = k + 1                                                              
      Enddo                                                                     
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      Subroutine Configuration(N, Mass_kind, RP, R01, R02, Bconst1,             
     1           Bconst2, Mas, Mas1, Mu_reduced, AK, AK2, Rt, R1t,              
     2           R2t, Bt)                                                       
                                                                                
      Implicit None                                                             
      Real*8 AK, AK2, B, Bconst1, Bconst2, Bt, Dh1, Eps, H1, H1_initial,        
     1       H1_final, Mas, Mas1, Mu_reduced, Pi, R, Rp, Rt, R01, R02,          
     2       R1, R1t, R2, R2t                                                   
      Integer*4 j, Mass_kind, N                                                 
c------------------------------------------------------------------------------ 
      Dimension Bt(N), Rt(N), R1t(N), R2t(N)                                    
      Pi = 3.141592653589793d0                                                  
      Do j = 1, N      ! inicializacion de la tabla de B(R) para cada (A,Ae)    
         Bt(j)  = 0.0d0                                                         
         Rt(j)  = 0.0d0                                                         
         R1t(j) = 0.0d0                                                         
         R2t(j) = 0.0d0                                                         
      Enddo                                                                     
c------------------------------------------------------------------------------ 
C    definicion de valor maximo de h1, que corresponde a R=0(momento inicial)   
                                                                                
      H1_initial = Rp*(1.d0 - 2.d0*Dcos(Pi/3.d0 +                               
     +             Dacos(1.d0 - 2.d0 * (Ak-Ak2)/Ak)/3.d0))                      
      H1_final = 0.d0                                                           
      dH1 = (H1_final - H1_initial)/N                                           
      Eps = 1.d-10                                                              
      j = 1                                                                     
      Do H1 = H1_initial-Eps, H1_final+Eps, dH1                                 
         Call Geom(Mass_kind, R01, R02, Bconst1, Bconst2, Mas,                  
     +        Mas1, Mu_reduced, H1, B, R, R1, R2)                               
         Bt(j)  = B                                                             
         Rt(j)  = R                                                             
         R1t(j) = R1                                                            
         R2t(j) = R2                                                            
         j = j + 1                                                              
      Enddo                                                                     
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      Subroutine Geom(Mass_kind, R01, R02, BCONST1, BCONST2, Mas, Mas1,         
     +                Mu_reduced, H1, B, R, R1, R2)                             
      Implicit None                                                             
      Real*8 Aux, B, BCONST1, BCONST2, D1, D2, H1, H2, Mas, Mas1,               
     +       Mu_reduced, R, RON2, RS1, RS2, R01, R02, R1, R2                    
      Real*8 Cubic, Root                                                        
      Integer*4 Mass_kind                                                       
                                                                                
c------------------------------------------------------------------------------ 
C     CALCULO DE R1                                                             
C  Calculo de volumen y densidad reajustadas con la conservacion                
c  de la densidad de carga en los fragmentos filhos.                            
c                                                                               
      RS1=R02                                                                   
      RS2=R01                                                                   
                                                                                
c     Calculo de la raiz real de la ecuacion de tercer grado para R1            
      Aux = 4.d0*RS1**6 - 2.d0*(RS1*H1)**3                                      
      If (Aux .le. 0.d0) Aux = 0.d0                                             
      Root = Sqrt( Aux )                                                        
      Cubic = ( - 0.125d0*H1**3 + 0.5d0*RS1**3 +                                
     +            0.25d0*Root )**(1.d0/3.d0)                                    
      R1 = Cubic +  0.25d0 * H1*H1/Cubic                                        
      D1 = R1 - H1                                                              
c------------------------------------------------------------------------------ 
C     Calculo de Rhon2                                                          
      Ron2 = R1*R1 - D1*D1                                                      
c------------------------------------------------------------------------------ 
C     Calculo de H2                                                             
      If (Ron2 .le. 0.d0) Then                                                  
         H2 = 1.e-10                                                            
         R2 = R02-1.e-10                                                        
         D2 = R2                                                                
         R = R01 + R02                                                          
      Else                                                                      
         Aux = 16.d0*Rs2**6 + Ron2**3                                           
         If (Aux .le. 0.d0) Aux = 0.d0                                          
         Root = Sqrt( Aux )                                                     
         Cubic = ( Ron2**3*(32.d0*Rs2**6 + Ron2**3)/(512.d0*Rs2**9) +           
     +             Ron2**3*Root/(64.d0*Rs2**6) )**(1.d0/3.d0)                   
         H2 = Cubic + Ron2**4/( 64.d0*Rs2**6*Cubic ) + Ron2**2/                 
     +        (8.d0*Rs2**3)                                                     
c------------------------------------------------------------------------------ 
C        Calculo de R2 e D2                                                     
         R2 = 0.5d0 * (H2 + Ron2/H2)                                            
         D2 = R2 - H2                                                           
c------------------------------------------------------------------------------ 
C        Calculo de R e B                                                       
         R = D1 + D2                                                            
      Endif                                                                     
                                                                                
      If (Mass_kind .eq. 1) Call BWWfl (Bconst1, Bconst2, Mas, Mas1,            
     +   H1, H2, D1, D2, R1, R2, R, B)                                          
      If (Mass_kind .eq. 2) Call BEfffl(Mu_reduced, D1, R1, R2, R, B)           
                                                                                
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      Subroutine BWWfl (Bconst1, Bconst2, Mas, Mas1, H1, H2, D1, D2,            
     +                  R1, R2, R, B)                                           
      Implicit None                                                             
      Real*8 B, D1, D2, H1, H2, Mas, Mas1, R, R1, R2                            
      Real*8 Bconst1, Bconst2, BFL11, BFL12, BFL13, BFL21, BFL22, BFL23,        
     +       B1R, B2R, Pi, R1R1, R2R2, V1til, V2til, Z1_prime, Z2_prime         
                                                                                
      Data Pi/3.141592653589793d0/                                              
                                                                                
      R1R1 = -0.5d0 * H1 * (R2+D2)/(R+R1+R2)                                    
      R2R2 = -0.5d0 * H2 * (R1+D1)/(R+R1+R2)                                    
                                                                                
      V1til = Pi*(R1 + D1)**2 * (R1 + H1) / 3.d0                                
      V2til = Pi*(R2 + D2)**2 * (R2 + H2) / 3.d0                                
      Z1_prime = - Mas1/Mas                                                     
      Z2_prime = Z1_prime + 1.d0                                                
      BFL11 = Z1_prime*Z1_prime*V1til/Pi                                        
      BFL12 = - 2.d0*Z1_prime*R1R1*(R1 + D1)*(R1 + D1)                          
      BFL13 = R1R1*R1R1*(2.d0*R1*R1/H1 - 4.5d0*R1 - 3.5d0*D1 +                  
     +        6.d0*R1*Dlog(2.d0*R1/H1))                                         
      BFL21 = Z2_prime*Z2_prime*V2til/Pi                                        
      BFL22 =   2.d0*Z2_prime*R2R2*(R2 + D2)*(R2 + D2)                          
      BFL23 = R2R2*R2R2*(2.d0*R2*R2/H2 - 4.5d0*R2 - 3.5d0*D2 +                  
     +        6.d0*R2*Dlog(2.d0*R2/H2))                                         
                                                                                
      B1R = Pi*Bconst2*( BFL11 + BFL12 + BFL13 )                                
      B2R = Pi*Bconst1*( BFL21 + BFL22 + BFL23 )                                
                                                                                
      B = B1R + B2R                                                             
      Return                                                                    
      End                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      Subroutine BEfffl(Mu_reduced, D1, R1ww, R2ww, R, B)                       
      Implicit None                                                             
      Real*8 B, D1, R, R1, R1ww, R2, R2ww                                       
      Real*8 Alpha2, Beta, Chi, Mu_reduced, V1til, V2til, Z,                    
     +       Z1_bar, Z2_bar                                                     
      Real*8 Betadot_zdot, Chidot_zdot, Rdot_zdot_CM, Par1, Par2, Par3,         
     +       R1dot_zdot, R2dot_zdot, Z1_bar_dot, Z2_bar_dot                     
                                                                                
      Z = R                                                                     
      R1 = R2ww                                                                 
      R2 = R1ww                                                                 
      Chi = D1                                                                  
c------------------------------------------------------------------------------ 
      Beta = Z - Chi                                                            
      V1til = (R1 + Beta)**2 * (2.d0*R1 - Beta)/3.                              
      V2til = (R2 + Chi )**2 * (2.d0*R2 - Chi )/3.                              
      Z1_bar = 0.25d0*( R1**2 - Beta**2 )**2/V1til                              
      Z2_bar = 0.25d0*( R2**2 - Chi**2  )**2/V2til                              
c------------------------------------------------------------------------------ 
      Par1 = 6.d0*R1 + 4.d0*Z - 4.d0*Chi                                        
      Par2 = 5.d0*R1 + 3.d0*Z - 3.d0*Chi                                        
      Par3 = (6.d0*R2 + 4.d0*Chi)/(5.d0*R2 + 3.d0*Chi)                          
      R2dot_zdot = - ( Beta*Par1/R1 + Par2 ) /                                  
     +               ( Par1*(R2 + Z*Par3)/R1 + Par2*Par3 )                      
      R1dot_zdot = ( Beta + (R2 + Z*Par3)*R2dot_zdot )/R1                       
      Chidot_zdot = - Par3*R2dot_zdot                                           
      Betadot_zdot = 1.d0 - Chidot_zdot                                         
c------------------------------------------------------------------------------ 
      Z1_bar_dot = (R1**2 - Beta**2)*(R1*R1dot_zdot - Beta*Betadot_zdot)        
     +             /V1til                                                       
      Z2_bar_dot = (R2**2 - Chi**2 )*(R2*R2dot_zdot - Chi *Chidot_zdot)         
     +             /V2til                                                       
c------------------------------------------------------------------------------ 
      Rdot_zdot_CM = 1.d0 + Z1_bar_dot + Z2_bar_dot                             
      Alpha2 = Rdot_zdot_CM**2                                                  
      B = Alpha2*Mu_reduced                                                     
                                                                                
      Return                                                                    
      End                                                                       
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCC                                                                          
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! %%%                          IN.DAT   --- next 11 lines                     %%% 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!     234      92            1             1       0                              
! --------------------------------------------------                              
!      Ap      Zp     Key-like     Mass-kind       L                              
! --------------------------------------------------                              
! Ap -> Parent nucleus mass number                                                
! Zp -> Parent nucleus atomic number                                              
! L -> Total orbital momentum                                                     
! Key like = 1 -> VMAS                                                            
!            2 -> CMAS                                                            
! Mass kind = 1 -> Werner-Wheeler inertial coefficient                            
!             2 -> Effective inertial coefficient                                 
                                                                                
                                                                                
                                                                                
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! %%%                           OUT.DAT  --- next 579 lines                   %%% 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!  234.  92.   1.   0. .1000+101   40.0000                                        
!  234.  92.   1.   1. .1000+101   40.0000                                        
!  234.  92.   2.   1. .1000+101   40.0000                                        
!  234.  92.   3.   1. .1000+101   40.0000                                        
!  234.  92.   3.   2. .1000+101   40.0000                                        
!  234.  92.   4.   1. .1000+101   40.0000                                        
!  234.  92.   4.   2. .1156E-35   13.4155                                        
!  234.  92.   4.   3. .1000+101   40.0000                                        
!  234.  92.   5.   1. .1000+101   40.0000                                        
!  234.  92.   5.   2. .1000+101   40.0000                                        
!  234.  92.   5.   3. .1000+101   40.0000                                        
!  234.  92.   5.   4. .1000+101   40.0000                                        
!  234.  92.   6.   1. .1000+101   40.0000                                        
!  234.  92.   6.   2. .1000+101   40.0000                                        
!  234.  92.   6.   3. .1000+101   40.0000                                        
!  234.  92.   6.   4. .1000+101   40.0000                                        
!  234.  92.   7.   2. .1000+101   40.0000                                        
!  234.  92.   7.   3. .1000+101   40.0000                                        
!  234.  92.   7.   4. .1000+101   40.0000                                        
!  234.  92.   7.   5. .1000+101   40.0000                                        
!  234.  92.   8.   2. .1000+101   40.0000                                        
!  234.  92.   8.   3. .1000+101   40.0000                                        
!  234.  92.   8.   4. .1000+101   40.0000                                        
!  234.  92.   8.   5. .1000+101   40.0000                                        
!  234.  92.   8.   6. .1000+101   40.0000                                        
!  234.  92.   9.   2. .1000+101   40.0000                                        
!  234.  92.   9.   3. .1000+101   40.0000                                        
!  234.  92.   9.   4. .1000+101   40.0000                                        
!  234.  92.   9.   5. .1000+101   40.0000                                        
!  234.  92.   9.   6. .1000+101   40.0000                                        
!  234.  92.  10.   2. .1000+101   40.0000                                        
!  234.  92.  10.   3. .1000+101   40.0000                                        
!  234.  92.  10.   4. .1000+101   40.0000                                        
!  234.  92.  10.   5. .1000+101   40.0000                                        
!  234.  92.  10.   6. .1000+101   40.0000                                        
!  234.  92.  11.   3. .1000+101   40.0000                                        
!  234.  92.  11.   4. .1000+101   40.0000                                        
!  234.  92.  11.   5. .1000+101   40.0000                                        
!  234.  92.  11.   6. .1000+101   40.0000                                        
!  234.  92.  11.   7. .1000+101   40.0000                                        
!  234.  92.  12.   3. .1000+101   40.0000                                        
!  234.  92.  12.   4. .1000+101   40.0000                                        
!  234.  92.  12.   5. .1000+101   40.0000                                        
!  234.  92.  12.   6. .1000+101   40.0000                                        
!  234.  92.  12.   7. .1000+101   40.0000                                        
!  234.  92.  13.   4. .1000+101   40.0000                                        
!  234.  92.  13.   5. .1000+101   40.0000                                        
!  234.  92.  13.   6. .1000+101   40.0000                                        
!  234.  92.  13.   7. .1000+101   40.0000                                        
!  234.  92.  14.   4. .1000+101   40.0000                                        
!  234.  92.  14.   5. .1000+101   40.0000                                        
!  234.  92.  14.   6. .3116E-58   35.9851                                        
!  234.  92.  14.   7. .1000+101   40.0000                                        
!  234.  92.  15.   5. .1000+101   40.0000                                        
!  234.  92.  15.   6. .1000+101   40.0000                                        
!  234.  92.  15.   7. .1000+101   40.0000                                        
!  234.  92.  16.   5. .1000+101   40.0000                                        
!  234.  92.  16.   6. .1000+101   40.0000                                        
!  234.  92.  16.   7. .1000+101   40.0000                                        
!  234.  92.  16.   8. .1000+101   40.0000                                        
!  234.  92.  17.   5. .1000+101   40.0000                                        
!  234.  92.  17.   6. .1000+101   40.0000                                        
!  234.  92.  17.   7. .1000+101   40.0000                                        
!  234.  92.  17.   8. .1000+101   40.0000                                        
!  234.  92.  18.   5. .1000+101   40.0000                                        
!  234.  92.  18.   6. .1000+101   40.0000                                        
!  234.  92.  18.   7. .1000+101   40.0000                                        
!  234.  92.  18.   8. .3837E-60   37.8946                                        
!  234.  92.  18.   9. .1000+101   40.0000                                        
!  234.  92.  19.   5. .1000+101   40.0000                                        
!  234.  92.  19.   6. .1000+101   40.0000                                        
!  234.  92.  19.   7. .1000+101   40.0000                                        
!  234.  92.  19.   8. .1000+101   40.0000                                        
!  234.  92.  19.   9. .1000+101   40.0000                                        
!  234.  92.  20.   6. .1000+101   40.0000                                        
!  234.  92.  20.   7. .1000+101   40.0000                                        
!  234.  92.  20.   8. .1117E-58   36.4307                                        
!  234.  92.  20.   9. .1000+101   40.0000                                        
!  234.  92.  20.  10. .1000+101   40.0000                                        
!  234.  92.  21.   6. .1000+101   40.0000                                        
!  234.  92.  21.   7. .1000+101   40.0000                                        
!  234.  92.  21.   8. .1000+101   40.0000                                        
!  234.  92.  21.   9. .1000+101   40.0000                                        
!  234.  92.  21.  10. .1000+101   40.0000                                        
!  234.  92.  22.   6. .1000+101   40.0000                                        
!  234.  92.  22.   7. .1000+101   40.0000                                        
!  234.  92.  22.   8. .5838E-60   37.7123                                        
!  234.  92.  22.   9. .1000+101   40.0000                                        
!  234.  92.  22.  10. .6294E-55   32.6797                                        
!  234.  92.  23.   7. .1000+101   40.0000                                        
!  234.  92.  23.   8. .1000+101   40.0000                                        
!  234.  92.  23.   9. .9773E-58   35.4886                                        
!  234.  92.  23.  10. .9759E-56   33.4892                                        
!  234.  92.  24.   7. .1000+101   40.0000                                        
!  234.  92.  24.   8. .1000+101   40.0000                                        
!  234.  92.  24.   9. .1823E-61   39.2178                                        
!  234.  92.  24.  10. .6518E-48   25.6645                                        
!  234.  92.  24.  11. .1000+101   40.0000                                        
!  234.  92.  25.   8. .1000+101   40.0000                                        
!  234.  92.  25.   9. .1000+101   40.0000                                        
!  234.  92.  25.  10. .2616E-50   28.0609                                        
!  234.  92.  25.  11. .2628E-55   33.0589                                        
!  234.  92.  26.   8. .1000+101   40.0000                                        
!  234.  92.  26.   9. .1000+101   40.0000                                        
!  234.  92.  26.  10. .2677E-48   26.0509                                        
!  234.  92.  26.  11. .5481E-55   32.7398                                        
!  234.  92.  26.  12. .9343E-56   33.5081                                        
!  234.  92.  27.   9. .1000+101   40.0000                                        
!  234.  92.  27.  10. .8977E-60   37.5255                                        
!  234.  92.  27.  11. .7274E-51   28.6169                                        
!  234.  92.  27.  12. .4011E-54   31.8754                                        
!  234.  92.  28.   9. .1000+101   40.0000                                        
!  234.  92.  28.  10. .1000+101   40.0000                                        
!  234.  92.  28.  11. .3795E-57   34.8995                                        
!  234.  92.  28.  12. .4141E-47   24.8615                                        
!  234.  92.  29.   9. .1000+101   40.0000                                        
!  234.  92.  29.  10. .1000+101   40.0000                                        
!  234.  92.  29.  11. .1004E-61   39.4768                                        
!  234.  92.  29.  12. .3742E-52   29.9055                                        
!  234.  92.  29.  13. .4730E-55   32.8037                                        
!  234.  92.  30.  10. .1000+101   40.0000                                        
!  234.  92.  30.  11. .1000+101   40.0000                                        
!  234.  92.  30.  12. .8027E-52   29.5741                                        
!  234.  92.  30.  13. .1376E-56   34.3400                                        
!  234.  92.  31.  10. .1000+101   40.0000                                        
!  234.  92.  31.  11. .1000+101   40.0000                                        
!  234.  92.  31.  12. .7998E-61   38.5756                                        
!  234.  92.  31.  13. .5757E-55   32.7184                                        
!  234.  92.  32.  10. .1000+101   40.0000                                        
!  234.  92.  32.  11. .1000+101   40.0000                                        
!  234.  92.  32.  12. .3955E-62   39.8814                                        
!  234.  92.  32.  13. .9354E-60   37.5076                                        
!  234.  92.  32.  14. .6230E-52   29.6841                                        
!  234.  92.  33.  11. .1000+101   40.0000                                        
!  234.  92.  33.  12. .1000+101   40.0000                                        
!  234.  92.  33.  13. .2491E-61   39.0822                                        
!  234.  92.  33.  14. .5654E-56   33.7263                                        
!  234.  92.  34.  11. .1000+101   40.0000                                        
!  234.  92.  34.  12. .1000+101   40.0000                                        
!  234.  92.  34.  13. .1000+101   40.0000                                        
!  234.  92.  34.  14. .3727E-53   30.9072                                        
!  234.  92.  35.  11. .1000+101   40.0000                                        
!  234.  92.  35.  12. .1000+101   40.0000                                        
!  234.  92.  35.  13. .1000+101   40.0000                                        
!  234.  92.  35.  14. .6021E-61   38.6989                                        
!  234.  92.  35.  15. .2515E-58   36.0780                                        
!  234.  92.  36.  12. .1000+101   40.0000                                        
!  234.  92.  36.  13. .1000+101   40.0000                                        
!  234.  92.  36.  14. .1197E-60   38.4005                                        
!  234.  92.  36.  15. .1000+101   40.0000                                        
!  234.  92.  37.  12. .1000+101   40.0000                                        
!  234.  92.  37.  13. .1000+101   40.0000                                        
!  234.  92.  37.  14. .1000+101   40.0000                                        
!  234.  92.  37.  15. .1000+101   40.0000                                        
!  234.  92.  38.  13. .1000+101   40.0000                                        
!  234.  92.  38.  14. .1000+101   40.0000                                        
!  234.  92.  38.  15. .1000+101   40.0000                                        
!  234.  92.  38.  16. .1211E-59   37.3954                                        
!  234.  92.  39.  13. .1000+101   40.0000                                        
!  234.  92.  39.  14. .1000+101   40.0000                                        
!  234.  92.  39.  15. .1000+101   40.0000                                        
!  234.  92.  39.  16. .1000+101   40.0000                                        
!  234.  92.  40.  14. .1000+101   40.0000                                        
!  234.  92.  40.  15. .1000+101   40.0000                                        
!  234.  92.  40.  16. .1860E-60   38.2091                                        
!  234.  92.  41.  14. .1000+101   40.0000                                        
!  234.  92.  41.  15. .1000+101   40.0000                                        
!  234.  92.  41.  16. .1000+101   40.0000                                        
!  234.  92.  42.  14. .1000+101   40.0000                                        
!  234.  92.  42.  15. .1000+101   40.0000                                        
!  234.  92.  42.  16. .1000+101   40.0000                                        
!  234.  92.  42.  17. .1000+101   40.0000                                        
!  234.  92.  43.  15. .1000+101   40.0000                                        
!  234.  92.  43.  16. .1000+101   40.0000                                        
!  234.  92.  43.  17. .1000+101   40.0000                                        
!  234.  92.  44.  15. .1000+101   40.0000                                        
!  234.  92.  44.  16. .1000+101   40.0000                                        
!  234.  92.  44.  17. .1000+101   40.0000                                        
!  234.  92.  44.  18. .9414E-62   39.5048                                        
!  234.  92.  45.  15. .1000+101   40.0000                                        
!  234.  92.  45.  16. .1000+101   40.0000                                        
!  234.  92.  45.  17. .1000+101   40.0000                                        
!  234.  92.  45.  18. .1000+101   40.0000                                        
!  234.  92.  46.  15. .1000+101   40.0000                                        
!  234.  92.  46.  16. .1000+101   40.0000                                        
!  234.  92.  46.  17. .1000+101   40.0000                                        
!  234.  92.  46.  18. .5062E-60   37.7743                                        
!  234.  92.  46.  19. .1000+101   40.0000                                        
!  234.  92.  47.  16. .1000+101   40.0000                                        
!  234.  92.  47.  17. .1000+101   40.0000                                        
!  234.  92.  47.  18. .1000+101   40.0000                                        
!  234.  92.  47.  19. .7133E-62   39.6253                                        
!  234.  92.  48.  16. .1000+101   40.0000                                        
!  234.  92.  48.  17. .1000+101   40.0000                                        
!  234.  92.  48.  18. .1000+101   40.0000                                        
!  234.  92.  48.  19. .1000+101   40.0000                                        
!  234.  92.  48.  20. .2647E-58   36.0558                                        
!  234.  92.  49.  16. .1000+101   40.0000                                        
!  234.  92.  49.  17. .1000+101   40.0000                                        
!  234.  92.  49.  18. .1000+101   40.0000                                        
!  234.  92.  49.  19. .1000+101   40.0000                                        
!  234.  92.  49.  20. .7154E-60   37.6241                                        
!  234.  92.  50.  17. .1000+101   40.0000                                        
!  234.  92.  50.  18. .1000+101   40.0000                                        
!  234.  92.  50.  19. .1000+101   40.0000                                        
!  234.  92.  50.  20. .2150E-58   36.1462                                        
!  234.  92.  50.  21. .1000+101   40.0000                                        
!  234.  92.  51.  17. .1000+101   40.0000                                        
!  234.  92.  51.  18. .1000+101   40.0000                                        
!  234.  92.  51.  19. .1000+101   40.0000                                        
!  234.  92.  51.  20. .3007E-61   39.0005                                        
!  234.  92.  51.  21. .1000+101   40.0000                                        
!  234.  92.  52.  18. .1000+101   40.0000                                        
!  234.  92.  52.  19. .1000+101   40.0000                                        
!  234.  92.  52.  20. .1000+101   40.0000                                        
!  234.  92.  52.  21. .1000+101   40.0000                                        
!  234.  92.  53.  18. .1000+101   40.0000                                        
!  234.  92.  53.  19. .1000+101   40.0000                                        
!  234.  92.  53.  20. .1000+101   40.0000                                        
!  234.  92.  53.  21. .1000+101   40.0000                                        
!  234.  92.  53.  22. .1000+101   40.0000                                        
!  234.  92.  54.  19. .1000+101   40.0000                                        
!  234.  92.  54.  20. .1000+101   40.0000                                        
!  234.  92.  54.  21. .1000+101   40.0000                                        
!  234.  92.  54.  22. .1935E-60   38.1920                                        
!  234.  92.  55.  19. .1000+101   40.0000                                        
!  234.  92.  55.  20. .1000+101   40.0000                                        
!  234.  92.  55.  21. .1000+101   40.0000                                        
!  234.  92.  55.  22. .1000+101   40.0000                                        
!  234.  92.  55.  23. .1000+101   40.0000                                        
!  234.  92.  56.  20. .1000+101   40.0000                                        
!  234.  92.  56.  21. .1000+101   40.0000                                        
!  234.  92.  56.  22. .1000+101   40.0000                                        
!  234.  92.  56.  23. .1000+101   40.0000                                        
!  234.  92.  57.  20. .1000+101   40.0000                                        
!  234.  92.  57.  21. .1000+101   40.0000                                        
!  234.  92.  57.  22. .1000+101   40.0000                                        
!  234.  92.  57.  23. .1000+101   40.0000                                        
!  234.  92.  57.  24. .1000+101   40.0000                                        
!  234.  92.  58.  21. .1000+101   40.0000                                        
!  234.  92.  58.  22. .1000+101   40.0000                                        
!  234.  92.  58.  23. .1000+101   40.0000                                        
!  234.  92.  58.  24. .1000+101   40.0000                                        
!  234.  92.  59.  21. .1000+101   40.0000                                        
!  234.  92.  59.  22. .1000+101   40.0000                                        
!  234.  92.  59.  23. .1000+101   40.0000                                        
!  234.  92.  59.  24. .1000+101   40.0000                                        
!  234.  92.  59.  25. .1000+101   40.0000                                        
!  234.  92.  60.  22. .1000+101   40.0000                                        
!  234.  92.  60.  23. .1000+101   40.0000                                        
!  234.  92.  60.  24. .1000+101   40.0000                                        
!  234.  92.  60.  25. .1000+101   40.0000                                        
!  234.  92.  61.  22. .1000+101   40.0000                                        
!  234.  92.  61.  23. .1000+101   40.0000                                        
!  234.  92.  61.  24. .1000+101   40.0000                                        
!  234.  92.  61.  25. .1000+101   40.0000                                        
!  234.  92.  61.  26. .1000+101   40.0000                                        
!  234.  92.  62.  23. .1000+101   40.0000                                        
!  234.  92.  62.  24. .1000+101   40.0000                                        
!  234.  92.  62.  25. .1000+101   40.0000                                        
!  234.  92.  62.  26. .1000+101   40.0000                                        
!  234.  92.  63.  23. .1000+101   40.0000                                        
!  234.  92.  63.  24. .1000+101   40.0000                                        
!  234.  92.  63.  25. .1000+101   40.0000                                        
!  234.  92.  63.  26. .1000+101   40.0000                                        
!  234.  92.  63.  27. .1000+101   40.0000                                        
!  234.  92.  64.  24. .1000+101   40.0000                                        
!  234.  92.  64.  25. .1000+101   40.0000                                        
!  234.  92.  64.  26. .1578E-61   39.2805                                        
!  234.  92.  64.  27. .1000+101   40.0000                                        
!  234.  92.  65.  24. .1000+101   40.0000                                        
!  234.  92.  65.  25. .1000+101   40.0000                                        
!  234.  92.  65.  26. .1000+101   40.0000                                        
!  234.  92.  65.  27. .1000+101   40.0000                                        
!  234.  92.  65.  28. .1000+101   40.0000                                        
!  234.  92.  66.  25. .1000+101   40.0000                                        
!  234.  92.  66.  26. .1425E-61   39.3247                                        
!  234.  92.  66.  27. .1000+101   40.0000                                        
!  234.  92.  66.  28. .1000+101   40.0000                                        
!  234.  92.  67.  25. .1000+101   40.0000                                        
!  234.  92.  67.  26. .1000+101   40.0000                                        
!  234.  92.  67.  27. .1000+101   40.0000                                        
!  234.  92.  67.  28. .1000+101   40.0000                                        
!  234.  92.  67.  29. .1000+101   40.0000                                        
!  234.  92.  68.  26. .1000+101   40.0000                                        
!  234.  92.  68.  27. .1000+101   40.0000                                        
!  234.  92.  68.  28. .6446E-60   37.6693                                        
!  234.  92.  68.  29. .1000+101   40.0000                                        
!  234.  92.  69.  26. .1000+101   40.0000                                        
!  234.  92.  69.  27. .1000+101   40.0000                                        
!  234.  92.  69.  28. .2259E-61   39.1246                                        
!  234.  92.  69.  29. .1000+101   40.0000                                        
!  234.  92.  69.  30. .1000+101   40.0000                                        
!  234.  92.  70.  27. .1000+101   40.0000                                        
!  234.  92.  70.  28. .1213E-58   36.3947                                        
!  234.  92.  70.  29. .1000+101   40.0000                                        
!  234.  92.  70.  30. .1000+101   40.0000                                        
!  234.  92.  71.  27. .1000+101   40.0000                                        
!  234.  92.  71.  28. .4489E-61   38.8265                                        
!  234.  92.  71.  29. .1000+101   40.0000                                        
!  234.  92.  71.  30. .1000+101   40.0000                                        
!  234.  92.  71.  31. .1000+101   40.0000                                        
!  234.  92.  72.  27. .1000+101   40.0000                                        
!  234.  92.  72.  28. .2720E-59   37.0440                                        
!  234.  92.  72.  29. .1000+101   40.0000                                        
!  234.  92.  72.  30. .1000+101   40.0000                                        
!  234.  92.  72.  31. .1000+101   40.0000                                        
!  234.  92.  73.  28. .1000+101   40.0000                                        
!  234.  92.  73.  29. .3063E-61   38.9924                                        
!  234.  92.  73.  30. .1000+101   40.0000                                        
!  234.  92.  73.  31. .1000+101   40.0000                                        
!  234.  92.  73.  32. .1000+101   40.0000                                        
!  234.  92.  74.  28. .1000+101   40.0000                                        
!  234.  92.  74.  29. .1000+101   40.0000                                        
!  234.  92.  74.  30. .5011E-59   36.7787                                        
!  234.  92.  74.  31. .1000+101   40.0000                                        
!  234.  92.  74.  32. .1000+101   40.0000                                        
!  234.  92.  75.  28. .1000+101   40.0000                                        
!  234.  92.  75.  29. .1000+101   40.0000                                        
!  234.  92.  75.  30. .6642E-61   38.6563                                        
!  234.  92.  75.  31. .1000+101   40.0000                                        
!  234.  92.  75.  32. .1000+101   40.0000                                        
!  234.  92.  75.  33. .1000+101   40.0000                                        
!  234.  92.  76.  28. .1000+101   40.0000                                        
!  234.  92.  76.  29. .1000+101   40.0000                                        
!  234.  92.  76.  30. .5930E-58   35.7056                                        
!  234.  92.  76.  31. .1000+101   40.0000                                        
!  234.  92.  76.  32. .1000+101   40.0000                                        
!  234.  92.  76.  33. .1000+101   40.0000                                        
!  234.  92.  77.  28. .1000+101   40.0000                                        
!  234.  92.  77.  29. .1000+101   40.0000                                        
!  234.  92.  77.  30. .2125E-60   38.1513                                        
!  234.  92.  77.  31. .9825E-61   38.4863                                        
!  234.  92.  77.  32. .1000+101   40.0000                                        
!  234.  92.  77.  33. .1000+101   40.0000                                        
!  234.  92.  77.  34. .1000+101   40.0000                                        
!  234.  92.  78.  28. .1000+101   40.0000                                        
!  234.  92.  78.  29. .1000+101   40.0000                                        
!  234.  92.  78.  30. .5155E-59   36.7664                                        
!  234.  92.  78.  31. .3916E-61   38.8858                                        
!  234.  92.  78.  32. .2557E-59   37.0709                                        
!  234.  92.  78.  33. .1000+101   40.0000                                        
!  234.  92.  78.  34. .1000+101   40.0000                                        
!  234.  92.  79.  29. .1000+101   40.0000                                        
!  234.  92.  79.  30. .3839E-62   39.8944                                        
!  234.  92.  79.  31. .2410E-59   37.0967                                        
!  234.  92.  79.  32. .2034E-59   37.1702                                        
!  234.  92.  79.  33. .1000+101   40.0000                                        
!  234.  92.  79.  34. .1000+101   40.0000                                        
!  234.  92.  79.  35. .1000+101   40.0000                                        
!  234.  92.  80.  29. .1000+101   40.0000                                        
!  234.  92.  80.  30. .3822E-62   39.8963                                        
!  234.  92.  80.  31. .7242E-62   39.6188                                        
!  234.  92.  80.  32. .4092E-56   33.8667                                        
!  234.  92.  80.  33. .1000+101   40.0000                                        
!  234.  92.  80.  34. .1000+101   40.0000                                        
!  234.  92.  80.  35. .1000+101   40.0000                                        
!  234.  92.  81.  30. .1000+101   40.0000                                        
!  234.  92.  81.  31. .1519E-60   38.2971                                        
!  234.  92.  81.  32. .5743E-58   35.7195                                        
!  234.  92.  81.  33. .2637E-60   38.0576                                        
!  234.  92.  81.  34. .1000+101   40.0000                                        
!  234.  92.  81.  35. .1000+101   40.0000                                        
!  234.  92.  81.  36. .1000+101   40.0000                                        
!  234.  92.  82.  30. .1000+101   40.0000                                        
!  234.  92.  82.  31. .1000+101   40.0000                                        
!  234.  92.  82.  32. .1454E-55   33.3161                                        
!  234.  92.  82.  33. .8899E-61   38.5293                                        
!  234.  92.  82.  34. .2511E-60   38.0787                                        
!  234.  92.  82.  35. .1000+101   40.0000                                        
!  234.  92.  82.  36. .1000+101   40.0000                                        
!  234.  92.  83.  31. .1000+101   40.0000                                        
!  234.  92.  83.  32. .3273E-60   37.9636                                        
!  234.  92.  83.  33. .1153E-57   35.4168                                        
!  234.  92.  83.  34. .1607E-60   38.2725                                        
!  234.  92.  83.  35. .1000+101   40.0000                                        
!  234.  92.  83.  36. .1000+101   40.0000                                        
!  234.  92.  83.  37. .1000+101   40.0000                                        
!  234.  92.  84.  31. .1000+101   40.0000                                        
!  234.  92.  84.  32. .3679E-60   37.9129                                        
!  234.  92.  84.  33. .6221E-61   38.6848                                        
!  234.  92.  84.  34. .1048E-55   33.4582                                        
!  234.  92.  84.  35. .1000+101   40.0000                                        
!  234.  92.  84.  36. .1000+101   40.0000                                        
!  234.  92.  84.  37. .1000+101   40.0000                                        
!  234.  92.  85.  32. .1000+101   40.0000                                        
!  234.  92.  85.  33. .1607E-60   38.2726                                        
!  234.  92.  85.  34. .8305E-58   35.5593                                        
!  234.  92.  85.  35. .7379E-60   37.6106                                        
!  234.  92.  85.  36. .1000+101   40.0000                                        
!  234.  92.  85.  37. .1000+101   40.0000                                        
!  234.  92.  86.  32. .1000+101   40.0000                                        
!  234.  92.  86.  33. .1000+101   40.0000                                        
!  234.  92.  86.  34. .7532E-56   33.6017                                        
!  234.  92.  86.  35. .4498E-61   38.8255                                        
!  234.  92.  86.  36. .1645E-60   38.2624                                        
!  234.  92.  86.  37. .1000+101   40.0000                                        
!  234.  92.  87.  33. .1000+101   40.0000                                        
!  234.  92.  87.  34. .1576E-58   36.2810                                        
!  234.  92.  87.  35. .1967E-58   36.1849                                        
!  234.  92.  87.  36. .1507E-59   37.3006                                        
!  234.  92.  87.  37. .1000+101   40.0000                                        
!  234.  92.  87.  38. .1000+101   40.0000                                        
!  234.  92.  88.  33. .1000+101   40.0000                                        
!  234.  92.  88.  34. .1354E-57   35.3471                                        
!  234.  92.  88.  35. .7006E-60   37.6332                                        
!  234.  92.  88.  36. .1447E-56   34.3180                                        
!  234.  92.  88.  37. .1000+101   40.0000                                        
!  234.  92.  88.  38. .1000+101   40.0000                                        
!  234.  92.  89.  33. .1000+101   40.0000                                        
!  234.  92.  89.  34. .3614E-61   38.9207                                        
!  234.  92.  89.  35. .4959E-58   35.7832                                        
!  234.  92.  89.  36. .1221E-56   34.3919                                        
!  234.  92.  89.  37. .1000+101   40.0000                                        
!  234.  92.  89.  38. .1000+101   40.0000                                        
!  234.  92.  90.  34. .4887E-61   38.7895                                        
!  234.  92.  90.  35. .1506E-60   38.3007                                        
!  234.  92.  90.  36. .1938E-54   32.1912                                        
!  234.  92.  90.  37. .2089E-61   39.1588                                        
!  234.  92.  90.  38. .1000+101   40.0000                                        
!  234.  92.  90.  39. .1000+101   40.0000                                        
!  234.  92.  91.  34. .1000+101   40.0000                                        
!  234.  92.  91.  35. .2156E-60   38.1451                                        
!  234.  92.  91.  36. .3141E-56   33.9815                                        
!  234.  92.  91.  37. .3479E-58   35.9371                                        
!  234.  92.  91.  38. .1676E-61   39.2544                                        
!  234.  92.  91.  39. .1000+101   40.0000                                        
!  234.  92.  92.  34. .1000+101   40.0000                                        
!  234.  92.  92.  35. .1000+101   40.0000                                        
!  234.  92.  92.  36. .1052E-54   32.4566                                        
!  234.  92.  92.  37. .2050E-58   36.1669                                        
!  234.  92.  92.  38. .9326E-57   34.5089                                        
!  234.  92.  92.  39. .1000+101   40.0000                                        
!  234.  92.  92.  40. .1000+101   40.0000                                        
!  234.  92.  93.  35. .1000+101   40.0000                                        
!  234.  92.  93.  36. .4088E-58   35.8671                                        
!  234.  92.  93.  37. .2397E-56   34.0989                                        
!  234.  92.  93.  38. .9668E-57   34.4933                                        
!  234.  92.  93.  39. .1000+101   40.0000                                        
!  234.  92.  93.  40. .1000+101   40.0000                                        
!  234.  92.  94.  35. .1000+101   40.0000                                        
!  234.  92.  94.  36. .2180E-57   35.1401                                        
!  234.  92.  94.  37. .3688E-58   35.9118                                        
!  234.  92.  94.  38. .7932E-53   30.5792                                        
!  234.  92.  94.  39. .7455E-61   38.6061                                        
!  234.  92.  94.  40. .1000+101   40.0000                                        
!  234.  92.  95.  36. .1416E-61   39.3276                                        
!  234.  92.  95.  37. .4200E-57   34.8553                                        
!  234.  92.  95.  38. .4195E-54   31.8559                                        
!  234.  92.  95.  39. .1254E-56   34.3805                                        
!  234.  92.  95.  40. .1000+101   40.0000                                        
!  234.  92.  95.  41. .1000+101   40.0000                                        
!  234.  92.  96.  36. .3228E-61   38.9696                                        
!  234.  92.  96.  37. .4830E-60   37.7947                                        
!  234.  92.  96.  38. .1718E-51   29.2436                                        
!  234.  92.  96.  39. .5584E-56   33.7317                                        
!  234.  92.  96.  40. .1385E-56   34.3372                                        
!  234.  92.  96.  41. .1000+101   40.0000                                        
!  234.  92.  97.  36. .1000+101   40.0000                                        
!  234.  92.  97.  37. .4024E-59   36.8740                                        
!  234.  92.  97.  38. .1027E-53   31.4671                                        
!  234.  92.  97.  39. .1477E-53   31.3092                                        
!  234.  92.  97.  40. .2676E-55   33.0511                                        
!  234.  92.  97.  41. .1000+101   40.0000                                        
!  234.  92.  97.  42. .1000+101   40.0000                                        
!  234.  92.  98.  37. .1000+101   40.0000                                        
!  234.  92.  98.  38. .1398E-51   29.3332                                        
!  234.  92.  98.  39. .1579E-54   32.2803                                        
!  234.  92.  98.  40. .1183E-51   29.4057                                        
!  234.  92.  98.  41. .9367E-62   39.5070                                        
!  234.  92.  98.  42. .1000+101   40.0000                                        
!  234.  92.  99.  37. .1000+101   40.0000                                        
!  234.  92.  99.  38. .7799E-57   34.5866                                        
!  234.  92.  99.  39. .3329E-52   29.9564                                        
!  234.  92.  99.  40. .8444E-52   29.5521                                        
!  234.  92.  99.  41. .3819E-57   34.8966                                        
!  234.  92.  99.  42. .1000+101   40.0000                                        
!  234.  92. 100.  37. .1000+101   40.0000                                        
!  234.  92. 100.  38. .4135E-57   34.8621                                        
!  234.  92. 100.  39. .2085E-55   33.1595                                        
!  234.  92. 100.  40. .6062E-48   25.6960                                        
!  234.  92. 100.  41. .6113E-55   32.6923                                        
!  234.  92. 100.  42. .3035E-57   34.9965                                        
!  234.  92. 100.  43. .1000+101   40.0000                                        
!  234.  92. 101.  37. .1000+101   40.0000                                        
!  234.  92. 101.  38. .1000+101   40.0000                                        
!  234.  92. 101.  39. .5817E-56   33.7139                                        
!  234.  92. 101.  40. .6668E-51   28.6546                                        
!  234.  92. 101.  41. .2060E-50   28.1647                                        
!  234.  92. 101.  42. .2494E-55   33.0816                                        
!  234.  92. 101.  43. .1000+101   40.0000                                        
!  234.  92. 102.  37. .1000+101   40.0000                                        
!  234.  92. 102.  38. .1000+101   40.0000                                        
!  234.  92. 102.  39. .9258E-60   37.5121                                        
!  234.  92. 102.  40. .2618E-50   28.0606                                        
!  234.  92. 102.  41. .1599E-52   30.2746                                        
!  234.  92. 102.  42. .8427E-49   26.5529                                        
!  234.  92. 102.  43. .1000+101   40.0000                                        
!  234.  92. 103.  38. .1000+101   40.0000                                        
!  234.  92. 103.  39. .1790E-61   39.2257                                        
!  234.  92. 103.  40. .3396E-54   31.9476                                        
!  234.  92. 103.  41. .4437E-51   28.8315                                        
!  234.  92. 103.  42. .5060E-51   28.7745                                        
!  234.  92. 103.  43. .4423E-58   35.8329                                        
!  234.  92. 104.  38. .1000+101   40.0000                                        
!  234.  92. 104.  39. .1000+101   40.0000                                        
!  234.  92. 104.  40. .4402E-54   31.8349                                        
!  234.  92. 104.  41. .3198E-54   31.9738                                        
!  234.  92. 104.  42. .2307E-48   26.1156                                        
!  234.  92. 104.  43. .1756E-58   36.2340                                        
!  234.  92. 104.  44. .1000+101   40.0000                                        
!  234.  92. 105.  39. .1000+101   40.0000                                        
!  234.  92. 105.  40. .3583E-59   36.9244                                        
!  234.  92. 105.  41. .3093E-53   30.9882                                        
!  234.  92. 105.  42. .2405E-51   29.0975                                        
!  234.  92. 105.  43. .3153E-55   32.9799                                        
!  234.  92. 105.  44. .1000+101   40.0000                                        
!  234.  92. 106.  39. .1000+101   40.0000                                        
!  234.  92. 106.  40. .5466E-60   37.7409                                        
!  234.  92. 106.  41. .6858E-58   35.6424                                        
!  234.  92. 106.  42. .1713E-49   27.2448                                        
!  234.  92. 106.  43. .1503E-56   34.3016                                        
!  234.  92. 106.  44. .1065E-57   35.4514                                        
!  234.  92. 107.  40. .1000+101   40.0000                                        
!  234.  92. 107.  41. .9509E-58   35.5005                                        
!  234.  92. 107.  42. .4218E-53   30.8535                                        
!  234.  92. 107.  43. .2877E-54   32.0197                                        
!  234.  92. 107.  44. .4398E-59   36.8353                                        
!  234.  92. 107.  45. .1000+101   40.0000                                        
!  234.  92. 108.  40. .1000+101   40.0000                                        
!  234.  92. 108.  41. .1000+101   40.0000                                        
!  234.  92. 108.  42. .1193E-52   30.4020                                        
!  234.  92. 108.  43. .5247E-57   34.7587                                        
!  234.  92. 108.  44. .6381E-55   32.6737                                        
!  234.  92. 108.  45. .1000+101   40.0000                                        
!  234.  92. 109.  41. .1000+101   40.0000                                        
!  234.  92. 109.  42. .6009E-57   34.6998                                        
!  234.  92. 109.  43. .1630E-55   33.2663                                        
!  234.  92. 109.  44. .5450E-57   34.7422                                        
!  234.  92. 109.  45. .1000+101   40.0000                                        
!  234.  92. 110.  41. .1000+101   40.0000                                        
!  234.  92. 110.  42. .2700E-56   34.0472                                        
!  234.  92. 110.  43. .8902E-59   36.5291                                        
!  234.  92. 110.  44. .6773E-54   31.6478                                        
!  234.  92. 110.  45. .1000+101   40.0000                                        
!  234.  92. 111.  42. .4463E-62   39.8290                                        
!  234.  92. 111.  43. .1404E-57   35.3313                                        
!  234.  92. 111.  44. .4287E-57   34.8464                                        
!  234.  92. 111.  45. .1360E-59   37.3450                                        
!  234.  92. 111.  46. .1000+101   40.0000                                        
!  234.  92. 112.  42. .4059E-62   39.8702                                        
!  234.  92. 112.  43. .4854E-62   39.7925                                        
!  234.  92. 112.  44. .1241E-54   32.3849                                        
!  234.  92. 112.  45. .2183E-60   38.1395                                        
!  234.  92. 112.  46. .2942E-60   38.0100                                        
!  234.  92. 113.  42. .1000+101   40.0000                                        
!  234.  92. 113.  43. .1215E-61   39.3942                                        
!  234.  92. 113.  44. .4137E-58   35.8620                                        
!  234.  92. 113.  45. .3969E-58   35.8799                                        
!  234.  92. 113.  46. .1336E-61   39.3530                                        
!  234.  92. 113.  47. .1000+101   40.0000                                        
!  234.  92. 114.  43. .1000+101   40.0000                                        
!  234.  92. 114.  44. .3876E-56   33.8902                                        
!  234.  92. 114.  45. .4679E-59   36.8084                                        
!  234.  92. 114.  46. .3357E-57   34.9526                                        
!  234.  92. 114.  47. .1000+101   40.0000                                        
!  234.  92. 115.  43. .1000+101   40.0000                                        
!  234.  92. 115.  44. .4505E-61   38.8249                                        
!  234.  92. 115.  45. .1189E-58   36.4034                                        
!  234.  92. 115.  46. .2387E-59   37.1007                                        
!  234.  92. 115.  47. .1000+101   40.0000                                        
!  234.  92. 116.  44. .9461E-60   37.5027                                        
!  234.  92. 116.  45. .2498E-60   38.0811                                        
!  234.  92. 116.  46. .9070E-56   33.5210                                        
!  234.  92. 116.  47. .1000+101   40.0000                                        
!  234.  92. 116.  48. .1000+101   40.0000                                        
!  234.  92. 117.  44. .1000+101   40.0000                                        
!  234.  92. 117.  45. .5762E-60   37.7180                                        
!  234.  92. 117.  46. .1598E-58   36.2750                                        
                                                                                
                                                                                
                                                                                
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! %%%                          MASS_EXC.DAT  --- next 2930 lines              %%% 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!     1    1       7.288969                                                       
!     2    1      13.135720                                                       
!     3    1      14.949794                                                       
!     3    2      14.931204                                                       
!     4    1      25.927784                                                       
!     4    2       2.424911                                                       
!     4    3      25.320173                                                       
!     5    1      36.833979                                                       
!     5    2      11.386234                                                       
!     5    3      11.678880                                                       
!     5    4      37.996000                                                       
!     6    1      41.863763                                                       
!     6    2      17.594123                                                       
!     6    3      14.086312                                                       
!     6    4      18.374465                                                       
!     7    2      26.110264                                                       
!     7    3      14.907673                                                       
!     7    4      15.769489                                                       
!     7    5      27.867864                                                       
!     8    2      31.597983                                                       
!     8    3      20.946195                                                       
!     8    4       4.941662                                                       
!     8    5      22.921002                                                       
!     8    6      35.094056                                                       
!     9    2      40.818362                                                       
!     9    3      24.953903                                                       
!     9    4      11.347584                                                       
!     9    5      12.415701                                                       
!     9    6      28.913650                                                       
!    10    2      48.810011                                                       
!    10    3      33.050226                                                       
!    10    4      12.606577                                                       
!    10    5      12.050761                                                       
!    10    6      15.698568                                                       
!    10    7      39.699000                                                       
!    11    3      40.795860                                                       
!    11    4      20.173970                                                       
!    11    5       8.667984                                                       
!    11    6      10.650531                                                       
!    11    7      24.960520                                                       
!    12    3      50.096000                                                       
!    12    4      25.076402                                                       
!    12    5      13.368901                                                       
!    12    6        .000000                                                       
!    12    7      17.338083                                                       
!    12    8      32.047838                                                       
!    13    4      33.658445                                                       
!    13    5      16.562209                                                       
!    13    6       3.125011                                                       
!    13    7       5.345456                                                       
!    13    8      23.110735                                                       
!    14    4      39.882396                                                       
!    14    5      23.663730                                                       
!    14    6       3.019892                                                       
!    14    7       2.863417                                                       
!    14    8       8.006456                                                       
!    14    9      33.608000                                                       
!    15    5      28.966936                                                       
!    15    6       9.873143                                                       
!    15    7        .101438                                                       
!    15    8       2.855388                                                       
!    15    9      16.777001                                                       
!    16    5      37.081686                                                       
!    16    6      13.694117                                                       
!    16    7       5.683432                                                       
!    16    8      -4.736998                                                       
!    16    9      10.680257                                                       
!    16   10      23.992401                                                       
!    17    5      43.716310                                                       
!    17    6      21.036589                                                       
!    17    7       7.870819                                                       
!    17    8       -.809002                                                       
!    17    9       1.951701                                                       
!    17   10      16.485173                                                       
!    18    5      52.322000                                                       
!    18    6      24.924036                                                       
!    18    7      13.117136                                                       
!    18    8       -.782064                                                       
!    18    9        .873431                                                       
!    18   10       5.306782                                                       
!    18   11      25.318000                                                       
!    19    5      59.364000                                                       
!    19    6      32.833383                                                       
!    19    7      15.860449                                                       
!    19    8       3.333565                                                       
!    19    9      -1.487405                                                       
!    19   10       1.751059                                                       
!    19   11      12.928622                                                       
!    20    6      37.560063                                                       
!    20    7      21.766492                                                       
!    20    8       3.796909                                                       
!    20    9       -.017396                                                       
!    20   10      -7.041930                                                       
!    20   11       6.844859                                                       
!    20   12      17.570530                                                       
!    21    6      45.960000                                                       
!    21    7      25.231909                                                       
!    21    8       8.061736                                                       
!    21    9       -.047580                                                       
!    21   10      -5.731720                                                       
!    21   11      -2.184261                                                       
!    21   12      10.911681                                                       
!    21   13      26.119000                                                       
!    22    6      52.583000                                                       
!    22    7      32.080890                                                       
!    22    8       9.284346                                                       
!    22    9       2.793783                                                       
!    22   10      -8.024344                                                       
!    22   11      -5.182104                                                       
!    22   12       -.396766                                                       
!    22   13      18.183000                                                       
!    22   14      32.164000                                                       
!    23    7      37.735000                                                       
!    23    8      14.616373                                                       
!    23    9       3.329518                                                       
!    23   10      -5.153641                                                       
!    23   11      -9.529485                                                       
!    23   12      -5.472666                                                       
!    23   13       6.767210                                                       
!    23   14      23.772000                                                       
!    24    7      47.040000                                                       
!    24    8      18.974457                                                       
!    24    9       7.544514                                                       
!    24   10      -5.947519                                                       
!    24   11      -8.417601                                                       
!    24   12     -13.933381                                                       
!    24   13       -.055041                                                       
!    24   14      10.754759                                                       
!    24   15      31.997000                                                       
!    25    8      27.144000                                                       
!    25    9      11.266384                                                       
!    25   10      -2.058696                                                       
!    25   11      -9.357459                                                       
!    25   12     -13.192726                                                       
!    25   13      -8.915742                                                       
!    25   14       3.825310                                                       
!    25   15      18.872000                                                       
!    26    8      35.164000                                                       
!    26    9      18.288165                                                       
!    26   10        .429882                                                       
!    26   11      -6.902465                                                       
!    26   12     -16.214476                                                       
!    26   13     -12.210339                                                       
!    26   14      -7.144618                                                       
!    26   15      10.973000                                                       
!    26   16      25.970000                                                       
!    27    9      25.050028                                                       
!    27   10       7.093512                                                       
!    27   11      -5.580857                                                       
!    27   12     -14.586503                                                       
!    27   13     -17.196829                                                       
!    27   14     -12.384431                                                       
!    27   15       -.752978                                                       
!    27   16      17.507000                                                       
!    28    9      33.226000                                                       
!    28   10      11.278595                                                       
!    28   11      -1.033576                                                       
!    28   12     -15.018752                                                       
!    28   13     -16.850552                                                       
!    28   14     -21.492793                                                       
!    28   15      -7.161018                                                       
!    28   16       4.073107                                                       
!    28   17      26.557000                                                       
!    29    9      40.296000                                                       
!    29   10      18.020589                                                       
!    29   11       2.618710                                                       
!    29   12     -10.661187                                                       
!    29   13     -18.215504                                                       
!    29   14     -21.895025                                                       
!    29   15     -16.951907                                                       
!    29   16      -3.158878                                                       
!    29   17      13.143000                                                       
!    30   10      22.236621                                                       
!    30   11       8.594416                                                       
!    30   12      -8.882233                                                       
!    30   13     -15.872372                                                       
!    30   14     -24.432881                                                       
!    30   15     -20.200556                                                       
!    30   16     -14.062806                                                       
!    30   17       4.443000                                                       
!    30   18      20.083000                                                       
!    31   10      30.842000                                                       
!    31   11      12.663760                                                       
!    31   12      -3.215089                                                       
!    31   13     -14.954181                                                       
!    31   14     -22.948958                                                       
!    31   15     -24.440991                                                       
!    31   16     -19.044932                                                       
!    31   17      -7.064436                                                       
!    31   18      11.296000                                                       
!    32   10      37.176000                                                       
!    32   11      18.303661                                                       
!    32   12       -.795599                                                       
!    32   13     -11.062068                                                       
!    32   14     -24.080859                                                       
!    32   15     -24.305318                                                       
!    32   16     -26.015981                                                       
!    32   17     -13.330694                                                       
!    32   18      -2.179081                                                       
!    32   19      20.418000                                                       
!    33   11      25.509891                                                       
!    33   12       5.204233                                                       
!    33   13      -8.504923                                                       
!    33   14     -20.492381                                                       
!    33   15     -26.337727                                                       
!    33   16     -26.586235                                                       
!    33   17     -21.003512                                                       
!    33   18      -9.381336                                                       
!    33   19       6.763000                                                       
!    34   11      32.509000                                                       
!    34   12       8.450922                                                       
!    34   13      -2.862243                                                       
!    34   14     -19.956562                                                       
!    34   15     -24.557549                                                       
!    34   16     -29.931850                                                       
!    34   17     -24.440566                                                       
!    34   18     -18.378264                                                       
!    34   19      -1.481000                                                       
!    34   20      13.153000                                                       
!    35   11      41.153000                                                       
!    35   12      17.390059                                                       
!    35   13       -.058079                                                       
!    35   14     -14.359762                                                       
!    35   15     -24.857613                                                       
!    35   16     -28.846371                                                       
!    35   17     -29.013512                                                       
!    35   18     -23.048208                                                       
!    35   19     -11.167107                                                       
!    35   20       4.439000                                                       
!    36   12      20.912000                                                       
!    36   13       5.916384                                                       
!    36   14     -12.400639                                                       
!    36   15     -20.250840                                                       
!    36   16     -30.663958                                                       
!    36   17     -29.521886                                                       
!    36   18     -30.230439                                                       
!    36   19     -17.425078                                                       
!    36   20      -6.439199                                                       
!    36   21      13.898000                                                       
!    37   12      29.100000                                                       
!    37   13       9.603702                                                       
!    37   14      -6.524192                                                       
!    37   15     -18.994708                                                       
!    37   16     -26.896218                                                       
!    37   17     -31.761519                                                       
!    37   18     -30.948034                                                       
!    37   19     -24.799240                                                       
!    37   20     -13.160607                                                       
!    37   21       2.841000                                                       
!    38   13      15.742000                                                       
!    38   14      -3.744605                                                       
!    38   15     -14.466100                                                       
!    38   16     -26.861076                                                       
!    38   17     -29.797976                                                       
!    38   18     -34.714763                                                       
!    38   19     -28.801691                                                       
!    38   20     -22.059044                                                       
!    38   21      -4.937000                                                       
!    38   22       9.101000                                                       
!    39   13      20.400000                                                       
!    39   14       2.142000                                                       
!    39   15     -12.649687                                                       
!    39   16     -23.161339                                                       
!    39   17     -29.800652                                                       
!    39   18     -33.241836                                                       
!    39   19     -33.806836                                                       
!    39   20     -27.276256                                                       
!    39   21     -14.168014                                                       
!    39   22       1.232000                                                       
!    40   14       5.403000                                                       
!    40   15      -8.336870                                                       
!    40   16     -22.849544                                                       
!    40   17     -27.557728                                                       
!    40   18     -35.039890                                                       
!    40   19     -33.535016                                                       
!    40   20     -34.846110                                                       
!    40   21     -20.526386                                                       
!    40   22      -8.850210                                                       
!    40   23      10.330000                                                       
!    41   14      11.830000                                                       
!    41   15      -4.843768                                                       
!    41   16     -18.601932                                                       
!    41   17     -27.339148                                                       
!    41   18     -33.067261                                                       
!    41   19     -35.558872                                                       
!    41   20     -35.137487                                                       
!    41   21     -28.642211                                                       
!    41   22     -15.713000                                                       
!    41   23       -.242000                                                       
!    42   14      14.997000                                                       
!    42   15        .084000                                                       
!    42   16     -17.241951                                                       
!    42   17     -24.987328                                                       
!    42   18     -34.422065                                                       
!    42   19     -35.021319                                                       
!    42   20     -38.546765                                                       
!    42   21     -32.120925                                                       
!    42   22     -25.120878                                                       
!    42   23      -8.169000                                                       
!    42   24       5.990000                                                       
!    43   15       3.083000                                                       
!    43   16     -12.482018                                                       
!    43   17     -24.029389                                                       
!    43   18     -31.977531                                                       
!    43   19     -36.593042                                                       
!    43   20     -38.408442                                                       
!    43   21     -36.187623                                                       
!    43   22     -29.320314                                                       
!    43   23     -18.024000                                                       
!    43   24      -2.136000                                                       
!    44   15       9.203000                                                       
!    44   16     -10.880000                                                       
!    44   17     -19.991058                                                       
!    44   18     -32.262039                                                       
!    44   19     -35.810214                                                       
!    44   20     -41.469087                                                       
!    44   21     -37.815811                                                       
!    44   22     -37.548299                                                       
!    44   23     -23.846000                                                       
!    44   24     -13.535000                                                       
!    44   25       6.399000                                                       
!    45   15      14.103000                                                       
!    45   16      -4.825000                                                       
!    45   17     -18.909325                                                       
!    45   18     -29.719331                                                       
!    45   19     -36.608028                                                       
!    45   20     -40.812530                                                       
!    45   21     -41.069338                                                       
!    45   22     -39.006912                                                       
!    45   23     -31.873591                                                       
!    45   24     -19.412000                                                       
!    45   25      -5.114000                                                       
!    45   26      13.563000                                                       
!    46   15      22.197000                                                       
!    46   16       -.401000                                                       
!    46   17     -14.792000                                                       
!    46   18     -29.720739                                                       
!    46   19     -35.418933                                                       
!    46   20     -43.134911                                                       
!    46   21     -41.758637                                                       
!    46   22     -44.125338                                                       
!    46   23     -37.073929                                                       
!    46   24     -29.470929                                                       
!    46   25     -12.370000                                                       
!    46   26        .755000                                                       
!    47   16       7.098000                                                       
!    47   17     -11.225000                                                       
!    47   18     -25.908348                                                       
!    47   19     -35.696887                                                       
!    47   20     -42.339694                                                       
!    47   21     -44.331631                                                       
!    47   22     -44.931732                                                       
!    47   23     -42.003929                                                       
!    47   24     -34.552357                                                       
!    47   25     -22.263000                                                       
!    47   26      -6.623000                                                       
!    48   16      12.100000                                                       
!    48   17      -4.797000                                                       
!    48   18     -23.222000                                                       
!    48   19     -32.124479                                                       
!    48   20     -44.214742                                                       
!    48   21     -44.492812                                                       
!    48   22     -48.487000                                                       
!    48   23     -44.474659                                                       
!    48   24     -42.815313                                                       
!    48   25     -28.997000                                                       
!    48   26     -18.108000                                                       
!    48   27       1.639000                                                       
!    49   16      20.502000                                                       
!    49   17       -.102000                                                       
!    49   18     -16.599000                                                       
!    49   19     -30.320047                                                       
!    49   20     -41.290047                                                       
!    49   21     -46.552277                                                       
!    49   22     -48.558040                                                       
!    49   23     -47.956179                                                       
!    49   24     -45.325434                                                       
!    49   25     -37.610541                                                       
!    49   26     -24.582000                                                       
!    49   27      -9.576000                                                       
!    50   17       7.200000                                                       
!    50   18     -13.097000                                                       
!    50   19     -25.352625                                                       
!    50   20     -39.571456                                                       
!    50   21     -44.537508                                                       
!    50   22     -51.425849                                                       
!    50   23     -49.217535                                                       
!    50   24     -50.254460                                                       
!    50   25     -42.621470                                                       
!    50   26     -34.471499                                                       
!    50   27     -17.195000                                                       
!    50   28      -3.791000                                                       
!    51   17      12.603000                                                       
!    51   18      -6.297000                                                       
!    51   19     -22.002000                                                       
!    51   20     -35.886511                                                       
!    51   21     -43.218800                                                       
!    51   22     -49.726852                                                       
!    51   23     -52.197493                                                       
!    51   24     -51.444760                                                       
!    51   25     -48.236956                                                       
!    51   26     -40.217307                                                       
!    51   27     -27.274000                                                       
!    51   28     -11.439000                                                       
!    52   18      -1.705000                                                       
!    52   19     -16.199000                                                       
!    52   20     -32.509136                                                       
!    52   21     -40.380259                                                       
!    52   22     -49.464024                                                       
!    52   23     -51.437409                                                       
!    52   24     -55.412797                                                       
!    52   25     -50.701137                                                       
!    52   26     -48.329137                                                       
!    52   27     -33.916000                                                       
!    52   28     -22.654000                                                       
!    52   29      -2.627000                                                       
!    53   18       5.800000                                                       
!    53   19     -11.998000                                                       
!    53   20     -27.898000                                                       
!    53   21     -38.964388                                                       
!    53   22     -46.824604                                                       
!    53   23     -51.844604                                                       
!    53   24     -55.280638                                                       
!    53   25     -54.683625                                                       
!    53   26     -50.941273                                                       
!    53   27     -42.639145                                                       
!    53   28     -29.379000                                                       
!    53   29     -13.460000                                                       
!    54   19      -5.598000                                                       
!    54   20     -23.585000                                                       
!    54   21     -34.465273                                                       
!    54   22     -45.764293                                                       
!    54   23     -49.886730                                                       
!    54   24     -56.928321                                                       
!    54   25     -55.551269                                                       
!    54   26     -56.248410                                                       
!    54   27     -48.005330                                                       
!    54   28     -39.206103                                                       
!    54   29     -21.694000                                                       
!    54   30      -6.567000                                                       
!    55   19       -.570000                                                       
!    55   20     -18.118000                                                       
!    55   21     -28.503712                                                       
!    55   22     -41.805444                                                       
!    55   23     -49.147298                                                       
!    55   24     -55.103298                                                       
!    55   25     -57.706384                                                       
!    55   26     -57.475007                                                       
!    55   27     -54.023711                                                       
!    55   28     -45.329911                                                       
!    55   29     -31.624000                                                       
!    55   30     -14.923000                                                       
!    56   20     -13.237000                                                       
!    56   21     -25.467000                                                       
!    56   22     -39.132057                                                       
!    56   23     -46.239355                                                       
!    56   24     -55.288596                                                       
!    56   25     -56.905551                                                       
!    56   26     -60.601003                                                       
!    56   27     -56.035003                                                       
!    56   28     -53.899645                                                       
!    56   29     -38.601000                                                       
!    56   30     -25.728000                                                       
!    56   31      -4.741000                                                       
!    57   20      -7.120000                                                       
!    57   21     -21.387000                                                       
!    57   22     -33.254331                                                       
!    57   23     -44.376367                                                       
!    57   24     -52.392990                                                       
!    57   25     -57.484854                                                       
!    57   26     -60.175708                                                       
!    57   27     -59.339666                                                       
!    57   28     -56.075475                                                       
!    57   29     -47.305268                                                       
!    57   30     -32.686000                                                       
!    57   31     -15.901000                                                       
!    58   21     -15.770000                                                       
!    58   22     -31.568000                                                       
!    58   23     -40.380259                                                       
!    58   24     -51.930783                                                       
!    58   25     -55.902253                                                       
!    58   26     -62.148844                                                       
!    58   27     -59.841428                                                       
!    58   28     -60.223014                                                       
!    58   29     -51.659966                                                       
!    58   30     -42.293114                                                       
!    58   31     -23.986000                                                       
!    58   32      -8.374000                                                       
!    59   21     -11.140000                                                       
!    59   22     -26.119000                                                       
!    59   23     -37.911800                                                       
!    59   24     -47.850840                                                       
!    59   25     -55.473100                                                       
!    59   26     -60.658421                                                       
!    59   27     -62.223609                                                       
!    59   28     -61.151125                                                       
!    59   29     -56.351547                                                       
!    59   30     -47.257409                                                       
!    59   31     -34.121000                                                       
!    59   32     -17.000000                                                       
!    60   22     -22.691000                                                       
!    60   23     -33.068032                                                       
!    60   24     -46.826196                                                       
!    60   25     -52.914442                                                       
!    60   26     -61.406923                                                       
!    60   27     -61.644218                                                       
!    60   28     -64.468100                                                       
!    60   29     -58.341209                                                       
!    60   30     -54.183106                                                       
!    60   31     -39.998000                                                       
!    60   32     -27.768000                                                       
!    60   33      -6.399000                                                       
!    61   22     -16.750000                                                       
!    61   23     -30.357000                                                       
!    61   24     -42.764883                                                       
!    61   25     -51.735169                                                       
!    61   26     -58.917489                                                       
!    61   27     -62.895042                                                       
!    61   28     -64.216775                                                       
!    61   29     -61.979570                                                       
!    61   30     -56.342425                                                       
!    61   31     -47.348000                                                       
!    61   32     -33.729000                                                       
!    61   33     -18.052000                                                       
!    62   23     -25.020000                                                       
!    62   24     -41.172029                                                       
!    62   25     -48.465626                                                       
!    62   26     -58.897896                                                       
!    62   27     -61.428097                                                       
!    62   28     -66.742688                                                       
!    62   29     -62.794517                                                       
!    62   30     -61.167353                                                       
!    62   31     -51.996353                                                       
!    62   32     -42.243000                                                       
!    62   33     -24.964000                                                       
!    63   23     -21.657000                                                       
!    63   24     -35.527000                                                       
!    63   25     -46.751677                                                       
!    63   26     -55.779304                                                       
!    63   27     -61.837017                                                       
!    63   28     -65.509217                                                       
!    63   29     -65.576162                                                       
!    63   30     -62.209293                                                       
!    63   31     -56.689293                                                       
!    63   32     -46.910000                                                       
!    63   33     -33.823000                                                       
!    64   24     -33.347000                                                       
!    64   25     -43.100221                                                       
!    64   26     -55.079232                                                       
!    64   27     -59.789309                                                       
!    64   28     -67.095900                                                       
!    64   29     -65.420802                                                       
!    64   30     -65.999528                                                       
!    64   31     -58.834729                                                       
!    64   32     -54.424729                                                       
!    64   33     -39.521000                                                       
!    65   24     -27.600000                                                       
!    65   25     -40.892580                                                       
!    65   26     -51.288052                                                       
!    65   27     -59.164222                                                       
!    65   28     -65.122587                                                       
!    65   29     -67.259719                                                       
!    65   30     -65.907775                                                       
!    65   31     -62.652909                                                       
!    65   32     -56.410559                                                       
!    65   33     -47.056000                                                       
!    65   34     -32.919000                                                       
!    66   25     -36.496000                                                       
!    66   26     -50.319298                                                       
!    66   27     -56.052259                                                       
!    66   28     -66.028726                                                       
!    66   29     -66.254326                                                       
!    66   30     -68.896301                                                       
!    66   31     -63.721301                                                       
!    66   32     -61.621301                                                       
!    66   33     -52.071301                                                       
!    66   34     -41.722000                                                       
!    67   25     -33.701000                                                       
!    67   26     -46.574693                                                       
!    67   27     -55.321420                                                       
!    67   28     -63.742462                                                       
!    67   29     -67.300158                                                       
!    67   30     -67.877158                                                       
!    67   31     -66.876681                                                       
!    67   32     -62.653753                                                       
!    67   33     -56.643753                                                       
!    67   34     -46.491000                                                       
!    67   35     -32.798000                                                       
!    68   26     -44.237000                                                       
!    68   27     -51.828318                                                       
!    68   28     -63.486027                                                       
!    68   29     -65.541887                                                       
!    68   30     -70.004030                                                       
!    68   31     -67.082930                                                       
!    68   32     -66.976955                                                       
!    68   33     -58.876955                                                       
!    68   34     -54.148000                                                       
!    68   35     -38.892000                                                       
!    69   26     -39.402000                                                       
!    69   27     -51.045864                                                       
!    69   28     -60.377721                                                       
!    69   29     -65.739917                                                       
!    69   30     -68.414928                                                       
!    69   31     -69.320923                                                       
!    69   32     -67.093638                                                       
!    69   33     -63.080621                                                       
!    69   34     -56.297481                                                       
!    69   35     -46.409000                                                       
!    69   36     -32.304000                                                       
!    70   27     -46.752000                                                       
!    70   28     -59.485198                                                       
!    70   29     -62.960334                                                       
!    70   30     -69.559425                                                       
!    70   31     -68.904705                                                       
!    70   32     -70.560320                                                       
!    70   33     -64.340320                                                       
!    70   34     -61.940000                                                       
!    70   35     -51.970000                                                       
!    70   36     -40.976000                                                       
!    71   27     -44.963000                                                       
!    71   28     -55.889632                                                       
!    71   29     -62.764225                                                       
!    71   30     -67.321674                                                       
!    71   31     -70.136821                                                       
!    71   32     -69.904897                                                       
!    71   33     -67.892187                                                       
!    71   34     -63.464187                                                       
!    71   35     -56.592000                                                       
!    71   36     -46.100000                                                       
!    71   37     -32.304000                                                       
!    72   27     -40.604000                                                       
!    72   28     -54.678690                                                       
!    72   29     -60.063000                                                       
!    72   30     -68.128416                                                       
!    72   31     -68.586498                                                       
!    72   32     -72.585557                                                       
!    72   33     -68.229459                                                       
!    72   34     -67.894432                                                       
!    72   35     -59.152770                                                       
!    72   36     -54.112770                                                       
!    72   37     -38.117000                                                       
!    73   28     -50.226000                                                       
!    73   29     -59.159000                                                       
!    73   30     -65.409994                                                       
!    73   31     -69.703842                                                       
!    73   32     -71.297136                                                       
!    73   33     -70.956276                                                       
!    73   34     -68.216276                                                       
!    73   35     -63.532642                                                       
!    73   36     -56.885291                                                       
!    73   37     -46.234000                                                       
!    73   38     -31.699000                                                       
!    74   28     -48.522000                                                       
!    74   29     -55.703000                                                       
!    74   30     -65.709197                                                       
!    74   31     -68.054011                                                       
!    74   32     -73.422011                                                       
!    74   33     -70.859599                                                       
!    74   34     -72.212607                                                       
!    74   35     -65.305961                                                       
!    74   36     -62.169554                                                       
!    74   37     -51.725504                                                       
!    74   38     -40.697000                                                       
!    75   28     -43.808000                                                       
!    75   29     -54.306000                                                       
!    75   30     -62.468419                                                       
!    75   31     -68.464198                                                       
!    75   32     -71.855908                                                       
!    75   33     -73.032457                                                       
!    75   34     -72.168818                                                       
!    75   35     -69.138818                                                       
!    75   36     -64.241597                                                       
!    75   37     -57.222414                                                       
!    75   38     -46.649000                                                       
!    76   28     -41.610000                                                       
!    76   29     -50.310000                                                       
!    76   30     -62.042888                                                       
!    76   31     -66.202888                                                       
!    76   32     -73.212888                                                       
!    76   33     -72.289575                                                       
!    76   34     -75.251563                                                       
!    76   35     -70.288688                                                       
!    76   36     -68.978700                                                       
!    76   37     -60.480547                                                       
!    76   38     -54.390000                                                       
!    77   28     -36.487000                                                       
!    77   29     -48.484000                                                       
!    77   30     -58.604138                                                       
!    77   31     -65.874138                                                       
!    77   32     -71.214138                                                       
!    77   33     -73.916177                                                       
!    77   34     -74.599049                                                       
!    77   35     -73.233933                                                       
!    77   36     -70.171407                                                       
!    77   37     -64.825826                                                       
!    77   38     -57.974770                                                       
!    77   39     -46.929000                                                       
!    78   28     -33.720000                                                       
!    78   29     -43.957000                                                       
!    78   30     -57.222063                                                       
!    78   31     -63.662063                                                       
!    78   32     -71.862063                                                       
!    78   33     -72.816201                                                       
!    78   34     -77.025673                                                       
!    78   35     -73.451896                                                       
!    78   36     -74.159700                                                       
!    78   37     -66.935766                                                       
!    78   38     -63.174508                                                       
!    78   39     -52.629000                                                       
!    79   29     -41.656000                                                       
!    79   30     -53.937989                                                       
!    79   31     -62.487989                                                       
!    79   32     -69.487989                                                       
!    79   33     -73.635989                                                       
!    79   34     -75.916934                                                       
!    79   35     -76.067980                                                       
!    79   36     -74.442202                                                       
!    79   37     -70.796589                                                       
!    79   38     -65.477427                                                       
!    79   39     -58.357427                                                       
!    79   40     -47.357000                                                       
!    80   29     -35.499000                                                       
!    80   30     -51.777345                                                       
!    80   31     -59.067745                                                       
!    80   32     -69.447745                                                       
!    80   33     -72.117967                                                       
!    80   34     -77.759404                                                       
!    80   35     -75.888849                                                       
!    80   36     -77.893342                                                       
!    80   37     -72.172776                                                       
!    80   38     -70.304883                                                       
!    80   39     -63.357975                                                       
!    80   40     -55.377000                                                       
!    81   30     -46.128000                                                       
!    81   31     -57.982741                                                       
!    81   32     -66.302741                                                       
!    81   33     -72.532741                                                       
!    81   34     -76.389081                                                       
!    81   35     -77.974363                                                       
!    81   36     -77.693650                                                       
!    81   37     -75.456438                                                       
!    81   38     -71.526530                                                       
!    81   39     -66.016164                                                       
!    81   40     -58.856164                                                       
!    81   41     -47.460000                                                       
!    82   30     -42.066000                                                       
!    82   31     -52.946000                                                       
!    82   32     -65.623439                                                       
!    82   33     -70.323439                                                       
!    82   34     -77.593439                                                       
!    82   35     -77.495943                                                       
!    82   36     -80.588563                                                       
!    82   37     -76.189034                                                       
!    82   38     -76.008727                                                       
!    82   39     -68.192736                                                       
!    82   40     -64.192736                                                       
!    82   41     -52.974000                                                       
!    83   31     -49.490000                                                       
!    83   32     -61.004000                                                       
!    83   33     -69.880088                                                       
!    83   34     -75.340088                                                       
!    83   35     -79.009106                                                       
!    83   36     -79.981833                                                       
!    83   37     -79.072696                                                       
!    83   38     -76.796985                                                       
!    83   39     -72.328103                                                       
!    83   40     -66.460103                                                       
!    83   41     -58.960103                                                       
!    83   42     -47.748000                                                       
!    84   31     -44.395000                                                       
!    84   32     -58.395000                                                       
!    84   33     -66.080000                                                       
!    84   34     -75.949796                                                       
!    84   35     -77.776304                                                       
!    84   36     -82.431033                                                       
!    84   37     -79.750149                                                       
!    84   38     -80.644288                                                       
!    84   39     -74.158306                                                       
!    84   40     -71.492000                                                       
!    84   41     -61.879000                                                       
!    84   42     -55.806000                                                       
!    85   32     -53.384000                                                       
!    85   33     -63.519000                                                       
!    85   34     -72.428605                                                       
!    85   35     -78.610605                                                       
!    85   36     -81.480605                                                       
!    85   37     -82.167687                                                       
!    85   38     -81.102665                                                       
!    85   39     -77.847665                                                       
!    85   40     -73.154665                                                       
!    85   41     -67.154665                                                       
!    85   42     -59.066000                                                       
!    85   43     -47.562000                                                       
!    86   32     -50.049000                                                       
!    86   33     -59.401000                                                       
!    86   34     -70.540945                                                       
!    86   35     -75.639945                                                       
!    86   36     -83.265945                                                       
!    86   37     -82.747319                                                       
!    86   38     -84.521563                                                       
!    86   39     -79.281563                                                       
!    86   40     -77.805027                                                       
!    86   41     -69.827027                                                       
!    86   42     -64.557027                                                       
!    86   43     -53.207000                                                       
!    87   33     -56.281000                                                       
!    87   34     -66.582484                                                       
!    87   35     -73.857484                                                       
!    87   36     -80.709984                                                       
!    87   37     -84.595045                                                       
!    87   38     -84.878358                                                       
!    87   39     -83.016752                                                       
!    87   40     -79.347835                                                       
!    87   41     -74.182835                                                       
!    87   42     -67.694611                                                       
!    87   43     -59.122000                                                       
!    87   44     -47.339000                                                       
!    88   33     -51.642000                                                       
!    88   34     -63.878140                                                       
!    88   35     -70.732140                                                       
!    88   36     -79.692140                                                       
!    88   37     -82.606221                                                       
!    88   38     -87.919663                                                       
!    88   39     -84.297063                                                       
!    88   40     -83.623763                                                       
!    88   41     -76.073763                                                       
!    88   42     -72.700555                                                       
!    88   43     -62.568000                                                       
!    88   44     -55.498000                                                       
!    89   33     -47.292000                                                       
!    89   34     -59.597000                                                       
!    89   35     -68.569815                                                       
!    89   36     -76.724815                                                       
!    89   37     -81.710698                                                       
!    89   38     -86.207050                                                       
!    89   39     -87.702101                                                       
!    89   40     -84.869416                                                       
!    89   41     -80.578408                                                       
!    89   42     -75.003363                                                       
!    89   43     -67.493363                                                       
!    89   44     -59.513000                                                       
!    89   45     -47.152000                                                       
!    90   34     -56.430000                                                       
!    90   35     -64.613083                                                       
!    90   36     -74.963083                                                       
!    90   37     -79.354948                                                       
!    90   38     -85.941863                                                       
!    90   39     -86.487861                                                       
!    90   40     -88.767938                                                       
!    90   41     -82.656938                                                       
!    90   42     -80.167938                                                       
!    90   43     -71.207275                                                       
!    90   44     -65.409000                                                       
!    90   45     -53.216000                                                       
!    91   34     -50.888000                                                       
!    91   35     -61.510918                                                       
!    91   36     -71.312918                                                       
!    91   37     -77.747918                                                       
!    91   38     -83.638978                                                       
!    91   39     -86.346300                                                       
!    91   40     -87.891133                                                       
!    91   41     -86.637743                                                       
!    91   42     -82.203631                                                       
!    91   43     -75.983631                                                       
!    91   44     -68.578969                                                       
!    91   45     -59.103000                                                       
!    91   46     -47.059000                                                       
!    92   34     -47.199000                                                       
!    92   35     -56.583354                                                       
!    92   36     -68.788258                                                       
!    92   37     -74.775258                                                       
!    92   38     -82.875106                                                       
!    92   39     -84.815467                                                       
!    92   40     -88.454559                                                       
!    92   41     -86.448952                                                       
!    92   42     -86.805466                                                       
!    92   43     -78.935111                                                       
!    92   44     -74.408000                                                       
!    92   45     -63.360000                                                       
!    92   46     -55.498000                                                       
!    93   35     -53.002000                                                       
!    93   36     -64.026001                                                       
!    93   37     -72.626001                                                       
!    93   38     -80.087598                                                       
!    93   39     -84.224201                                                       
!    93   40     -87.117379                                                       
!    93   41     -87.208744                                                       
!    93   42     -86.803852                                                       
!    93   43     -83.602997                                                       
!    93   44     -77.265997                                                       
!    93   45     -69.173000                                                       
!    93   46     -59.699000                                                       
!    94   35     -47.804000                                                       
!    94   36     -61.141000                                                       
!    94   37     -68.551124                                                       
!    94   38     -78.841774                                                       
!    94   39     -82.349639                                                       
!    94   40     -87.266289                                                       
!    94   41     -86.364891                                                       
!    94   42     -88.410338                                                       
!    94   43     -84.154593                                                       
!    94   44     -82.568017                                                       
!    94   45     -72.938000                                                       
!    94   46     -66.350000                                                       
!    94   47     -53.300000                                                       
!    95   36     -56.039000                                                       
!    95   37     -65.838676                                                       
!    95   38     -75.117330                                                       
!    95   39     -81.204179                                                       
!    95   40     -85.657624                                                       
!    95   41     -86.782460                                                       
!    95   42     -87.708077                                                       
!    95   43     -86.017446                                                       
!    95   44     -83.449993                                                       
!    95   45     -78.339993                                                       
!    95   46     -70.151000                                                       
!    95   47     -60.100000                                                       
!    96   36     -53.030000                                                       
!    96   37     -61.214086                                                       
!    96   38     -72.954158                                                       
!    96   39     -78.340695                                                       
!    96   40     -85.440646                                                       
!    96   41     -85.604215                                                       
!    96   42     -88.791015                                                       
!    96   43     -85.817781                                                       
!    96   44     -86.072192                                                       
!    96   45     -79.625760                                                       
!    96   46     -76.175760                                                       
!    96   47     -64.571000                                                       
!    96   48     -56.104000                                                       
!    97   36     -47.916000                                                       
!    97   37     -58.364738                                                       
!    97   38     -68.791979                                                       
!    97   39     -76.260455                                                       
!    97   40     -82.948861                                                       
!    97   41     -85.606946                                                       
!    97   42     -87.540831                                                       
!    97   43     -87.220574                                                       
!    97   44     -86.112372                                                       
!    97   45     -82.589372                                                       
!    97   46     -77.799372                                                       
!    97   47     -70.794000                                                       
!    97   48     -60.603000                                                       
!    98   37     -54.302779                                                       
!    98   38     -66.628659                                                       
!    98   39     -72.452035                                                       
!    98   40     -81.276224                                                       
!    98   41     -83.526412                                                       
!    98   42     -88.112011                                                       
!    98   43     -86.428013                                                       
!    98   44     -88.224474                                                       
!    98   45     -83.167096                                                       
!    98   46     -81.300084                                                       
!    98   47     -72.880084                                                       
!    98   48     -67.460000                                                       
!    98   49     -53.803000                                                       
!    99   37     -50.840361                                                       
!    99   38     -62.116633                                                       
!    99   39     -70.202280                                                       
!    99   40     -77.769413                                                       
!    99   41     -82.327417                                                       
!    99   42     -85.966080                                                       
!    99   43     -87.323307                                                       
!    99   44     -87.616958                                                       
!    99   45     -85.574384                                                       
!    99   46     -82.187793                                                       
!    99   47     -76.757793                                                       
!    99   48     -69.853000                                                       
!    99   49     -60.910000                                                       
!   100   37     -46.696000                                                       
!   100   38     -60.219468                                                       
!   100   39     -67.294468                                                       
!   100   40     -76.604468                                                       
!   100   41     -79.939468                                                       
!   100   42     -86.184468                                                       
!   100   43     -86.016384                                                       
!   100   44     -89.218795                                                       
!   100   45     -85.588795                                                       
!   100   46     -85.227407                                                       
!   100   47     -78.180851                                                       
!   100   48     -74.305049                                                       
!   100   49     -64.134253                                                       
!   100   50     -56.864000                                                       
!   101   37     -43.597645                                                       
!   101   38     -55.407645                                                       
!   101   39     -64.912645                                                       
!   101   40     -73.457645                                                       
!   101   41     -78.942645                                                       
!   101   42     -83.511645                                                       
!   101   43     -86.336086                                                       
!   101   44     -87.949584                                                       
!   101   45     -87.408099                                                       
!   101   46     -85.428099                                                       
!   101   47     -81.224275                                                       
!   101   48     -75.747738                                                       
!   101   49     -68.409000                                                       
!   101   50     -59.560000                                                       
!   102   37     -37.996000                                                       
!   102   38     -53.077643                                                       
!   102   39     -61.892643                                                       
!   102   40     -71.742643                                                       
!   102   41     -76.347643                                                       
!   102   42     -83.557643                                                       
!   102   43     -84.567591                                                       
!   102   44     -89.097851                                                       
!   102   45     -86.775318                                                       
!   102   46     -87.925833                                                       
!   102   47     -81.971463                                                       
!   102   48     -79.384463                                                       
!   102   49     -70.134463                                                       
!   102   50     -64.748000                                                       
!   103   38     -47.553000                                                       
!   103   39     -58.740000                                                       
!   103   40     -68.374386                                                       
!   103   41     -75.319386                                                       
!   103   42     -80.849386                                                       
!   103   43     -84.599386                                                       
!   103   44     -87.258919                                                       
!   103   45     -88.022274                                                       
!   103   46     -87.479193                                                       
!   103   47     -84.791601                                                       
!   103   48     -80.649714                                                       
!   103   49     -74.599714                                                       
!   103   50     -66.946000                                                       
!   103   51     -55.778000                                                       
!   104   38     -44.404000                                                       
!   104   39     -54.539000                                                       
!   104   40     -66.341000                                                       
!   104   41     -72.228533                                                       
!   104   42     -80.333533                                                       
!   104   43     -82.488533                                                       
!   104   44     -88.091239                                                       
!   104   45     -86.950001                                                       
!   104   46     -89.390890                                                       
!   104   47     -85.112244                                                       
!   104   48     -83.975950                                                       
!   104   49     -76.067264                                                       
!   104   50     -71.552264                                                       
!   104   51     -59.348000                                                       
!   105   39     -51.148000                                                       
!   105   40     -62.364000                                                       
!   105   41     -70.854991                                                       
!   105   42     -77.339991                                                       
!   105   43     -82.289991                                                       
!   105   44     -85.929991                                                       
!   105   45     -87.846909                                                       
!   105   46     -88.413628                                                       
!   105   47     -87.068376                                                       
!   105   48     -84.330172                                                       
!   105   49     -79.481172                                                       
!   105   50     -73.224351                                                       
!   105   51     -63.780695                                                       
!   106   39     -46.370000                                                       
!   106   40     -59.699000                                                       
!   106   41     -66.891000                                                       
!   106   42     -76.257412                                                       
!   106   43     -79.777412                                                       
!   106   44     -86.324412                                                       
!   106   45     -86.363812                                                       
!   106   46     -89.904912                                                       
!   106   47     -86.939647                                                       
!   106   48     -87.133792                                                       
!   106   49     -80.610422                                                       
!   106   50     -77.425330                                                       
!   106   51     -66.896970                                                       
!   106   52     -58.030000                                                       
!   107   40     -55.089000                                                       
!   107   41     -64.916000                                                       
!   107   42     -72.940884                                                       
!   107   43     -79.100884                                                       
!   107   44     -83.920884                                                       
!   107   45     -86.861300                                                       
!   107   46     -88.372264                                                       
!   107   47     -88.405269                                                       
!   107   48     -86.988269                                                       
!   107   49     -83.562269                                                       
!   107   50     -78.555949                                                       
!   107   51     -70.654000                                                       
!   107   52     -60.513000                                                       
!   108   40     -51.903000                                                       
!   108   41     -60.538000                                                       
!   108   42     -70.817904                                                       
!   108   43     -75.935404                                                       
!   108   44     -83.655404                                                       
!   108   45     -85.016729                                                       
!   108   46     -89.521729                                                       
!   108   47     -87.603546                                                       
!   108   48     -89.252572                                                       
!   108   49     -84.095561                                                       
!   108   50     -82.003746                                                       
!   108   51     -72.507000                                                       
!   108   52     -65.682577                                                       
!   108   53     -52.824000                                                       
!   109   41     -58.097000                                                       
!   109   42     -67.245000                                                       
!   109   43     -74.537209                                                       
!   109   44     -80.852209                                                       
!   109   45     -85.012209                                                       
!   109   46     -87.603705                                                       
!   109   47     -88.719654                                                       
!   109   48     -88.505359                                                       
!   109   49     -86.485406                                                       
!   109   50     -82.635728                                                       
!   109   51     -76.255728                                                       
!   109   52     -67.573840                                                       
!   109   53     -57.574091                                                       
!   110   41     -53.393000                                                       
!   110   42     -65.456000                                                       
!   110   43     -71.362000                                                       
!   110   44     -80.139971                                                       
!   110   45     -82.949971                                                       
!   110   46     -88.349971                                                       
!   110   47     -87.457530                                                       
!   110   48     -90.349708                                                       
!   110   49     -86.471708                                                       
!   110   50     -85.834656                                                       
!   110   51     -76.816656                                                       
!   110   52     -72.277250                                                       
!   110   53     -60.887859                                                       
!   110   54     -51.721000                                                       
!   111   42     -61.004000                                                       
!   111   43     -69.815000                                                       
!   111   44     -76.792000                                                       
!   111   45     -82.288000                                                       
!   111   46     -86.029092                                                       
!   111   47     -88.217425                                                       
!   111   48     -89.254225                                                       
!   111   49     -88.388821                                                       
!   111   50     -85.943905                                                       
!   111   51     -81.473905                                                       
!   111   52     -73.475687                                                       
!   111   53     -64.947000                                                       
!   111   54     -54.369000                                                       
!   112   42     -58.833000                                                       
!   112   43     -65.913000                                                       
!   112   44     -75.617113                                                       
!   112   45     -80.137113                                                       
!   112   46     -86.337113                                                       
!   112   47     -86.625080                                                       
!   112   48     -90.581047                                                       
!   112   49     -87.995116                                                       
!   112   50     -88.658831                                                       
!   112   51     -81.603855                                                       
!   112   52     -77.256594                                                       
!   112   53     -67.096000                                                       
!   112   54     -59.927331                                                       
!   112   55     -46.266000                                                       
!   113   42     -53.999000                                                       
!   113   43     -63.966000                                                       
!   113   44     -72.154000                                                       
!   113   45     -78.786000                                                       
!   113   46     -83.693469                                                       
!   113   47     -87.033469                                                       
!   113   48     -89.049931                                                       
!   113   49     -89.366381                                                       
!   113   50     -88.330421                                                       
!   113   51     -84.413891                                                       
!   113   52     -78.755430                                                       
!   113   53     -71.124917                                                       
!   113   54     -62.053479                                                       
!   113   55     -51.664830                                                       
!   114   43     -59.727000                                                       
!   114   44     -70.894153                                                       
!   114   45     -76.994153                                                       
!   114   46     -83.494153                                                       
!   114   47     -84.944875                                                       
!   114   48     -90.021317                                                       
!   114   49     -88.569457                                                       
!   114   50     -90.558142                                                       
!   114   51     -84.676634                                                       
!   114   52     -81.507552                                                       
!   114   53     -72.796000                                                       
!   114   54     -66.932000                                                       
!   114   55     -55.105948                                                       
!   114   56     -45.698000                                                       
!   115   43     -57.492000                                                       
!   115   44     -66.779000                                                       
!   115   45     -74.403374                                                       
!   115   46     -80.403374                                                       
!   115   47     -84.987374                                                       
!   115   48     -88.090859                                                       
!   115   49     -89.536747                                                       
!   115   50     -90.032633                                                       
!   115   51     -87.002633                                                       
!   115   52     -82.363966                                                       
!   115   53     -76.129000                                                       
!   115   54     -68.018583                                                       
!   115   55     -59.672000                                                       
!   115   56     -48.708000                                                       
!   116   44     -65.056000                                                       
!   116   45     -71.961027                                                       
!   116   46     -79.961027                                                       
!   116   47     -82.568027                                                       
!   116   48     -88.719729                                                       
!   116   49     -88.249733                                                       
!   116   50     -91.524722                                                       
!   116   51     -86.817803                                                       
!   116   52     -85.305972                                                       
!   116   53     -77.560823                                                       
!   116   54     -73.220823                                                       
!   116   55     -62.490055                                                       
!   116   56     -54.325000                                                       
!   117   44     -60.743000                                                       
!   117   45     -69.536000                                                       
!   117   46     -76.532000                                                       
!   117   47     -82.265638                                                       
!   117   48     -86.425638                                                       
!   117   49     -88.943010                                                       
!   117   50     -90.397972                                                       
!   117   51     -88.641338                                                       
!   117   52     -85.106700                                                       
!   117   53     -80.436645                                                       
!   117   54     -73.993816                                                       
!   117   55     -66.471883                                                       
!   117   56     -58.031854                                                       
!   117   57     -46.565000                                                       
!   118   44     -58.656000                                                       
!   118   45     -65.736000                                                       
!   118   46     -75.465986                                                       
!   118   47     -79.565986                                                       
!   118   48     -86.708904                                                       
!   118   49     -87.230094                                                       
!   118   50     -91.653102                                                       
!   118   51     -87.996470                                                       
!   118   52     -87.723260                                                       
!   118   53     -80.690442                                                       
!   118   54     -77.713676                                                       
!   118   55     -68.413676                                                       
!   118   56     -62.000000                                                       
!   118   57     -49.770000                                                       
!   119   45     -63.938000                                                       
!   119   46     -72.023000                                                       
!   119   47     -78.556562                                                       
!   119   48     -83.906562                                                       
!   119   49     -87.703562                                                       
!   119   50     -90.067185                                                       
!   119   51     -89.473283                                                       
!   119   52     -87.180271                                                       
!   119   53     -83.665999                                                       
!   119   54     -78.660656                                                       
!   119   55     -72.311049                                                       
!   119   56     -64.224707                                                       
!   119   57     -54.967000                                                       
!   119   58     -44.004000                                                       
!   120   45     -59.821000                                                       
!   120   46     -70.766000                                                       
!   120   47     -75.647912                                                       
!   120   48     -83.972912                                                       
!   120   49     -85.733293                                                       
!   120   50     -91.103293                                                       
!   120   51     -88.422693                                                       
!   120   52     -89.404882                                                       
!   120   53     -83.789882                                                       
!   120   54     -81.829882                                                       
!   120   55     -73.887752                                                       
!   120   56     -68.887752                                                       
!   120   57     -57.687000                                                       
!   120   58     -49.705000                                                       
!   121   45     -57.678000                                                       
!   121   46     -66.900000                                                       
!   121   47     -74.658233                                                       
!   121   48     -81.058233                                                       
!   121   49     -85.838233                                                       
!   121   50     -89.202770                                                       
!   121   51     -89.592901                                                       
!   121   52     -88.557294                                                       
!   121   53     -86.287943                                                       
!   121   54     -82.542934                                                       
!   121   55     -77.142934                                                       
!   121   56     -70.340913                                                       
!   121   57     -62.401000                                                       
!   121   58     -52.471000                                                       
!   121   59     -41.579000                                                       
!   122   46     -65.391000                                                       
!   122   47     -71.427000                                                       
!   122   48     -80.574000                                                       
!   122   49     -83.576327                                                       
!   122   50     -89.944918                                                       
!   122   51     -88.328519                                                       
!   122   52     -90.311064                                                       
!   122   53     -86.077064                                                       
!   122   54     -85.186607                                                       
!   122   55     -78.131893                                                       
!   122   56     -74.277000                                                       
!   122   57     -64.543000                                                       
!   122   58     -57.743000                                                       
!   122   59     -45.038000                                                       
!   123   46     -61.236000                                                       
!   123   47     -69.955000                                                       
!   123   48     -77.310567                                                       
!   123   49     -83.425567                                                       
!   123   50     -87.819470                                                       
!   123   51     -89.222491                                                       
!   123   52     -89.169158                                                       
!   123   53     -87.934936                                                       
!   123   54     -85.258936                                                       
!   123   55     -81.049124                                                       
!   123   56     -75.591000                                                       
!   123   57     -68.707000                                                       
!   123   58     -60.072000                                                       
!   123   59     -50.338000                                                       
!   124   47     -66.574000                                                       
!   124   48     -76.710101                                                       
!   124   49     -80.876101                                                       
!   124   50     -88.236101                                                       
!   124   51     -87.618618                                                       
!   124   52     -90.523071                                                       
!   124   53     -87.363484                                                       
!   124   54     -87.657508                                                       
!   124   55     -81.742563                                                       
!   124   56     -79.094600                                                       
!   124   57     -70.300000                                                       
!   124   58     -64.720000                                                       
!   124   59     -53.132000                                                       
!   125   47     -64.702000                                                       
!   125   48     -73.357777                                                       
!   125   49     -80.479777                                                       
!   125   50     -85.897777                                                       
!   125   51     -88.261089                                                       
!   125   52     -89.027789                                                       
!   125   53     -88.842019                                                       
!   125   54     -87.189468                                                       
!   125   55     -84.090728                                                       
!   125   56     -79.530728                                                       
!   125   57     -73.895000                                                       
!   125   58     -66.565000                                                       
!   125   59     -57.911000                                                       
!   126   47     -61.013000                                                       
!   126   48     -72.326776                                                       
!   126   49     -77.812776                                                       
!   126   50     -86.019776                                                       
!   126   51     -86.397776                                                       
!   126   52     -90.070293                                                       
!   126   53     -87.914962                                                       
!   126   54     -89.172962                                                       
!   126   55     -84.348676                                                       
!   126   56     -82.675533                                                       
!   126   57     -75.106000                                                       
!   126   58     -70.700000                                                       
!   126   59     -60.258000                                                       
!   126   60     -53.030000                                                       
!   127   47     -58.796000                                                       
!   127   48     -68.525512                                                       
!   127   49     -76.993512                                                       
!   127   50     -83.507512                                                       
!   127   51     -86.708512                                                       
!   127   52     -88.289512                                                       
!   127   53     -88.987080                                                       
!   127   54     -88.324638                                                       
!   127   55     -86.239937                                                       
!   127   56     -82.789937                                                       
!   127   57     -78.096000                                                       
!   127   58     -71.958000                                                       
!   127   59     -64.431000                                                       
!   127   60     -55.424000                                                       
!   128   48     -67.290542                                                       
!   128   49     -74.360542                                                       
!   128   50     -83.336142                                                       
!   128   51     -84.610075                                                       
!   128   52     -88.993635                                                       
!   128   53     -87.741827                                                       
!   128   54     -89.860807                                                       
!   128   55     -85.932247                                                       
!   128   56     -85.409725                                                       
!   128   57     -78.759725                                                       
!   128   58     -75.572000                                                       
!   128   59     -66.322000                                                       
!   128   60     -60.184000                                                       
!   128   61     -48.195000                                                       
!   129   48     -63.099000                                                       
!   129   49     -72.975131                                                       
!   129   50     -80.630131                                                       
!   129   51     -84.626131                                                       
!   129   52     -87.005631                                                       
!   129   53     -88.503573                                                       
!   129   54     -88.697350                                                       
!   129   55     -87.501395                                                       
!   129   56     -85.069842                                                       
!   129   57     -81.349842                                                       
!   129   58     -75.749842                                                       
!   129   59     -69.992000                                                       
!   129   60     -62.983000                                                       
!   129   61     -52.946000                                                       
!   130   48     -61.497000                                                       
!   130   49     -69.997161                                                       
!   130   50     -80.246161                                                       
!   130   51     -82.393930                                                       
!   130   52     -87.352930                                                       
!   130   53     -86.932580                                                       
!   130   54     -89.881796                                                       
!   130   55     -86.902636                                                       
!   130   56     -87.271214                                                       
!   130   57     -81.673000                                                       
!   130   58     -79.792359                                                       
!   130   59     -71.371000                                                       
!   130   60     -66.341000                                                       
!   130   61     -55.470000                                                       
!   130   62     -47.851000                                                       
!   131   49     -68.215710                                                       
!   131   50     -77.389307                                                       
!   131   51     -82.021307                                                       
!   131   52     -85.211307                                                       
!   131   53     -87.444761                                                       
!   131   54     -88.415608                                                       
!   131   55     -88.063214                                                       
!   131   56     -86.693390                                                       
!   131   57     -83.733390                                                       
!   131   58     -79.713390                                                       
!   131   59     -74.463390                                                       
!   131   60     -67.903390                                                       
!   131   61     -59.802000                                                       
!   131   62     -50.403000                                                       
!   132   49     -62.485535                                                       
!   132   50     -76.620535                                                       
!   132   51     -79.723535                                                       
!   132   52     -85.209535                                                       
!   132   53     -85.702535                                                       
!   132   54     -89.279535                                                       
!   132   55     -87.160068                                                       
!   132   56     -88.439611                                                       
!   132   57     -83.731611                                                       
!   132   58     -82.447000                                                       
!   132   59     -75.339000                                                       
!   132   60     -71.613000                                                       
!   132   61     -61.711000                                                       
!   132   62     -55.126000                                                       
!   132   63     -42.700000                                                       
!   133   49     -57.436000                                                       
!   133   50     -70.966712                                                       
!   133   51     -78.956712                                                       
!   133   52     -82.959712                                                       
!   133   53     -85.877712                                                       
!   133   54     -87.648300                                                       
!   133   55     -88.075660                                                       
!   133   56     -87.558217                                                       
!   133   57     -85.328217                                                       
!   133   58     -82.391000                                                       
!   133   59     -78.059000                                                       
!   133   60     -72.461000                                                       
!   133   61     -65.465000                                                       
!   133   62     -57.073000                                                       
!   133   63     -47.599000                                                       
!   134   49     -51.549000                                                       
!   134   50     -66.635740                                                       
!   134   51     -74.005740                                                       
!   134   52     -82.399438                                                       
!   134   53     -83.949438                                                       
!   134   54     -88.124438                                                       
!   134   55     -86.895877                                                       
!   134   56     -88.954546                                                       
!   134   57     -85.241369                                                       
!   134   58     -84.741369                                                       
!   134   59     -78.551000                                                       
!   134   60     -75.781000                                                       
!   134   61     -66.611000                                                       
!   134   62     -61.460000                                                       
!   134   63     -50.003000                                                       
!   135   50     -60.799000                                                       
!   135   51     -69.705584                                                       
!   135   52     -77.825584                                                       
!   135   53     -83.787584                                                       
!   135   54     -86.435645                                                       
!   135   55     -87.586595                                                       
!   135   56     -87.855940                                                       
!   135   57     -86.655940                                                       
!   135   58     -84.630357                                                       
!   135   59     -80.910357                                                       
!   135   60     -76.159000                                                       
!   135   61     -70.219000                                                       
!   135   62     -63.016000                                                       
!   135   63     -54.287000                                                       
!   136   50     -56.504000                                                       
!   136   51     -64.590000                                                       
!   136   52     -74.423420                                                       
!   136   53     -79.498246                                                       
!   136   54     -86.424442                                                       
!   136   55     -86.344134                                                       
!   136   56     -88.892358                                                       
!   136   57     -86.022358                                                       
!   136   58     -86.495190                                                       
!   136   59     -81.368844                                                       
!   136   60     -79.157844                                                       
!   136   61     -71.307844                                                       
!   136   62     -66.788000                                                       
!   136   63     -56.355000                                                       
!   136   64     -49.304000                                                       
!   137   50     -50.496000                                                       
!   137   51     -60.258000                                                       
!   137   52     -69.559519                                                       
!   137   53     -76.501119                                                       
!   137   54     -82.378579                                                       
!   137   55     -86.551145                                                       
!   137   56     -87.726774                                                       
!   137   57     -87.126667                                                       
!   137   58     -85.904567                                                       
!   137   59     -83.202567                                                       
!   137   60     -79.512567                                                       
!   137   61     -73.856000                                                       
!   137   62     -67.955543                                                       
!   137   63     -60.351000                                                       
!   137   64     -51.558000                                                       
!   138   51     -54.995000                                                       
!   138   52     -65.931000                                                       
!   138   53     -72.299139                                                       
!   138   54     -80.119139                                                       
!   138   55     -82.893139                                                       
!   138   56     -88.267172                                                       
!   138   57     -86.529421                                                       
!   138   58     -87.573860                                                       
!   138   59     -83.136860                                                       
!   138   60     -81.116860                                                       
!   138   61     -74.116860                                                       
!   138   62     -71.222000                                                       
!   138   63     -61.991000                                                       
!   138   64     -55.918000                                                       
!   138   65     -43.901000                                                       
!   139   51     -50.571000                                                       
!   139   52     -60.799000                                                       
!   139   53     -68.843542                                                       
!   139   54     -75.649542                                                       
!   139   55     -80.706565                                                       
!   139   56     -84.919280                                                       
!   139   57     -87.236114                                                       
!   139   58     -86.958114                                                       
!   139   59     -84.829114                                                       
!   139   60     -82.042114                                                       
!   139   61     -77.537722                                                       
!   139   62     -72.375210                                                       
!   139   63     -66.295210                                                       
!   139   64     -57.678000                                                       
!   139   65     -48.410000                                                       
!   140   52     -57.101000                                                       
!   140   53     -64.077000                                                       
!   140   54     -72.995897                                                       
!   140   55     -77.055897                                                       
!   140   56     -83.276031                                                       
!   140   57     -84.325762                                                       
!   140   58     -88.087615                                                       
!   140   59     -84.699615                                                       
!   140   60     -84.477342                                                       
!   140   61     -78.430247                                                       
!   140   62     -75.459386                                                       
!   140   63     -66.989386                                                       
!   140   64     -62.189386                                                       
!   140   65     -50.889386                                                       
!   140   66     -43.044000                                                       
!   141   52     -51.800000                                                       
!   141   53     -60.705000                                                       
!   141   54     -68.328538                                                       
!   141   55     -74.478538                                                       
!   141   56     -79.729877                                                       
!   141   57     -82.942993                                                       
!   141   58     -85.444905                                                       
!   141   59     -86.025576                                                       
!   141   60     -84.202574                                                       
!   141   61     -80.474888                                                       
!   141   62     -75.946081                                                       
!   141   63     -69.968353                                                       
!   141   64     -63.146000                                                       
!   141   65     -54.809000                                                       
!   141   66     -45.466000                                                       
!   142   52     -47.972000                                                       
!   142   53     -55.722000                                                       
!   142   54     -65.481242                                                       
!   142   55     -70.521242                                                       
!   142   56     -77.828012                                                       
!   142   57     -80.039086                                                       
!   142   58     -84.542631                                                       
!   142   59     -83.797314                                                       
!   142   60     -85.959517                                                       
!   142   61     -81.085853                                                       
!   142   62     -78.996944                                                       
!   142   63     -71.352399                                                       
!   142   64     -67.152399                                                       
!   142   65     -56.752399                                                       
!   142   66     -49.652399                                                       
!   142   67     -37.390000                                                       
!   143   53     -52.098000                                                       
!   143   54     -60.650000                                                       
!   143   55     -67.691387                                                       
!   143   56     -73.944606                                                       
!   143   57     -78.190856                                                       
!   143   58     -81.616413                                                       
!   143   59     -83.077858                                                       
!   143   60     -84.011780                                                       
!   143   61     -82.970420                                                       
!   143   62     -79.527634                                                       
!   143   63     -74.252511                                                       
!   143   64     -68.242511                                                       
!   143   65     -60.780000                                                       
!   143   66     -52.322000                                                       
!   143   67     -42.206000                                                       
!   144   53     -46.938000                                                       
!   144   54     -57.538000                                                       
!   144   55     -63.316084                                                       
!   144   56     -71.780481                                                       
!   144   57     -74.899869                                                       
!   144   58     -80.441308                                                       
!   144   59     -80.759963                                                       
!   144   60     -83.757479                                                       
!   144   61     -81.425820                                                       
!   144   62     -81.976369                                                       
!   144   63     -75.661412                                                       
!   144   64     -71.361412                                                       
!   144   65     -62.848000                                                       
!   144   66     -56.756000                                                       
!   144   67     -45.047000                                                       
!   144   68     -36.710000                                                       
!   145   54     -52.471000                                                       
!   145   55     -60.185471                                                       
!   145   56     -68.070025                                                       
!   145   57     -72.993377                                                       
!   145   58     -77.101730                                                       
!   145   59     -79.636301                                                       
!   145   60     -81.441582                                                       
!   145   61     -81.278541                                                       
!   145   62     -80.662142                                                       
!   145   63     -78.002099                                                       
!   145   64     -72.947615                                                       
!   145   65     -66.248000                                                       
!   145   66     -58.948000                                                       
!   145   67     -49.481000                                                       
!   145   68     -39.626000                                                       
!   146   54     -49.090000                                                       
!   146   55     -55.738704                                                       
!   146   56     -65.105231                                                       
!   146   57     -69.209858                                                       
!   146   58     -75.740025                                                       
!   146   59     -76.766257                                                       
!   146   60     -80.935509                                                       
!   146   61     -79.463724                                                       
!   146   62     -81.005724                                                       
!   146   63     -77.127958                                                       
!   146   64     -76.098070                                                       
!   146   65     -67.830797                                                       
!   146   66     -62.670797                                                       
!   146   67     -52.071000                                                       
!   146   68     -44.600000                                                       
!   146   69     -31.210000                                                       
!   147   54     -43.771000                                                       
!   147   55     -52.289934                                                       
!   147   56     -61.485563                                                       
!   147   57     -67.235563                                                       
!   147   58     -72.180563                                                       
!   147   59     -75.470563                                                       
!   147   60     -78.156253                                                       
!   147   61     -79.052253                                                       
!   147   62     -79.276392                                                       
!   147   63     -77.555055                                                       
!   147   64     -75.367684                                                       
!   147   65     -70.758904                                                       
!   147   66     -64.386257                                                       
!   147   67     -56.039000                                                       
!   147   68     -47.217000                                                       
!   147   69     -36.253000                                                       
!   148   55     -47.599766                                                       
!   148   56     -58.048483                                                       
!   148   57     -63.163483                                                       
!   148   58     -70.425837                                                       
!   148   59     -72.485837                                                       
!   148   60     -77.417837                                                       
!   148   61     -76.878251                                                       
!   148   62     -79.346590                                                       
!   148   63     -76.239256                                                       
!   148   64     -76.280245                                                       
!   148   65     -70.515356                                                       
!   148   66     -67.833356                                                       
!   148   67     -58.433000                                                       
!   148   68     -51.754000                                                       
!   148   69     -39.542000                                                       
!   148   70     -30.963000                                                       
!   149   55     -44.041000                                                       
!   149   56     -53.598000                                                       
!   149   57     -61.134000                                                       
!   149   58     -66.798164                                                       
!   149   59     -70.988164                                                       
!   149   60     -74.385196                                                       
!   149   61     -76.075853                                                       
!   149   62     -77.146768                                                       
!   149   63     -76.451500                                                       
!   149   64     -75.137623                                                       
!   149   65     -71.499950                                                       
!   149   66     -67.687950                                                       
!   149   67     -61.674261                                                       
!   149   68     -53.299588                                                       
!   149   69     -44.106000                                                       
!   149   70     -34.018000                                                       
!   150   55     -39.151000                                                       
!   150   56     -50.655000                                                       
!   150   57     -57.222000                                                       
!   150   58     -64.993681                                                       
!   150   59     -68.003681                                                       
!   150   60     -73.693681                                                       
!   150   61     -73.607132                                                       
!   150   62     -77.061132                                                       
!   150   63     -74.800546                                                       
!   150   64     -75.771944                                                       
!   150   65     -71.115684                                                       
!   150   66     -69.322027                                                       
!   150   67     -62.632796                                                       
!   150   68     -58.524796                                                       
!   150   69     -46.882000                                                       
!   150   70     -39.132000                                                       
!   150   71     -25.460000                                                       
!   151   55     -35.397000                                                       
!   151   56     -45.923000                                                       
!   151   57     -54.437000                                                       
!   151   58     -61.441000                                                       
!   151   59     -66.855300                                                       
!   151   60     -70.956788                                                       
!   151   61     -73.399207                                                       
!   151   62     -74.586250                                                       
!   151   63     -74.662939                                                       
!   151   64     -74.198821                                                       
!   151   65     -71.633583                                                       
!   151   66     -68.763221                                                       
!   151   67     -63.638923                                                       
!   151   68     -58.256000                                                       
!   151   69     -51.379000                                                       
!   151   70     -42.235827                                                       
!   151   71     -30.602000                                                       
!   152   56     -42.700000                                                       
!   152   57     -50.198000                                                       
!   152   58     -59.262000                                                       
!   152   59     -63.714000                                                       
!   152   60     -70.157856                                                       
!   152   61     -71.268076                                                       
!   152   62     -74.772647                                                       
!   152   63     -72.898338                                                       
!   152   64     -74.717096                                                       
!   152   65     -70.727096                                                       
!   152   66     -70.128564                                                       
!   152   67     -63.583214                                                       
!   152   68     -60.474024                                                       
!   152   69     -51.884000                                                       
!   152   70     -46.419000                                                       
!   152   71     -33.897000                                                       
!   153   56     -37.623000                                                       
!   153   57     -47.087000                                                       
!   153   58     -55.349000                                                       
!   153   59     -61.805000                                                       
!   153   60     -67.352098                                                       
!   153   61     -70.688098                                                       
!   153   62     -72.569048                                                       
!   153   63     -73.377294                                                       
!   153   64     -72.892857                                                       
!   153   65     -71.323686                                                       
!   153   66     -69.153298                                                       
!   153   67     -65.023389                                                       
!   153   68     -60.460355                                                       
!   153   69     -54.000905                                                       
!   153   70     -47.311000                                                       
!   153   71     -38.480000                                                       
!   154   57     -42.476000                                                       
!   154   58     -52.797000                                                       
!   154   59     -58.321000                                                       
!   154   60     -65.685878                                                       
!   154   61     -68.421001                                                       
!   154   62     -72.465282                                                       
!   154   63     -71.747955                                                       
!   154   64     -73.716308                                                       
!   154   65     -70.154308                                                       
!   154   66     -70.400400                                                       
!   154   67     -64.649151                                                       
!   154   68     -62.617538                                                       
!   154   69     -55.114326                                                       
!   154   70     -50.625625                                                       
!   154   71     -39.961000                                                       
!   154   72     -33.301000                                                       
!   155   57     -39.002000                                                       
!   155   58     -48.400000                                                       
!   155   59     -55.899000                                                       
!   155   60     -62.755159                                                       
!   155   61     -66.977159                                                       
!   155   62     -70.201159                                                       
!   155   63     -71.828023                                                       
!   155   64     -72.080112                                                       
!   155   65     -71.258898                                                       
!   155   66     -69.164398                                                       
!   155   67     -66.062398                                                       
!   155   68     -62.219810                                                       
!   155   69     -56.642687                                                       
!   155   70     -50.494000                                                       
!   155   71     -43.182822                                                       
!   155   72     -34.689000                                                       
!   156   58     -45.401000                                                       
!   156   59     -52.052000                                                       
!   156   60     -60.361000                                                       
!   156   61     -64.216854                                                       
!   156   62     -69.371854                                                       
!   156   63     -70.094116                                                       
!   156   64     -72.545159                                                       
!   156   65     -70.100736                                                       
!   156   66     -70.534324                                                       
!   156   67     -66.134324                                                       
!   156   68     -64.259103                                                       
!   156   69     -56.814703                                                       
!   156   70     -53.237567                                                       
!   156   71     -43.866000                                                       
!   156   72     -37.961000                                                       
!   156   73     -26.371000                                                       
!   157   58     -40.669000                                                       
!   157   59     -49.211000                                                       
!   157   60     -56.570000                                                       
!   157   61     -62.224000                                                       
!   157   62     -66.737338                                                       
!   157   63     -69.471338                                                       
!   157   64     -70.833880                                                       
!   157   65     -70.773828                                                       
!   157   66     -69.432382                                                       
!   157   67     -66.892382                                                       
!   157   68     -63.392333                                                       
!   157   69     -58.911333                                                       
!   157   70     -53.413115                                                       
!   157   71     -46.480113                                                       
!   157   72     -39.004000                                                       
!   157   73     -29.673000                                                       
!   158   59     -44.917000                                                       
!   158   60     -54.148000                                                       
!   158   61     -58.973000                                                       
!   158   62     -65.215806                                                       
!   158   63     -67.214806                                                       
!   158   64     -70.699888                                                       
!   158   65     -69.479885                                                       
!   158   66     -70.416616                                                       
!   158   67     -66.186616                                                       
!   158   68     -65.287000                                                       
!   158   69     -58.687000                                                       
!   158   70     -56.022000                                                       
!   158   71     -47.899295                                                       
!   158   72     -42.797225                                                       
!   158   73     -31.328000                                                       
!   158   74     -24.276000                                                       
!   159   59     -41.703000                                                       
!   159   60     -49.937000                                                       
!   159   61     -56.700000                                                       
!   159   62     -62.224000                                                       
!   159   63     -66.057353                                                       
!   159   64     -68.571851                                                       
!   159   65     -69.542412                                                       
!   159   66     -69.176777                                                       
!   159   67     -67.339054                                                       
!   159   68     -64.570486                                                       
!   159   69     -60.725048                                                       
!   159   70     -55.746428                                                       
!   159   71     -49.727694                                                       
!   159   72     -42.846000                                                       
!   159   73     -35.095811                                                       
!   159   74     -25.821000                                                       
!   160   60     -47.143000                                                       
!   160   61     -53.104000                                                       
!   160   62     -60.417000                                                       
!   160   63     -63.844211                                                       
!   160   64     -67.951903                                                       
!   160   65     -67.846267                                                       
!   160   66     -69.681592                                                       
!   160   67     -66.391592                                                       
!   160   68     -66.062547                                                       
!   160   69     -60.462547                                                       
!   160   70     -58.163000                                                       
!   160   71     -50.876000                                                       
!   160   72     -45.909990                                                       
!   160   73     -35.995000                                                       
!   160   74     -29.464000                                                       
!   160   75     -17.247000                                                       
!   161   60     -42.541000                                                       
!   161   61     -50.431000                                                       
!   161   62     -56.979000                                                       
!   161   63     -61.777000                                                       
!   161   64     -65.515980                                                       
!   161   65     -67.471557                                                       
!   161   66     -68.064634                                                       
!   161   67     -67.205734                                                       
!   161   68     -65.203315                                                       
!   161   69     -62.039315                                                       
!   161   70     -58.350900                                                       
!   161   71     -53.050900                                                       
!   161   72     -46.266506                                                       
!   161   73     -38.775302                                                       
!   161   74     -30.656000                                                       
!   161   75     -20.809000                                                       
!   162   61     -46.305000                                                       
!   162   62     -54.753000                                                       
!   162   63     -58.647000                                                       
!   162   64     -64.290578                                                       
!   162   65     -65.684473                                                       
!   162   66     -68.190258                                                       
!   162   67     -66.049975                                                       
!   162   68     -66.345722                                                       
!   162   69     -61.506402                                                       
!   162   70     -59.848000                                                       
!   162   71     -52.888000                                                       
!   162   72     -49.180103                                                       
!   162   73     -40.467304                                                       
!   162   74     -34.697864                                                       
!   162   75     -22.629000                                                       
!   162   76     -15.072000                                                       
!   163   61     -43.296000                                                       
!   163   62     -50.897000                                                       
!   163   63     -56.626000                                                       
!   163   64     -61.488000                                                       
!   163   65     -64.604742                                                       
!   163   66     -66.389866                                                       
!   163   67     -66.387301                                                       
!   163   68     -65.177303                                                       
!   163   69     -62.738303                                                       
!   163   70     -59.368303                                                       
!   163   71     -54.768303                                                       
!   163   72     -49.316000                                                       
!   163   73     -42.553760                                                       
!   163   74     -34.901000                                                       
!   163   75     -26.663200                                                       
!   163   76     -16.722000                                                       
!   164   62     -48.177000                                                       
!   164   63     -53.104000                                                       
!   164   64     -59.746000                                                       
!   164   65     -62.086625                                                       
!   164   66     -65.976625                                                       
!   164   67     -64.989789                                                       
!   164   68     -65.952563                                                       
!   164   69     -61.990011                                                       
!   164   70     -60.994000                                                       
!   164   71     -54.758000                                                       
!   164   72     -51.770000                                                       
!   164   73     -43.249000                                                       
!   164   74     -38.206329                                                       
!   164   75     -27.647000                                                       
!   164   76     -20.561000                                                       
!   165   62     -43.799000                                                       
!   165   63     -50.561000                                                       
!   165   64     -56.467000                                                       
!   165   65     -60.659000                                                       
!   165   66     -63.621191                                                       
!   165   67     -64.907266                                                       
!   165   68     -64.531286                                                       
!   165   69     -62.938746                                                       
!   165   70     -60.176746                                                       
!   165   71     -56.256746                                                       
!   165   72     -51.661000                                                       
!   165   73     -45.813000                                                       
!   165   74     -38.809795                                                       
!   165   75     -30.692473                                                       
!   165   76     -21.914000                                                       
!   165   77     -11.569000                                                       
!   166   63     -46.603000                                                       
!   166   64     -54.399000                                                       
!   166   65     -57.706000                                                       
!   166   66     -62.593368                                                       
!   166   67     -63.079584                                                       
!   166   68     -64.934464                                                       
!   166   69     -61.894849                                                       
!   166   70     -61.590725                                                       
!   166   71     -56.110725                                                       
!   166   72     -53.794000                                                       
!   166   73     -46.137000                                                       
!   166   74     -41.898691                                                       
!   166   75     -32.405393                                                       
!   166   76     -26.142458                                                       
!   166   77     -13.501000                                                       
!   167   63     -43.734000                                                       
!   167   64     -50.701000                                                       
!   167   65     -55.843000                                                       
!   167   66     -59.942538                                                       
!   167   67     -62.292538                                                       
!   167   68     -63.299248                                                       
!   167   69     -62.550890                                                       
!   167   70     -60.596598                                                       
!   167   71     -57.466598                                                       
!   167   72     -53.468000                                                       
!   167   73     -47.843000                                                       
!   167   74     -42.223000                                                       
!   167   75     -34.872000                                                       
!   167   76     -26.497000                                                       
!   167   77     -17.743489                                                       
!   168   64     -48.102000                                                       
!   168   65     -52.499000                                                       
!   168   66     -58.470000                                                       
!   168   67     -60.084685                                                       
!   168   68     -62.998997                                                       
!   168   69     -61.319891                                                       
!   168   70     -61.576900                                                       
!   168   71     -57.101900                                                       
!   168   72     -55.303000                                                       
!   168   73     -48.636000                                                       
!   168   74     -44.839000                                                       
!   168   75     -35.761000                                                       
!   168   76     -29.963446                                                       
!   168   77     -18.662000                                                       
!   168   78     -11.146000                                                       
!   169   64     -43.901000                                                       
!   169   65     -50.096000                                                       
!   169   66     -55.606785                                                       
!   169   67     -58.806785                                                       
!   169   68     -60.930800                                                       
!   169   69     -61.281940                                                       
!   169   70     -60.372800                                                       
!   169   71     -58.079800                                                       
!   169   72     -54.810434                                                       
!   169   73     -50.375000                                                       
!   169   74     -44.936000                                                       
!   169   75     -38.350000                                                       
!   169   76     -30.668313                                                       
!   169   77     -21.991762                                                       
!   169   78     -12.649000                                                       
!   170   65     -46.342000                                                       
!   170   66     -53.403000                                                       
!   170   67     -56.248302                                                       
!   170   68     -60.118302                                                       
!   170   69     -59.803883                                                       
!   170   70     -60.771916                                                       
!   170   71     -57.312778                                                       
!   170   72     -56.216000                                                       
!   170   73     -50.217000                                                       
!   170   74     -47.996164                                                       
!   170   75     -38.971000                                                       
!   170   76     -33.934586                                                       
!   170   77     -23.807402                                                       
!   170   78     -17.013811                                                       
!   171   65     -43.501000                                                       
!   171   66     -49.854000                                                       
!   171   67     -54.528508                                                       
!   171   68     -57.728508                                                       
!   171   69     -59.218961                                                       
!   171   70     -59.315389                                                       
!   171   71     -57.836544                                                       
!   171   72     -55.433000                                                       
!   171   73     -51.735000                                                       
!   171   74     -47.078000                                                       
!   171   75     -41.408000                                                       
!   171   76     -34.428000                                                       
!   171   77     -26.288000                                                       
!   171   78     -17.465000                                                       
!   171   79      -8.213000                                                       
!   172   66     -47.404000                                                       
!   172   67     -51.400000                                                       
!   172   68     -56.493101                                                       
!   172   69     -57.383638                                                       
!   172   70     -59.263786                                                       
!   172   71     -56.744520                                                       
!   172   72     -56.394520                                                       
!   172   73     -51.474520                                                       
!   172   74     -48.224520                                                       
!   172   75     -41.651000                                                       
!   172   76     -37.187000                                                       
!   172   77     -27.346000                                                       
!   172   78     -21.073989                                                       
!   172   79      -9.213000                                                       
!   173   66     -43.370000                                                       
!   173   67     -49.099000                                                       
!   173   68     -53.654000                                                       
!   173   69     -56.261916                                                       
!   173   70     -57.560027                                                       
!   173   71     -56.889216                                                       
!   173   72     -55.284000                                                       
!   173   73     -51.614000                                                       
!   173   74     -47.614000                                                       
!   173   75     -43.722000                                                       
!   173   76     -37.454000                                                       
!   173   77     -30.080000                                                       
!   173   78     -21.890439                                                       
!   173   79     -12.670051                                                       
!   174   67     -45.503000                                                       
!   174   68     -51.847000                                                       
!   174   69     -53.873303                                                       
!   174   70     -56.953303                                                       
!   174   71     -55.578958                                                       
!   174   72     -55.852224                                                       
!   174   73     -52.007224                                                       
!   174   74     -50.152000                                                       
!   174   75     -44.606000                                                       
!   174   76     -40.699105                                                       
!   174   77     -30.922000                                                       
!   174   78     -25.326129                                                       
!   174   79     -14.600364                                                       
!   175   67     -42.802000                                                       
!   175   68     -48.503000                                                       
!   175   69     -52.319311                                                       
!   175   70     -54.704311                                                       
!   175   71     -55.174334                                                       
!   175   72     -54.489605                                                       
!   175   73     -52.490000                                                       
!   175   74     -49.583000                                                       
!   175   75     -45.277000                                                       
!   175   76     -39.980000                                                       
!   175   77     -33.274000                                                       
!   175   78     -25.825000                                                       
!   175   79     -17.185000                                                       
!   175   80      -8.000000                                                       
!   176   68     -46.305000                                                       
!   176   69     -49.377174                                                       
!   176   70     -53.497174                                                       
!   176   71     -53.390993                                                       
!   176   72     -54.583838                                                       
!   176   73     -51.473838                                                       
!   176   74     -50.683000                                                       
!   176   75     -45.112000                                                       
!   176   76     -43.074194                                                       
!   176   77     -33.989000                                                       
!   176   78     -28.876000                                                       
!   176   79     -18.379000                                                       
!   176   80     -11.724482                                                       
!   177   68     -42.504000                                                       
!   177   69     -47.469000                                                       
!   177   70     -50.992651                                                       
!   177   71     -52.391884                                                       
!   177   72     -52.890209                                                       
!   177   73     -51.724209                                                       
!   177   74     -49.723000                                                       
!   177   75     -46.323000                                                       
!   177   76     -41.875000                                                       
!   177   77     -36.170000                                                       
!   177   78     -29.386000                                                       
!   177   79     -21.224000                                                       
!   177   80     -12.727118                                                       
!   177   81      -2.905000                                                       
!   178   69     -44.116000                                                       
!   178   70     -49.701349                                                       
!   178   71     -50.345971                                                       
!   178   72     -52.445216                                                       
!   178   73     -50.533216                                                       
!   178   74     -50.441916                                                       
!   178   75     -45.781916                                                       
!   178   76     -43.455842                                                       
!   178   77     -37.180630                                                       
!   178   78     -32.700815                                                       
!   178   79     -22.379000                                                       
!   178   80     -16.323195                                                       
!   178   81      -4.995000                                                       
!   179   69     -41.601000                                                       
!   179   70     -46.416000                                                       
!   179   71     -49.067170                                                       
!   179   72     -50.472927                                                       
!   179   73     -50.362041                                                       
!   179   74     -49.302356                                                       
!   179   75     -46.592356                                                       
!   179   76     -42.894000                                                       
!   179   77     -38.052000                                                       
!   179   78     -32.160000                                                       
!   179   79     -24.767000                                                       
!   179   80     -16.969000                                                       
!   179   81      -7.949000                                                       
!   180   70     -44.404000                                                       
!   180   71     -46.686502                                                       
!   180   72     -49.789502                                                       
!   180   73     -48.935420                                                       
!   180   74     -49.643282                                                       
!   180   75     -45.840974                                                       
!   180   76     -44.424580                                                       
!   180   77     -38.605190                                                       
!   180   78     -35.374663                                                       
!   180   79     -25.713000                                                       
!   180   80     -20.193000                                                       
!   180   81      -9.135000                                                       
!   181   70     -40.846000                                                       
!   181   71     -44.740000                                                       
!   181   72     -47.413853                                                       
!   181   73     -48.441085                                                       
!   181   74     -48.253195                                                       
!   181   75     -46.514522                                                       
!   181   76     -43.524522                                                       
!   181   77     -39.456072                                                       
!   181   78     -34.300000                                                       
!   181   79     -27.993000                                                       
!   181   80     -20.674000                                                       
!   181   81     -12.199000                                                       
!   181   82      -3.061000                                                       
!   182   71     -41.722000                                                       
!   182   72     -46.059677                                                       
!   182   73     -46.432720                                                       
!   182   74     -48.246241                                                       
!   182   75     -45.446241                                                       
!   182   76     -44.538241                                                       
!   182   77     -39.003801                                                       
!   182   78     -36.078959                                                       
!   182   79     -29.228959                                                       
!   182   80     -24.278959                                                       
!   182   81     -13.404000                                                       
!   182   82      -6.822167                                                       
!   183   71     -39.523000                                                       
!   183   72     -43.285577                                                       
!   183   73     -45.295577                                                       
!   183   74     -46.365611                                                       
!   183   75     -45.809611                                                       
!   183   76     -43.678000                                                       
!   183   77     -40.228000                                                       
!   183   78     -35.650000                                                       
!   183   79     -30.161000                                                       
!   183   80     -23.696000                                                       
!   183   81     -16.118000                                                       
!   183   82      -7.517000                                                       
!   184   71     -36.170000                                                       
!   184   72     -41.500026                                                       
!   184   73     -42.840026                                                       
!   184   74     -45.706026                                                       
!   184   75     -44.223334                                                       
!   184   76     -44.254521                                                       
!   184   77     -39.692521                                                       
!   184   78     -37.397909                                                       
!   184   79     -30.947909                                                       
!   184   80     -27.287909                                                       
!   184   81     -16.990000                                                       
!   184   82     -10.993000                                                       
!   185   72     -38.396000                                                       
!   185   73     -41.396439                                                       
!   185   74     -43.388439                                                       
!   185   75     -43.821433                                                       
!   185   76     -42.808642                                                       
!   185   77     -40.436000                                                       
!   185   78     -36.557611                                                       
!   185   79     -31.850611                                                       
!   185   80     -26.098000                                                       
!   185   81     -19.468000                                                       
!   185   82     -11.569000                                                       
!   185   83      -2.135000                                                       
!   186   72     -36.403000                                                       
!   186   73     -38.610327                                                       
!   186   74     -42.511327                                                       
!   186   75     -41.929771                                                       
!   186   76     -42.999289                                                       
!   186   77     -39.168289                                                       
!   186   78     -37.788521                                                       
!   186   79     -31.672961                                                       
!   186   80     -28.447803                                                       
!   186   81     -20.912148                                                       
!   186   82     -15.383260                                                       
!   186   83      -3.279000                                                       
!   187   73     -36.878000                                                       
!   187   74     -39.906720                                                       
!   187   75     -41.217871                                                       
!   187   76     -41.220534                                                       
!   187   77     -39.718125                                                       
!   187   78     -36.610000                                                       
!   187   79     -33.010000                                                       
!   187   80     -28.145000                                                       
!   187   81     -22.197000                                                       
!   187   82     -14.876000                                                       
!   187   83      -6.094000                                                       
!   188   73     -33.804000                                                       
!   188   74     -38.669148                                                       
!   188   75     -39.018148                                                       
!   188   76     -41.138501                                                       
!   188   77     -38.329144                                                       
!   188   78     -37.822659                                                       
!   188   79     -32.302659                                                       
!   188   80     -30.262659                                                       
!   188   81     -22.430000                                                       
!   188   82     -18.751803                                                       
!   188   83      -7.291000                                                       
!   189   74     -35.478534                                                       
!   189   75     -37.978534                                                       
!   189   76     -38.987800                                                       
!   189   77     -38.455352                                                       
!   189   78     -36.484845                                                       
!   189   79     -33.324845                                                       
!   189   80     -29.124845                                                       
!   189   81     -23.948287                                                       
!   189   82     -17.810000                                                       
!   189   83      -9.776000                                                       
!   190   74     -34.298032                                                       
!   190   75     -35.568032                                                       
!   190   76     -38.708032                                                       
!   190   77     -36.708032                                                       
!   190   78     -37.324891                                                       
!   190   79     -32.882891                                                       
!   190   80     -30.777891                                                       
!   190   81     -23.777891                                                       
!   190   82     -20.325188                                                       
!   190   83     -11.625037                                                       
!   190   84      -5.315209                                                       
!   191   75     -34.350148                                                       
!   191   76     -36.395374                                                       
!   191   77     -36.709063                                                       
!   191   78     -35.690511                                                       
!   191   79     -33.860511                                                       
!   191   80     -30.680511                                                       
!   191   81     -25.837821                                                       
!   191   82     -20.307000                                                       
!   191   83     -12.991000                                                       
!   191   84      -4.980000                                                       
!   192   75     -31.708000                                                       
!   192   76     -35.882031                                                       
!   192   77     -34.835823                                                       
!   192   78     -36.295511                                                       
!   192   79     -32.779170                                                       
!   192   80     -31.743846                                                       
!   192   81     -25.363846                                                       
!   192   82     -22.616847                                                       
!   192   83     -13.629000                                                       
!   192   84      -9.007182                                                       
!   193   76     -33.395840                                                       
!   193   77     -34.536346                                                       
!   193   78     -34.479707                                                       
!   193   79     -33.411060                                                       
!   193   80     -31.070752                                                       
!   193   81     -27.434000                                                       
!   193   82     -22.281000                                                       
!   193   83     -15.218534                                                       
!   193   84      -8.289000                                                       
!   193   85        .175000                                                       
!   194   76     -32.435255                                                       
!   194   77     -32.531855                                                       
!   194   78     -34.778645                                                       
!   194   79     -32.286611                                                       
!   194   80     -32.246611                                                       
!   194   81     -26.964000                                                       
!   194   82     -23.615180                                                       
!   194   83     -15.434680                                                       
!   194   84     -10.913283                                                       
!   194   85      -1.889000                                                       
!   195   76     -29.692402                                                       
!   195   77     -31.692402                                                       
!   195   78     -32.812385                                                       
!   195   79     -32.585585                                                       
!   195   80     -31.075585                                                       
!   195   81     -28.271000                                                       
!   195   82     -22.430112                                                       
!   195   83     -17.580112                                                       
!   195   84     -11.136000                                                       
!   195   85      -3.210000                                                       
!   196   76     -28.296408                                                       
!   196   77     -29.453923                                                       
!   196   78     -32.662940                                                       
!   196   79     -31.157245                                                       
!   196   80     -31.843261                                                       
!   196   81     -27.466000                                                       
!   196   82     -25.420000                                                       
!   196   83     -17.478451                                                       
!   196   84     -13.535135                                                       
!   196   85      -4.003000                                                       
!   196   86       1.040729                                                       
!   197   77     -28.283417                                                       
!   197   78     -30.438051                                                       
!   197   79     -31.156971                                                       
!   197   80     -30.557346                                                       
!   197   81     -28.376843                                                       
!   197   82     -24.796000                                                       
!   197   83     -19.622582                                                       
!   197   84     -13.445000                                                       
!   197   85      -5.690223                                                       
!   197   86       1.547000                                                       
!   198   77     -25.821000                                                       
!   198   78     -29.923300                                                       
!   198   79     -29.597990                                                       
!   198   80     -30.970466                                                       
!   198   81     -27.510466                                                       
!   198   82     -26.100000                                                       
!   198   83     -19.538646                                                       
!   198   84     -14.881188                                                       
!   198   85      -6.116487                                                       
!   198   86      -1.136443                                                       
!   199   77     -24.417100                                                       
!   199   78     -27.408077                                                       
!   199   79     -29.111031                                                       
!   199   80     -29.563297                                                       
!   199   81     -28.118225                                                       
!   199   82     -25.234743                                                       
!   199   83     -20.886930                                                       
!   199   84     -13.930950                                                       
!   199   85      -8.375463                                                       
!   199   86      -1.575000                                                       
!   200   78     -26.618475                                                       
!   200   79     -27.276109                                                       
!   200   80     -29.520227                                                       
!   200   81     -27.064187                                                       
!   200   82     -26.253633                                                       
!   200   83     -20.360609                                                       
!   200   84     -17.014000                                                       
!   200   85      -8.457141                                                       
!   200   86      -4.066765                                                       
!   200   87       6.054000                                                       
!   201   78     -23.756385                                                       
!   201   79     -26.416385                                                       
!   201   80     -27.679084                                                       
!   201   81     -27.196109                                                       
!   201   82     -25.293236                                                       
!   201   83     -21.451632                                                       
!   201   84     -16.572000                                                       
!   201   85     -10.724374                                                       
!   201   86      -4.159000                                                       
!   201   87       4.272688                                                       
!   202   78     -22.598000                                                       
!   202   79     -24.415916                                                       
!   202   80     -27.362070                                                       
!   202   81     -25.997463                                                       
!   202   82     -25.947892                                                       
!   202   83     -20.796062                                                       
!   202   84     -17.975000                                                       
!   202   85     -10.760031                                                       
!   202   86      -5.682688                                                       
!   202   87       3.696959                                                       
!   203   79     -23.159493                                                       
!   203   80     -25.283448                                                       
!   203   81     -25.775284                                                       
!   203   82     -24.800566                                                       
!   203   83     -21.547206                                                       
!   203   84     -17.313804                                                       
!   203   85     -12.251737                                                       
!   203   86      -4.876191                                                       
!   203   87       1.326288                                                       
!   203   88       8.579000                                                       
!   204   79     -20.767000                                                       
!   204   80     -24.707278                                                       
!   204   81     -24.359825                                                       
!   204   82     -25.123544                                                       
!   204   83     -20.674356                                                       
!   204   84     -18.343803                                                       
!   204   85     -11.865779                                                       
!   204   86      -8.044000                                                       
!   204   87       1.137549                                                       
!   204   88       5.993709                                                       
!   205   79     -18.993000                                                       
!   205   80     -22.303593                                                       
!   205   81     -23.834813                                                       
!   205   82     -23.783728                                                       
!   205   83     -21.075339                                                       
!   205   84     -17.544318                                                       
!   205   85     -13.007052                                                       
!   205   86      -7.761000                                                       
!   205   87      -1.244512                                                       
!   205   88       5.763000                                                       
!   206   80     -20.959849                                                       
!   206   81     -22.267062                                                       
!   206   82     -23.800598                                                       
!   206   83     -20.043090                                                       
!   206   84     -18.196508                                                       
!   206   85     -12.482724                                                       
!   206   86      -9.166000                                                       
!   206   87      -1.409457                                                       
!   206   88       4.158015                                                       
!   207   80     -16.229395                                                       
!   207   81     -21.044395                                                       
!   207   82     -22.467068                                                       
!   207   83     -20.068833                                                       
!   207   84     -17.159768                                                       
!   207   85     -13.249698                                                       
!   207   86      -8.637905                                                       
!   207   87      -2.925464                                                       
!   207   88       4.821978                                                       
!   207   89      11.615499                                                       
!   208   80     -13.097000                                                       
!   208   81     -16.762555                                                       
!   208   82     -21.763563                                                       
!   208   83     -18.884455                                                       
!   208   84     -17.483154                                                       
!   208   85     -12.498312                                                       
!   208   86      -9.658440                                                       
!   208   87      -2.669802                                                       
!   208   88       1.654000                                                       
!   208   89      11.283138                                                       
!   209   81     -13.647200                                                       
!   209   82     -17.628707                                                       
!   209   83     -18.272891                                                       
!   209   84     -16.379586                                                       
!   209   85     -12.893106                                                       
!   209   86      -8.964106                                                       
!   209   87      -3.804760                                                       
!   209   88       1.811000                                                       
!   209   89       8.913219                                                       
!   210   81      -9.253857                                                       
!   210   82     -14.742634                                                       
!   210   83     -14.806147                                                       
!   210   84     -15.968230                                                       
!   210   85     -11.987107                                                       
!   210   86      -9.613145                                                       
!   210   87      -3.354936                                                       
!   210   88        .416000                                                       
!   210   89       8.622654                                                       
!   210   90      14.635523                                                       
!   211   82     -10.496563                                                       
!   211   83     -11.868965                                                       
!   211   84     -12.447675                                                       
!   211   85     -11.661552                                                       
!   211   86      -8.769633                                                       
!   211   87      -4.164400                                                       
!   211   88        .832752                                                       
!   211   89       7.124247                                                       
!   211   90      15.189789                                                       
!   212   82      -7.556749                                                       
!   212   83      -8.130505                                                       
!   212   84     -10.384522                                                       
!   212   85      -8.630610                                                       
!   212   86      -8.673233                                                       
!   212   87      -3.544346                                                       
!   212   88       -.201676                                                       
!   212   89       7.276309                                                       
!   212   90      12.032000                                                       
!   213   82      -3.260000                                                       
!   213   83      -5.239806                                                       
!   213   84      -6.667147                                                       
!   213   85      -6.593906                                                       
!   213   86      -5.711591                                                       
!   213   87      -3.563108                                                       
!   213   88        .322155                                                       
!   213   89       6.123351                                                       
!   213   90      12.074000                                                       
!   213   91      19.732030                                                       
!   214   82       -.188025                                                       
!   214   83      -1.212186                                                       
!   214   84      -4.484259                                                       
!   214   85      -3.393980                                                       
!   214   86      -4.334917                                                       
!   214   87       -.973653                                                       
!   214   88        .084898                                                       
!   214   89       6.420855                                                       
!   214   90      10.666000                                                       
!   214   91      19.318465                                                       
!   215   83       1.706821                                                       
!   215   84       -.545288                                                       
!   215   85      -1.265672                                                       
!   215   86      -1.183747                                                       
!   215   87        .303694                                                       
!   215   88       2.518941                                                       
!   215   89       6.008911                                                       
!   215   90      10.923253                                                       
!   215   91      17.789308                                                       
!   216   83       5.775000                                                       
!   216   84       1.774680                                                       
!   216   85       2.243818                                                       
!   216   86        .240468                                                       
!   216   87       2.969484                                                       
!   216   88       3.277370                                                       
!   216   89       8.123808                                                       
!   216   90      10.293904                                                       
!   216   91      17.800520                                                       
!   217   84       5.825000                                                       
!   217   85       4.386981                                                       
!   217   86       3.646384                                                       
!   217   87       4.300196                                                       
!   217   88       5.874010                                                       
!   217   89       8.693330                                                       
!   217   90      12.171056                                                       
!   217   91      17.035691                                                       
!   218   84       8.351563                                                       
!   218   85       8.086725                                                       
!   218   86       5.203618                                                       
!   218   87       7.045191                                                       
!   218   88       6.635913                                                       
!   218   89      10.828658                                                       
!   218   90      12.358822                                                       
!   218   91      18.637242                                                       
!   218   92      21.878000                                                       
!   219   85      10.522600                                                       
!   219   86       8.825747                                                       
!   219   87       8.607788                                                       
!   219   88       9.379013                                                       
!   219   89      11.555105                                                       
!   219   90      14.457952                                                       
!   219   91      18.518422                                                       
!   219   92      23.208564                                                       
!   220   85      14.253000                                                       
!   220   86      10.604265                                                       
!   220   87      11.469463                                                       
!   220   88      10.260096                                                       
!   220   89      13.741495                                                       
!   220   90      14.655310                                                       
!   220   91      20.377819                                                       
!   220   92      23.019000                                                       
!   221   85      16.897000                                                       
!   221   86      14.396000                                                       
!   221   87      13.269739                                                       
!   221   88      12.954995                                                       
!   221   89      14.508712                                                       
!   221   90      16.926640                                                       
!   221   91      20.365941                                                       
!   221   92      24.546000                                                       
!   222   85      20.800000                                                       
!   222   86      16.366787                                                       
!   222   87      16.342089                                                       
!   222   88      14.309441                                                       
!   222   89      16.607466                                                       
!   222   90      17.189909                                                       
!   222   91      22.100000                                                       
!   222   92      24.284000                                                       
!   223   85      23.604000                                                       
!   223   86      20.297000                                                       
!   223   87      18.379037                                                       
!   223   88      17.229973                                                       
!   223   89      17.815779                                                       
!   223   90      19.370558                                                       
!   223   91      22.322084                                                       
!   223   92      25.823763                                                       
!   224   86      22.440000                                                       
!   224   87      21.643737                                                       
!   224   88      18.818043                                                       
!   224   89      20.221274                                                       
!   224   90      19.989160                                                       
!   224   91      23.860080                                                       
!   224   92      25.700045                                                       
!   225   86      26.492000                                                       
!   225   87      23.852683                                                       
!   225   88      21.987412                                                       
!   225   89      21.629824                                                       
!   225   90      22.301306                                                       
!   225   91      24.326123                                                       
!   225   92      27.371360                                                       
!   225   93      31.577352                                                       
!   226   86      28.774000                                                       
!   226   87      27.333218                                                       
!   226   88      23.662324                                                       
!   226   89      24.302533                                                       
!   226   90      23.185517                                                       
!   226   91      26.019186                                                       
!   226   92      27.329797                                                       
!   226   93      32.723000                                                       
!   227   86      32.981000                                                       
!   227   87      29.652399                                                       
!   227   88      27.172306                                                       
!   227   89      25.846141                                                       
!   227   90      25.801317                                                       
!   227   91      26.820644                                                       
!   227   92      29.006783                                                       
!   227   93      32.563406                                                       
!   228   86      35.475000                                                       
!   228   87      32.393712                                                       
!   228   88      28.936019                                                       
!   228   89      28.890119                                                       
!   228   90      26.763075                                                       
!   228   91      28.910721                                                       
!   228   92      29.217569                                                       
!   228   93      33.701000                                                       
!   228   94      36.074602                                                       
!   229   87      35.793000                                                       
!   229   88      32.434904                                                       
!   229   89      30.674904                                                       
!   229   90      29.579904                                                       
!   229   91      29.890335                                                       
!   229   92      31.201446                                                       
!   229   93      33.763734                                                       
!   229   94      37.389171                                                       
!   230   87      39.598000                                                       
!   230   88      34.544239                                                       
!   230   89      33.557201                                                       
!   230   90      30.857201                                                       
!   230   91      32.166869                                                       
!   230   92      31.603158                                                       
!   230   93      35.222197                                                       
!   230   94      36.929636                                                       
!   231   87      42.296000                                                       
!   231   88      38.396000                                                       
!   231   89      35.910488                                                       
!   231   90      33.810488                                                       
!   231   91      33.420973                                                       
!   231   92      33.803128                                                       
!   231   93      35.613955                                                       
!   231   94      38.432000                                                       
!   231   95      42.439000                                                       
!   232   87      46.253000                                                       
!   232   88      40.700000                                                       
!   232   89      39.143677                                                       
!   232   90      35.443677                                                       
!   232   91      35.938635                                                       
!   232   92      34.601531                                                       
!   232   93      37.352000                                                       
!   232   94      38.358400                                                       
!   232   95      43.398000                                                       
!   233   88      44.707000                                                       
!   233   89      41.498000                                                       
!   233   90      38.728649                                                       
!   233   91      37.483532                                                       
!   233   92      36.913421                                                       
!   233   93      37.941934                                                       
!   233   94      40.042657                                                       
!   233   95      43.289000                                                       
!   233   96      47.320000                                                       
!   234   88      47.085000                                                       
!   234   89      45.103000                                                       
!   234   90      40.608938                                                       
!   234   91      40.335850                                                       
!   234   92      38.140580                                                       
!   234   93      39.950426                                                       
!   234   94      40.338044                                                       
!   234   95      44.520000                                                       
!   234   96      46.798000                                                       
!   235   89      47.601000                                                       
!   235   90      44.250076                                                       
!   235   91      42.324062                                                       
!   235   92      40.914062                                                       
!   235   93      41.037778                                                       
!   235   94      42.179440                                                       
!   235   95      44.739000                                                       
!   235   96      48.057000                                                       
!   235   97      52.704000                                                       
!   236   89      51.398000                                                       
!   236   90      46.305000                                                       
!   236   91      45.340627                                                       
!   236   92      42.440627                                                       
!   236   93      43.370097                                                       
!   236   94      42.893512                                                       
!   236   95      46.174000                                                       
!   236   96      47.883000                                                       
!   236   97      53.403000                                                       
!   237   90      50.202000                                                       
!   237   91      47.636065                                                       
!   237   92      45.386065                                                       
!   237   93      44.867501                                                       
!   237   94      45.087818                                                       
!   237   95      46.547445                                                       
!   237   96      49.268000                                                       
!   237   97      53.214000                                                       
!   237   98      57.818000                                                       
!   238   90      52.390000                                                       
!   238   91      50.763664                                                       
!   238   92      47.303664                                                       
!   238   93      47.450729                                                       
!   238   94      46.158688                                                       
!   238   95      48.417037                                                       
!   238   96      49.384356                                                       
!   238   97      54.274000                                                       
!   238   98      57.203000                                                       
!   239   91      53.216000                                                       
!   239   92      50.568731                                                       
!   239   93      49.305273                                                       
!   239   94      48.583478                                                       
!   239   95      49.386389                                                       
!   239   96      51.186000                                                       
!   239   97      54.364000                                                       
!   239   98      58.292000                                                       
!   240   91      56.802000                                                       
!   240   92      52.709264                                                       
!   240   93      52.320918                                                       
!   240   94      50.121319                                                       
!   240   95      51.500271                                                       
!   240   96      51.715651                                                       
!   240   97      55.656000                                                       
!   240   98      58.027000                                                       
!   240   99      64.199000                                                       
!   241   92      56.197000                                                       
!   241   93      54.256039                                                       
!   241   94      52.951039                                                       
!   241   95      52.930224                                                       
!   241   96      53.697581                                                       
!   241   97      56.098000                                                       
!   241   98      59.351000                                                       
!   241   99      63.959000                                                       
!   242   92      58.614000                                                       
!   242   93      57.413000                                                       
!   242   94      54.713012                                                       
!   242   95      55.463975                                                       
!   242   96      54.799156                                                       
!   242   97      57.799000                                                       
!   242   98      59.325645                                                       
!   242   99      64.924000                                                       
!   242  100      68.400000                                                       
!   243   93      59.870000                                                       
!   243   94      57.749837                                                       
!   243   95      57.168280                                                       
!   243   96      57.177189                                                       
!   243   97      58.685576                                                       
!   243   98      60.939000                                                       
!   243   99      64.861000                                                       
!   243  100      69.406000                                                       
!   244   93      63.202000                                                       
!   244   94      59.799717                                                       
!   244   95      59.875894                                                       
!   244   96      58.447839                                                       
!   244   97      60.703482                                                       
!   244   98      61.469643                                                       
!   244   99      66.107000                                                       
!   244  100      69.002000                                                       
!   245   94      63.098143                                                       
!   245   95      61.893480                                                       
!   245   96      60.999422                                                       
!   245   97      61.809635                                                       
!   245   98      63.378000                                                       
!   245   99      66.432000                                                       
!   245  100      70.212000                                                       
!   245  101      75.467000                                                       
!   246   94      65.389406                                                       
!   246   95      64.988872                                                       
!   246   96      62.612736                                                       
!   246   97      63.962736                                                       
!   246   98      64.085667                                                       
!   246   99      67.966000                                                       
!   246  100      70.124380                                                       
!   246  101      76.320000                                                       
!   247   94      68.996000                                                       
!   247   95      67.148000                                                       
!   247   96      65.527622                                                       
!   247   97      65.482652                                                       
!   247   98      66.128652                                                       
!   247   99      68.604000                                                       
!   247  100      71.556000                                                       
!   247  101      76.200000                                                       
!   248   95      70.556000                                                       
!   248   96      67.386359                                                       
!   248   97      68.074000                                                       
!   248   98      67.233439                                                       
!   248   99      70.289000                                                       
!   248  100      71.896805                                                       
!   248  101      77.229000                                                       
!   249   95      73.104000                                                       
!   249   96      70.744222                                                       
!   249   97      69.843351                                                       
!   249   98      69.719351                                                       
!   249   99      71.170000                                                       
!   249  100      73.611000                                                       
!   249  101      77.316000                                                       
!   249  102      81.807000                                                       
!   250   96      72.983184                                                       
!   250   97      72.945777                                                       
!   250   98      71.166085                                                       
!   250   99      73.266000                                                       
!   250  100      74.067510                                                       
!   250  101      78.700000                                                       
!   250  102      81.499000                                                       
!   251   96      76.641333                                                       
!   251   97      75.221333                                                       
!   251   98      74.128333                                                       
!   251   99      74.504225                                                       
!   251  100      75.978663                                                       
!   251  101      79.101000                                                       
!   251  102      82.866000                                                       
!   251  103      87.896000                                                       
!   252   96      79.056000                                                       
!   252   97      78.528000                                                       
!   252   98      76.028139                                                       
!   252   99      77.288139                                                       
!   252  100      76.811050                                                       
!   252  101      80.695000                                                       
!   252  102      82.871198                                                       
!   252  103      88.799000                                                       
!   253   97      80.929000                                                       
!   253   98      79.295083                                                       
!   253   99      79.007422                                                       
!   253  100      79.341163                                                       
!   253  101      81.301000                                                       
!   253  102      84.439000                                                       
!   253  103      88.732000                                                       
!   253  104      93.782000                                                       
!   254   97      84.393000                                                       
!   254   98      81.334503                                                       
!   254   99      81.986388                                                       
!   254  100      80.898188                                                       
!   254  101      83.578000                                                       
!   254  102      84.718198                                                       
!   254  103      89.971000                                                       
!   254  104      93.304000                                                       
!   255   98      84.803000                                                       
!   255   99      84.082584                                                       
!   255  100      83.792964                                                       
!   255  101      84.835986                                                       
!   255  102      86.845454                                                       
!   255  103      90.140000                                                       
!   255  104      94.539000                                                       
!   255  105     100.041000                                                       
!   256   98      87.039000                                                       
!   256   99      87.180000                                                       
!   256  100      85.479952                                                       
!   256  101      87.609566                                                       
!   256  102      87.817402                                                       
!   256  103      91.997000                                                       
!   256  104      94.248151                                                       
!   256  105     100.704000                                                       
!   257   99      89.403000                                                       
!   257  100      88.583794                                                       
!   257  101      88.989933                                                       
!   257  102      90.217768                                                       
!   257  103      92.782000                                                       
!   257  104      96.011000                                                       
!   257  105     100.469000                                                       
!   258  100      90.419000                                                       
!   258  101      91.682582                                                       
!   258  102      91.473000                                                       
!   258  103      94.903000                                                       
!   258  104      96.473000                                                       
!   258  105     101.941000                                                       
!   258  106     105.399000                                                       
!   259  100      93.697000                                                       
!   259  101      93.617000                                                       
!   259  102      94.103000                                                       
!   259  103      95.935000                                                       
!   259  104      98.392000                                                       
!   259  105     102.205000                                                       
!   259  106     106.798000                                                       
!   260  101      96.545000                                                       
!   260  102      95.605000                                                       
!   260  103      98.340000                                                       
!   260  104      99.142000                                                       
!   260  105     103.794000                                                       
!   260  106     106.595916                                                       
!   260  107     113.459000                                                       
!   261  102      98.499000                                                       
!   261  103      99.615000                                                       
!   261  104     101.302000                                                       
!   261  105     104.426000                                                       
!   261  106     108.238000                                                       
!   261  107     113.456000                                                       
!   262  102     100.154000                                                       
!   262  103     102.177000                                                       
!   262  104     102.388000                                                       
!   262  105     106.333000                                                       
!   262  106     108.498000                                                       
!   262  107     114.582000                                                       
!   263  103     103.762000                                                       
!   263  104     104.830000                                                       
!   263  105     107.194000                                                       
!   263  106     110.208000                                                       
!   263  107     114.710000                                                       
!   263  108     119.893000                                                       
!   264  104     106.170000                                                       
!   264  105     109.425000                                                       
!   264  106     110.777000                                                       
!   264  107     116.186000                                                       
!   264  108     119.611504                                                       
!   265  105     110.530000                                                       
!   265  106     112.772000                                                       
!   265  107     116.621000                                                       
!   265  108     121.096000                                                       
!   265  109     127.211000                                                       
!   266  106     113.575000                                                       
!   266  107     118.308000                                                       
!   266  108     121.133000                                                       
!   266  109     128.490000                                                       
!   267  107     118.989000                                                       
!   267  108     122.746000                                                       
!   267  109     128.105000                                                       
!   267  110     134.094000                                                       
!   268  108     123.102000                                                       
!   268  109     129.306000                                                       
!   268  110     133.696000                                                       
!   269  108     124.927000                                                       
!   269  109     129.576000                                                       
!   269  110     135.200000                                                       
!   270  109     131.083000                                                       
!   270  110     134.718000                                                       
!   271  109     131.554000                                                       
!   271  110     136.070000                                                       
!   272  110     136.287000                                                       
!   272  111     142.963000                                                       
!   273  110     139.021000                                                       
!                                                                         ADGZ****
