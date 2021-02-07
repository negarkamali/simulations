c***********************************************************
      subroutine uelmat(rhs,amatrx,svars,energy,ndofel,nrhs,
     1     nsvars,props,nprops,coords,mcrd,nnode,u,du, 
     2     v,a,jtype,time,dtime,kstep,kinc,jelem,params, 
     3     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx, 
     4     ddlmag,mdload,pnewdt,jprops,njpro,period,
     5     materiallib)
c
      include 'aba_param.inc'
C
      dimension rhs(mlvarx,*), amatrx(ndofel, ndofel), props(*),
     1  svars(*), energy(*), coords(mcrd, nnode), u(ndofel),
     2  du(mlvarx,*), v(ndofel), a(ndofel), time(2), params(*),
     3  jdltyp(mdload,*), adlmag(mdload,*), ddlmag(mdload,*),
     4  predef(2, npredf, nnode), lflags(*), jprops(*)
      parameter (zero=0.d0, dmone=-1.0d0, one=1.d0, four=4.0d0, 
     1     fourth=0.25d0,gaussCoord=0.577350269d0,pi=3.1415927)

      parameter (ndim=2, ndof=10, nshr=1,nnodemax=4,
     1     ntens=4, ninptx=6, ninpty=6, ninpt=36, nsvint=3, en_n=4)
c
c      ndim  ... number of spatial dimensions
c      ndof  ... number of degrees of freedom per node
c      nshr  ... number of shear stress component
c      ntens ... total number of stress tensor components
c                (=ndi+nshr)
c      ninptx ... number of integration points in x direction
c      nsvint... number of state variables per integration pt 
c                (strain)
c     n_en......number of enrichment functions
c
      dimension stiff(ndofel,ndofel),rhs1(ndofel,1),CMAS(ndofel,ndofel),
     1 force(ndofel), shape(nnodemax), dshape(ndim,nnodemax),
     2 xjac(ndim,ndim),xjaci(ndim,ndim),bmat(ndofel),S(NDOFEL,3),
     3 statevLocal(nsvint),stress(4), ddsdde(ntens, ntens),
     4 stran(4), dstran(ntens),Arhs(ndofel,1),D(3,3),B(3,ndofel),
     5 AM(ndofel,ndofel),bmat1(ndim*nnodemax),B_cos_enrich_xx(nnodemax),
     6 AN(2,ndofel),aN_cos_enrich_x(nnodemax),aN_cos_enrich_y(nnodemax),
     7 AMASS(ndofel,ndofel),aN_sin_enrich_x(nnodemax),BB(2,ndofel/2),
     8 aN_sin_enrich_y(nnodemax),Bx_en_x(nnodemax),
     9 B_cos_enrich_yy(nnodemax),By_en_x(nnodemax)
       DIMENSION STRAIN1(3),STRESS1(3),CCBB(ndofel,ndofel),
     1 B_sin_enrich_xx(nnodemax),B_sin_enrich_xy(nnodemax),
     2 B_sin_enrich_yy(nnodemax),B_sin_enrich_yx(nnodemax),
     3 B_cos_enrich_xy(nnodemax),
     4 B_cos_enrich_yx(nnodemax)
c
      dimension predef_loc(npredf),dpredef_loc(npredf),
     1     defGrad(3,3),utmp(4),xdu(4),stiff_p(4,4),force_p(4)
      dimension coord24(2,4),coords_ip(3)
      data  coord24 /dmone, dmone,
     2                 one, dmone,
     3                 one,   one,
     4               dmone,   one/
c
      !data wght /one, one, one, one/
C GAUSS INTEGRATION VARIABLES (3 INTEG POINT)
      DIMENSION GAUSS3(3), WEIGHT3(3)
      DIMENSION GAUSS2(2), WEIGHT2(2)
      DIMENSION GAUSS6(6), WEIGHT6(6)
      DIMENSION GAUSS12(12), WEIGHT12(12)

c
c*************************************************************
c
c     U1 = first-order, plane strain, full integration
c
c     State variables: each integration point has nsvint SDVs
c
c       isvinc=(npt-1)*nsvint    ... integration point counter
c       statev(1+isvinc        ) ... strain
c
c*************************************************************
C           write(*,*)('--')

          Rho = PROPS(1)
          E = PROPS(2)
          ENU = PROPS(3)
         nnodex = PROPS(4)
         nnodey = PROPS(5)
            ELX = PROPS(6)
            ELY = PROPS(7)
			   Omega = PROPS(8)
                 ndofn = ndofel/4
                
       Elambda = (Enu*E)/((1+Enu)*(1-(2*Enu)))
         EMU = E/(2*(1+Enu))
       
      if (lflags(3).eq.4) then
        do i=1, ndofel
          do j=1, ndofel
            amatrx(i,j) = zero
            
          end do
          amatrx(i,i) = one
       
        end do
        goto 999
      end if
      !      write(*,*) 'pnewdt'
      !WRITE(*,*)(pnewdt)
c
c     PRELIMINARIES
c
      !pnewdt = 0.5
      pnewdtLocal = pnewdt
      if(jtype .ne. 1) then
        write(7,*)'Incorrect element type'
        call xit
      endif 
      if(nsvars .lt. ninpt*nsvint) then
        write(7,*)'Increase the number of SDVs to', ninpt*nsvint
        call xit
      endif 
      thickness = one
c
c     INITIALIZE RHS AND LHS
c
      do k1=1, ndof*nnode
        rhs(k1, 1)= zero
        do k2=1, ndof*nnode 
          amatrx(k1, k2)= zero
          AM(k1, k2)= zero
          AMASS(k1, k2)= zero
          CMAS(k1, k2)= zero
        end do
      end do
c Gaussian integration (2 gauss points)
      GAUSS2(1) = -0.577350269d0
      GAUSS2(2) = 0.577350269d0
    
      
      WEIGHT2(1) = one
      WEIGHT2(2) = one
    
c Gaussian integration (3 gauss points)
      GAUSS3(1) = -SQRT(0.6)
      GAUSS3(2) = ZERO
      GAUSS3(3) = SQRT(0.6)
      
      WEIGHT3(1) = 0.55555555555555
      WEIGHT3(2) = 0.88888888888888
      WEIGHT3(3) = 0.55555555555555
c Gaussian integration (6 gauss points)
      GAUSS6(1) = -0.932469514203152
      GAUSS6(2) = -0.6612093864662646
      GAUSS6(3) = -0.2386191860831968
      GAUSS6(4) = 0.2386191860831968
      GAUSS6(5) = 0.6612093864662646
      GAUSS6(6) = 0.932469514203152
      
      WEIGHT6(1) = 0.1713244923791709
      WEIGHT6(2) = 0.3607615730481379
      WEIGHT6(3) = 0.4679139345726913
      WEIGHT6(4) = 0.4679139345726913
      WEIGHT6(5) = 0.3607615730481379
      WEIGHT6(6) = 0.1713244923791709
c Gaussian integration (12 gauss points)
      GAUSS12(1) = -0.981560634246732
      GAUSS12(2) = -0.904117256370452
      GAUSS12(3) = -0.7699026741943177
      GAUSS12(4) = -0.5873179542866143
      GAUSS12(5) = -0.3678314989981804
      GAUSS12(6) = -0.12523340851114688
      GAUSS12(7) = 0.12523340851114688
      GAUSS12(8) = 0.3678314989981804
      GAUSS12(9) = 0.5873179542866143
      GAUSS12(10) = 0.7699026741943177
      GAUSS12(11) = 0.904117256370452
      GAUSS12(12) = 0.981560634246732
      
      WEIGHT12(1) = 0.04717533638647547
      WEIGHT12(2) = 0.1069393259953637
      WEIGHT12(3) = 0.1600783285433586
      WEIGHT12(4) = 0.2031674267230672
      WEIGHT12(5) = 0.2334925365383534
      WEIGHT12(6) = 0.2491470458134027
      WEIGHT12(7) = 0.2491470458134027
      WEIGHT12(8) = 0.2334925365383534
      WEIGHT12(9) = 0.2031674267230672
      WEIGHT12(10) = 0.1600783285433586
      WEIGHT12(11) = 0.1069393259953637
      WEIGHT12(12) = 0.04717533638647547
c
C LOOP OVER INTEGRATION POINTS
C 						DO 100 iintp = 1,ninpt
C 						IF (NINTP.EQ.3.AND.INTS.EQ.1) THEN
C 						POINT = GAUSS3(IINTP)
C 						WEIGHT = WEIGHT3(IINTP)
C 						ELSE IF (NINTP.EQ.6.AND.INTS.EQ.1) THEN
C 						POINT = GAUSS6(IINTP)
C 						WEIGHT = WEIGHT6(IINTP)
C 						ELSE IF (NINTP.EQ.12.AND.INTS.EQ.1) THEN
C 						POINT = GAUSS12(IINTP)
C 						WEIGHT = WEIGHT12(IINTP)
C 						ELSE IF (NINTP.EQ.3.AND.INTS.EQ.2) THEN
C 						POINT = COTNEW(IINTP)
C 						WEIGHT = CWEIGHT(IINTP)
C 						ELSE IF (NINTP.EQ.6.AND.INTS.EQ.2) THEN
C 						POINT = COTNEW6(IINTP)
C 						WEIGHT = CW6(IINTP)
C 						ELSE IF (NINTP.EQ.12.AND.INTS.EQ.2) THEN
C 						POINT = COTNEW12(IINTP)
C 						WEIGHT = CW12(IINTP)
C 						ELSE
C 						WRITE(7,*) ’Unspecified integration required’
C 						CALL FLUSH_(7)
C 						CALL XIT
C 						END IF

c     LOOP OVER INTEGRATION POINTS
c
	     if (kinc.eq.1) then
            write(7,*)('')
		    endif
      do kintkx = 1, ninptx
          do kintky = 1, ninpty
              

              
      !        write(*,*) 'kintkx'
      !WRITE(*,*)(kintkx)
      !
      !write(*,*) 'kintky'
      !WRITE(*,*)(kintky)
              
              kintk = ((kintkx - one)* ninptx )+ kintky 
              
              
      !        write(*,*) 'kintk'
      !WRITE(*,*)(kintk)
c
c       EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES
c
c       determine (g,h)
c
        g = GAUSS6(kintkx)
        h = GAUSS6(kintky)
        
      !  write(*,*) 'g'
      !WRITE(*,*)(g)
      !
      !write(*,*) 'h'
      !WRITE(*,*)(h)
c
c       shape functions
        shape(1) = (one - g)*(one - h)/four;
        shape(2) = (one + g)*(one - h)/four;
        shape(3) = (one + g)*(one + h)/four;
        shape(4) = (one - g)*(one + h)/four;
              !
              !write(*,*) 'shape'
              !WRITE(*,*)(shape)
             
c
c       derivative d(Ni)/d(g)
        dshape(1,1) = -(one - h)/four;
        dshape(1,2) =  (one - h)/four;
        dshape(1,3) =  (one + h)/four;
        dshape(1,4) = -(one + h)/four;
c
c       derivative d(Ni)/d(h)
        dshape(2,1) = -(one - g)/four;
        dshape(2,2) = -(one + g)/four;
        dshape(2,3) =  (one + g)/four;
        dshape(2,4) =  (one - g)/four;
c
c       compute coordinates at the integration point
c
        do k1=1, 3
          coords_ip(k1) = zero
        end do
        do k1=1,nnode
          do k2=1,mcrd
            coords_ip(k2)=coords_ip(k2)+shape(k1)*coords(k2,k1)
          end do
        end do       
c
c       INTERPOLATE FIELD VARIABLES
c
c ------ ENRICHMENTS------------------------------------

c Shape Functions
      ! cutoff numbers
       bn = one
       bm = one
      ! if 2D 
         xx  = coords_ip(1)
         yy  = coords_ip(2) 
       ! if 2D and rectangular length in x and y directions
        ALX = ELX/(nnodex-one)
		   	ALY = ELY/(nnodey-one)
C 			Frequencies:
C           Bathe's 
            wx = ((pi)/(ALX))*bn
			wy = ((pi)/(ALY))*bm
C          Ours
            wx = Omega
			wy = Omega
			
        cos_enrich_x = cos(wx*xx) 
        sin_enrich_x = sin(wx*xx) 
     
        aN_cos_enrich_x = shape * cos_enrich_x 
        aN_sin_enrich_x = shape * sin_enrich_x 
        
      
        
        cos_enrich_y = cos(wy*yy) 
        sin_enrich_y = sin(wy*yy) 
            
            aN_cos_enrich_y = shape * cos_enrich_y 
            aN_sin_enrich_y = shape * sin_enrich_y
            
! en_n => number of enrichment functions per node
     
      DO I=1,2
           DO K=1,ndofel
             AN(I,K)=0.0D0
          ENDDO                    
         ENDDO
         !
             AN(1,1)          =shape(1)
             AN(1,1+ndofn)    =shape(2)
             AN(1,1+2*(ndofn))=shape(3)
             AN(1,1+3*(ndofn))=shape(4)
             
             AN(2,2)          =shape(1)
             AN(2,2+ndofn)    =shape(2)
             AN(2,2+2*(ndofn))=shape(3)
             AN(2,2+3*(ndofn))=shape(4)
             
! first row in x direction  
             
                          AN(1,3)=aN_cos_enrich_x(1)
             AN(1,3 + 2*(en_n+1))=aN_cos_enrich_x(2)
             AN(1,3 + 4*(en_n+1))=aN_cos_enrich_x(3)
             AN(1,3 + 6*(en_n+1))=aN_cos_enrich_x(4)
             
                          AN(1,5)=aN_sin_enrich_x(1)
             AN(1,5 + 2*(en_n+1))=aN_sin_enrich_x(2)
             AN(1,5 + 4*(en_n+1))=aN_sin_enrich_x(3)
             AN(1,5 + 6*(en_n+1))=aN_sin_enrich_x(4)
             
                          AN(1,7)=aN_cos_enrich_y(1)
             AN(1,7 + 2*(en_n+1))=aN_cos_enrich_y(2)
             AN(1,7 + 4*(en_n+1))=aN_cos_enrich_y(3)
             AN(1,7 + 6*(en_n+1))=aN_cos_enrich_y(4)
             
                          AN(1,9)=aN_sin_enrich_y(1)
             AN(1,9 + 2*(en_n+1))=aN_sin_enrich_y(2)
             AN(1,9 + 4*(en_n+1))=aN_sin_enrich_y(3)
             AN(1,9 + 6*(en_n+1))=aN_sin_enrich_y(4)
             
      ! second row in y direction             
                          AN(2,4)=aN_cos_enrich_x(1)
             AN(2,4 + 2*(en_n+1))=aN_cos_enrich_x(2)
             AN(2,4 + 4*(en_n+1))=aN_cos_enrich_x(3)
             AN(2,4 + 6*(en_n+1))=aN_cos_enrich_x(4)
             
                          AN(2,6)=aN_sin_enrich_x(1)
             AN(2,6 + 2*(en_n+1))=aN_sin_enrich_x(2)
             AN(2,6 + 4*(en_n+1))=aN_sin_enrich_x(3)
             AN(2,6 + 6*(en_n+1))=aN_sin_enrich_x(4)
             
                          AN(2,8)=aN_cos_enrich_y(1)
             AN(2,8 + 2*(en_n+1))=aN_cos_enrich_y(2)
             AN(2,8 + 4*(en_n+1))=aN_cos_enrich_y(3)
             AN(2,8 + 6*(en_n+1))=aN_cos_enrich_y(4)
             
                         AN(2,10)=aN_sin_enrich_y(1)
             AN(2,10+ 2*(en_n+1))=aN_sin_enrich_y(2)
             AN(2,10+ 4*(en_n+1))=aN_sin_enrich_y(3)
             AN(2,10+ 6*(en_n+1))=aN_sin_enrich_y(4)

c----------------------------------------------------------
        if(npredf.gt.0) then

          do k1=1,npredf
            predef_loc(k1) = zero
            dpredef_loc(k1) = zero
            do k2=1,nnode
              predef_loc(k1) = 
     1             predef_loc(k1)+
     2             (predef(1,k1,k2)-predef(2,k1,k2))*shape(k2)
             dpredef_loc(k1) = 
     1             dpredef_loc(k1)+predef(2,k1,k2)*shape(k2)
            end do
          end do
        end if
c
c       FORM B-MATRIX for standard FEM
c
        djac = one
c
        do i = 1, ndim
          do j = 1, ndim
            xjac(i,j)  = zero
            xjaci(i,j) = zero
          end do
        end do
c     
        do inod= 1, nnode
          do idim = 1, ndim
            do jdim = 1, ndim
              xjac(jdim,idim) = xjac(jdim,idim) + 
     1             dshape(jdim,inod)*coords(idim,inod)
            end do
          end do 
        end do
        djac = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
        if (djac .gt. zero) then
        ! jacobian is positive - o.k.
          xjaci(1,1) =  xjac(2,2)/djac
          xjaci(2,2) =  xjac(1,1)/djac
          xjaci(1,2) = -xjac(1,2)/djac
          xjaci(2,1) = -xjac(2,1)/djac
        else
          ! negative or zero jacobian
          write(7,*)'WARNING: element',jelem,'has neg. 
     1         Jacobian'
          pnewdt = fourth
        endif
        

        if (pnewdt .lt. pnewdtLocal) pnewdtLocal = pnewdt
c
        do i = 1, nnode*ndim
          bmat1(i) = zero
        end do

        do inod = 1, nnode
          do ider = 1, ndim
            do idim = 1, ndim
              irow = idim + (inod - 1)*ndim
              bmat1(irow) = bmat1(irow) + 
     1             xjaci(idim,ider)*dshape(ider,inod)      
            end do
          end do
        end do 
c ------------------- FORM ENRICHED B---------------------
                do i = 1, ndofel
          bmat(i) = zero
        end do
        
c****** change this part to change the shape functions only ****
        bmat(1)         = bmat1(1)
        bmat(1+ndofn)   = bmat1(3)
        bmat(1+2*ndofn) = bmat1(5)
        bmat(1+3*ndofn) = bmat1(7)
        
        bmat(2)         = bmat1(2)
        bmat(2+ndofn)   = bmat1(4)
        bmat(2+2*ndofn) = bmat1(6)
        bmat(2+3*ndofn) = bmat1(8)
        
c******* Sin x ENRICHMENT B MAT Contribution     
       B_sin_enrich_xx(1)=bmat1(1)*sin_enrich_x
     &  +shape(1)*wx*cos_enrich_x
        
                B_sin_enrich_xx(2)=bmat1(3)*sin_enrich_x
     &  +shape(2)*wx*cos_enrich_x
     
      B_sin_enrich_xx(3)=bmat1(5)*sin_enrich_x
     &  +shape(3)*wx*cos_enrich_x
     
      B_sin_enrich_xx(4)=bmat1(7)*sin_enrich_x
     &  +shape(4)*wx*cos_enrich_x
     
          ! in y direction
       B_sin_enrich_xy(1) = bmat1(2)*sin_enrich_x  
       B_sin_enrich_xy(2) = bmat1(4)*sin_enrich_x
       B_sin_enrich_xy(3) = bmat1(6)*sin_enrich_x
       B_sin_enrich_xy(4) = bmat1(8)*sin_enrich_x


        B_sin_enrich_yx(1) = bmat1(1)*sin_enrich_y  
       B_sin_enrich_yx(2) = bmat1(3)*sin_enrich_y
       B_sin_enrich_yx(3) = bmat1(5)*sin_enrich_y
       B_sin_enrich_yx(4) = bmat1(7)*sin_enrich_y
       
c******* Sin y ENRICHMENT B MAT Contribution  
    
              B_sin_enrich_yy(1)=bmat1(2)*sin_enrich_y
     &  +shape(1)*wy*cos_enrich_y
        
                B_sin_enrich_yy(2)=bmat1(4)*sin_enrich_y
     & +shape(2)*wy*cos_enrich_y
                
                        B_sin_enrich_yy(3)=bmat1(6)*sin_enrich_y
     &  +shape(3)*wy*cos_enrich_y
                        
             B_sin_enrich_yy(4)=bmat1(8)*sin_enrich_y
     &   +shape(4)*wy*cos_enrich_y 
 
       
c******* cos x ENRICHMENT B MAT Contribution  
       
              B_cos_enrich_xx(1)=bmat1(1)*cos_enrich_x
     &  -shape(1)*wx*sin_enrich_x
        
       B_cos_enrich_xx(2)=bmat1(3)*cos_enrich_x
     &  -shape(2)*wx*sin_enrich_x
     
      B_cos_enrich_xx(3)=bmat1(5)*cos_enrich_x
     &  -shape(3)*wx*sin_enrich_x
     
      B_cos_enrich_xx(4)=bmat1(7)*cos_enrich_x
     &  -shape(4)*wx*sin_enrich_x
     
          ! in y direction
       B_cos_enrich_xy(1) = bmat1(2)*cos_enrich_x  
       B_cos_enrich_xy(2) = bmat1(4)*cos_enrich_x
       B_cos_enrich_xy(3) = bmat1(6)*cos_enrich_x
       B_cos_enrich_xy(4) = bmat1(8)*cos_enrich_x


c******* cos y ENRICHMENT B MAT Contribution 
       
         B_cos_enrich_yy(1)=bmat1(2)*cos_enrich_y
     &  -shape(1)*wy*sin_enrich_y
        
      B_cos_enrich_yy(2)=bmat1(4)*cos_enrich_y
     & -shape(2)*wy*sin_enrich_y
                
      B_cos_enrich_yy(3)=bmat1(6)*cos_enrich_y
     &  -shape(3)*wy*sin_enrich_y
                        
        B_cos_enrich_yy(4)=bmat1(8)*cos_enrich_y
     &   -shape(4)*wy*sin_enrich_y 
	 
c******* cos y ENRICHMENT B MAT Contribution 
 
       B_cos_enrich_yx(1) = bmat1(1)*cos_enrich_y  
       B_cos_enrich_yx(2) = bmat1(3)*cos_enrich_y
       B_cos_enrich_yx(3) = bmat1(5)*cos_enrich_y
       B_cos_enrich_yx(4) = bmat1(7)*cos_enrich_y
c -------------
       bmat(3)              = B_cos_enrich_xx(1)
       bmat(3 + 2*(en_n+1)) = B_cos_enrich_xx(2)
       bmat(3 + 4*(en_n+1)) = B_cos_enrich_xx(3)
       bmat(3 + 6*(en_n+1)) = B_cos_enrich_xx(4)
c---------------
       bmat(4)              = B_cos_enrich_xy(1)
       bmat(4 + 2*(en_n+1)) = B_cos_enrich_xy(2)
       bmat(4 + 4*(en_n+1)) = B_cos_enrich_xy(3)
       bmat(4 + 6*(en_n+1)) = B_cos_enrich_xy(4) 
c----------------
c----------------------
       bmat(5)              = B_sin_enrich_xx(1)
       bmat(5 + 2*(en_n+1)) = B_sin_enrich_xx(2)
       bmat(5 + 4*(en_n+1)) = B_sin_enrich_xx(3)
       bmat(5 + 6*(en_n+1)) = B_sin_enrich_xx(4)
c---------------
       bmat(6)              = B_sin_enrich_xy(1)
       bmat(6 + 2*(en_n+1)) = B_sin_enrich_xy(2)
       bmat(6 + 4*(en_n+1)) = B_sin_enrich_xy(3)
       bmat(6 + 6*(en_n+1)) = B_sin_enrich_xy(4) 
c -------------
       bmat(7)              = B_cos_enrich_yx(1)
       bmat(7 + 2*(en_n+1)) = B_cos_enrich_yx(2)
       bmat(7 + 4*(en_n+1)) = B_cos_enrich_yx(3)
       bmat(7 + 6*(en_n+1)) = B_cos_enrich_yx(4)
c---------------
       bmat(8)              = B_cos_enrich_yy(1)
       bmat(8 + 2*(en_n+1)) = B_cos_enrich_yy(2)
       bmat(8 + 4*(en_n+1)) = B_cos_enrich_yy(3)
       bmat(8 + 6*(en_n+1)) = B_cos_enrich_yy(4) 
c----------------

c----------------------
       bmat(9)              =  B_sin_enrich_yx(1)
       bmat(9 + 2*(en_n+1)) =  B_sin_enrich_yx(2)
       bmat(9 + 4*(en_n+1)) =  B_sin_enrich_yx(3)
       bmat(9 + 6*(en_n+1)) =  B_sin_enrich_yx(4)
c---------------
       bmat(10)             = B_sin_enrich_yy(1)
       bmat(10+ 2*(en_n+1)) = B_sin_enrich_yy(2)
       bmat(10+ 4*(en_n+1)) = B_sin_enrich_yy(3)
       bmat(10+ 6*(en_n+1)) = B_sin_enrich_yy(4)  
             !
C               write(*,*) 'bmat '
C                WRITE(*,*)(bmat )
C               write(*,*) 'cos_enrich_y'
C               WRITE(*,*)(cos_enrich_y)    
       
c just to check the B matrix with MATLAB
       
       gg = ndofel*0.5D0       
       !             
           DO K=1,gg
             BB(1,K)=bmat(K*2-1)
             BB(2,K)=bmat(K*2)
             
      ENDDO
       
       
         !
              B=0.0D0
        ff = ndofel*0.5D0       
            !             
           DO K=1,ff
             B(1,K*2-1)=BB(1,K)
             B(2,K*2)=BB(2,K)
             B(3,K*2-1)=BB(2,K)
             B(3,K*2)=BB(1,K)
      ENDDO


           !WRITE(*,*) "The Stress for 4 Gauss points are"            
           !DO K1=1,ndofel
           ! WRITE(*,*)(S(K1,K2),K2=1,3) 
           !ENDDO         
    
            
c-------------------CALCULATE Stress and Strain                          
c============================================================
      !   DO I=1,3
      !     STRAIN1(I)=0.0D0
      !       STRESS1(I)=0.0D0
      !ENDDO
      !   
      !!WRITE(*,*)('kintk')
      !!      WRITE(*,*)(kintk)
      !   
      !
      !      STRAIN1=MATMUL(U,B)
      !      
      !
      !      STRESS1=MATMUL(D,STRAIN1)
          !  
      !                   WRITE(*,*)"STRAIN1"   
      !     DO K1=1,3
      !    WRITE(*,*)(STRAIN1(K1))
      !ENDDO
      !     WRITE(*,*)"STRESS1"   
      !     DO K1=1,3
      !    WRITE(*,*)(STRESS1(K1))
      !     ENDDO
      !      WRITE(*,*)('U')
      !      WRITE(*,*)(U)
            
           !
           !DO I=1,3
           !   SVARS(I)=STRAIN(I)
           !   SVARS(I+4)=STRESS(I)
           !ENDDO
           !SVARS(4)=ENU*(STRESS(1)+STRESS(2))
           ! WRITE(*,*) "The Stress Vector is" 
           !DO K1=1,7
           ! WRITE(*,100),SVARS(K1) 
           !ENDDO          
c ==============================================================
       
       !IF (KINC. EQ. 1) THEN
       !   write(*,*) 'BB'
       !       WRITE(*,*)(BB)
       !
       !ENDIf
       
        do nodi = 1, nnode
           
          incr_row = (nodi - 1)*ndof

          do i = 1, ndof
            xdu(i)= du(i + incr_row,1)
            utmp(i) = u(i + incr_row)
          end do
c***********************************************************
c******* correct from here**********************************

            dNidx   = bmat(1 + (nodi-1)*(ndim + 2*en_n))
            dNidy   = bmat(2 + (nodi-1)*(ndim + 2*en_n))
            dNicxdx = bmat(3 + (nodi-1)*(ndim + 2*en_n))
            dNicxdy = bmat(4 + (nodi-1)*(ndim + 2*en_n))
            dNisxdx = bmat(5 + (nodi-1)*(ndim + 2*en_n))
            dNisxdy = bmat(6 + (nodi-1)*(ndim + 2*en_n))
		    dNicydx = bmat(7 + (nodi-1)*(ndim + 2*en_n))
		    dNicydy = bmat(8 + (nodi-1)*(ndim + 2*en_n))
		    dNisydx = bmat(9 + (nodi-1)*(ndim + 2*en_n))
		    dNisydy = bmat(10 +(nodi-1)*(ndim + 2*en_n))

          dstran(1) = dstran(1) + dNidx*xdu(1) + dNicxdx*xdu(3)
     1    + dNisxdx*xdu(5) + dNicydx*xdu(7) + dNisydx*xdu(9)
          dstran(2) = dstran(2) + dNidy*xdu(2) + dNicxdy*xdu(4)
     1    + dNisxdy*xdu(6) + dNicydy*xdu(8) + dNisydy*xdu(10)
          dstran(4) = dstran(4) + 
     1         dNidy*xdu(1) + dNicxdy*xdu(3) + dNisxdy*xdu(5) +
     2         dNidx*xdu(2) + dNicxdx*xdu(4) + dNisxdx*xdu(6) + 
     3         dNicydy*xdu(7) + dNicydx*xdu(8) + dNisydy*xdu(9) +
     4          dNisydx*xdu(10)

c        deformation gradient

          defGrad(1,1) = defGrad(1,1) + dNidx*utmp(1)+ dNicxdx*utmp(3)
     1    + dNisxdx*utmp(5) + dNicydx*utmp(7) + dNisydx*utmp(9)
          defGrad(1,2) = defGrad(1,2) + dNidy*utmp(1)+ dNicxdy*utmp(3)
     1    + dNisxdy*utmp(5) + dNicydy*utmp(7) + dNisydy*utmp(9) 
          defGrad(2,1) = defGrad(2,1) + dNidx*utmp(2)+ dNicxdx*utmp(4)
     1   + dNisxdx*utmp(6) + dNicydx*utmp(8) + dNisydx*utmp(10)
          defGrad(2,2) = defGrad(2,2) + dNidy*utmp(2) + dNicxdy*utmp(4)
     1    + dNisxdy*utmp(6)  + dNicydy*utmp(8)  + dNisydy*utmp(10)
          
        end do
        !WRITE(*,*)"defGrad*****************************************"   
        !  
        !  WRITE(*,*)(defGrad)
     
c        
c       CALL CONSTITUTIVE ROUTINE
c

        isvinc= (kintk-1)*nsvint  ! integration point increment
c
c       prepare arrays for entry into material routines
c
        do i = 1, nsvint
          statevLocal(i)=svars(i+isvinc)
        end do
c
c       state variables
c
!DEC$ NOVECTOR
        do k1=1,ntens
           stran(k1) = statevLocal(k1)
           stress(k1) = zero
        end do
c
        do i=1, ntens
!DEC$ NOVECTOR
          do j=1, ntens
            ddsdde(i,j) = zero
          end do
          ddsdde(i,j) = one
      enddo
c
c       compute characteristic element length
c
        celent = sqrt(djac*dble(ninpt))
        dvmat  = djac*thickness
        !      dvmat = zero
        !celent = one
      
        dvdv0 = one
        call material_lib_mech(materiallib,stress,ddsdde,
     1       stran,dstran,kintk,dvdv0,dvmat,defGrad,
     2       predef_loc,dpredef_loc,npredf,celent,coords_ip)
c
        !
        !WRITE(*,*)"ddsdde*****************************************"   
        ! 
        !  WRITE(*,*)(ddsdde)
        !
      
      !  WRITE(*,*)"STRESS*****************************************"   
      !     DO K1=1,4
      !    WRITE(*,*)(stress(K1))
      !ENDDO
           D(1,1)=ddsdde(1,1)
           D(2,2)=ddsdde(2,2)
           D(1,2)=ddsdde(1,2)
           D(2,1)=ddsdde(2,1)
           D(1,3)=ddsdde(1,4)
           D(2,3)=ddsdde(2,4)
           D(3,1)=ddsdde(4,1)
           D(3,2)=ddsdde(4,2)
           D(3,3)=ddsdde(4,4)
		   
		       STRESS1(1)=STRESS(1)
			   STRESS1(2)=STRESS(2)
			   STRESS1(3)=STRESS(4)
			 
        do k1=1,ntens
           statevLocal(k1) = stran(k1) + dstran(k1)
        end do
c
        isvinc= (kintk-1)*nsvint  ! integration point increment
c
c       update element state variables 
c
        do i = 1, nsvint
          svars(i+isvinc)=statevLocal(i)
      end do
c        
c       form stiffness matrix and internal force vector
c

      
        dNjdx = zero
        dNjdy = zero
        do i = 1, ndof*nnode
          force(i) = zero
          do j = 1, ndof*nnode
            stiff(j,i) = zero
          end do
        end do

        dvol= WEIGHT6(kintkx)*WEIGHT6(kintky)*djac
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------
         DO I=1,NDOFEL
            DO K=1,3
                 S(I,K)=0.0D0
            ENDDO
          ENDDO
            
               S=MATMUL(TRANSPOSE(B),D)		
			   CCBB = zero
               CCBB=MATMUL(S,B)
		       AM = MATMUL(TRANSPOSE(AN),AN)

          
          CMAS = CMAS + (AM * Rho * dvol)
c       assemble rhs and lhs
c

       DO K1=1,NDOFEL
            DO K2=1,3
       RHS(K1,1)=RHS(K1,1)-B(K2,K1)*STRESS1(K2)*dvol
            
            ENDDO
       ENDDO
      
       AMATRX=AMATRX+(CCBB*dvol)
C             write(*,*)('amatrx**************************************')
C              do k1=1,ndofel
C 			    write(*,*)(k1)
C              write(*,*)(amatrx(k1,k2),k2=1,ndofel)
C        enddo 
        end do
        end do 
c		Dynamic ----------------------------------------------------       
      IF (LFLAGS(1).EQ.11 .OR. LFLAGS(1).EQ.12) THEN
	     ALPHA = PARAMS(1)
          BETA  = PARAMS(2)
          GAMMA = PARAMS(3)
          
          DADU = ONE/(BETA*dtime**2)
          DVDU = GAMMA/(BETA*dtime)
!c IF Lumped MASS--------------------------------------------    
!      do k1 = 1, ndof*nnode
!        do  k2 = 1, ndof*nnode
!          AMASS(k1,k1)= CMAS(k1,k2) + AMASS(k1,k1)
!          enddo
!      enddo
!c IF CONSISTENT MASS--------------------------------------------
      do k1 = 1, ndofel
        do  k2 = 1, ndofel
          AMASS(k1,k2)= CMAS(k1,k2) 
        enddo
      enddo
	      amatrx = amatrx + AMASS * DADU
!          
      do k1=1, ndofel
          do k2=1, ndofel
       rhs(k1,1) = rhs(k1,1)- AMASS(k1,k2)*a(k2) 
          enddo
      end do
      endif
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------
C         do nodj = 1, nnode

C           incr_col = (nodj - 1)*ndof

C       !    write(*,*) 'nodj'
C       !WRITE(*,*)(nodj)
C       !
C c***********************************************************
C c******* change here this 
C       
C           dNjdx = bmat(1 + (nodj-1)*(ndim + 2*en_n) )
C           dNjdy = bmat(2 + (nodj-1)*(ndim + 2*en_n) )
C           dNjcxdx = bmat(3 + (nodj-1)*(ndim + 2*en_n))
C           dNjcxdy = bmat(4 + (nodj-1)*(ndim + 2*en_n))
C           dNjcydx = bmat(5 + (nodj-1)*(ndim + 2*en_n))
C           dNjcydy = bmat(6 + (nodj-1)*(ndim + 2*en_n))
C           
C           force_p(1) = dNjdx*stress(1) + dNjdy*stress(4)
C           force_p(2) = dNjdy*stress(2) + dNjdx*stress(4)
C           force_p(3) = dNjcxdx*stress(1) + dNjcxdy*stress(4)
C           force_p(4) = dNjcxdy*stress(2) + dNjcxdx*stress(4)
C           force_p(5) = dNjcydx*stress(1) + dNjcydy*stress(4)
C           force_p(6) = dNjcydy*stress(2) + dNjcydx*stress(4)
C           
C           do jdof = 1, ndof

C             jcol = jdof + incr_col

C             force(jcol) = force(jcol) +
C      1           force_p(jdof)*dvol
C             
C           end do
C           
C              !    write(*,*) 'force'
C              !WRITE(*,*)(force)

C           do nodi = 1, nnode

C             incr_row = (nodi -1)*ndof


C           dNidx = bmat(1 + (nodi-1)*(ndim + 2*en_n))
C           dNidy = bmat(2 + (nodi-1)*(ndim + 2*en_n))
C           dNicxdx = bmat(3 + (nodi-1)*(ndim + 2*en_n))
C           dNicxdy = bmat(4 + (nodi-1)*(ndim + 2*en_n))
C           dNicydx = bmat(5 + (nodi-1)*(ndim + 2*en_n))
C           dNicydy = bmat(6 + (nodi-1)*(ndim + 2*en_n))
C           
C            !IF (nodi .EQ. 2) THEN
C            !                write(*,*) "dNidx"    
C            !   write(*,*)(dNidx)
C            !   write(*,*) "dNicxdx"    
C            !   write(*,*)(dNicxdx)
C            !   write(*,*) "dNicxdy"    
C            !   write(*,*)(dNicxdy)
C            !   write(*,*)('ddsdde(1,1)')
C            !   write(*,*)(ddsdde(1,1))
C            !   
C            !   end if
C           
C           !IF (KINC .EQ. 20) THEN
C           !    WRITE(*,*) "stress"    
C           !    WRITE(*,*)(stress)
C           !    ENDIF
C c ONE =================================================================           
C             stiff_p(1,1) = dNidx*ddsdde(1,1)*dNjdx
C      1           + dNidy*ddsdde(4,4)*dNjdy
C      2           + dNidx*ddsdde(1,4)*dNjdy
C      3           + dNidy*ddsdde(4,1)*dNjdx
C              
C             stiff_p(1,2) = dNidx*ddsdde(1,2)*dNjdy
C      1           + dNidy*ddsdde(4,4)*dNjdx
C      2           + dNidx*ddsdde(1,4)*dNjdx
C      3           + dNidy*ddsdde(4,2)*dNjdy
C             
C            stiff_p(1,3) = dNidx*ddsdde(1,1)*dNjcxdx
C      1           + dNidy*ddsdde(4,4)*dNjcxdy
C      2           + dNidx*ddsdde(1,4)*dNjcxdy
C      3           + dNidy*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(1,4) = dNidx*ddsdde(1,2)*dNjcxdy
C      1           + dNidy*ddsdde(4,4)*dNjcxdx
C      2          + dNidx*ddsdde(1,4)*dNjcxdx
C      3          + dNidy*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(1,5) = dNidx*ddsdde(1,1)*dNjcydx
C      1           + dNidy*ddsdde(4,4)*dNjcydy
C      2           + dNidx*ddsdde(1,4)*dNjcydy
C      3           + dNidy*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(1,6) = dNidx*ddsdde(1,2)*dNjcydy
C      1                + dNidy*ddsdde(4,4)*dNjcydx
C      2           + dNidx*ddsdde(1,4)*dNjcydx
C      3           + dNidy*ddsdde(4,2)*dNjcydy
C c TWO =================================================================           
C             stiff_p(2,1) = dNidy*ddsdde(2,1)*dNjdx
C      1                + dNidx*ddsdde(4,4)*dNjdy
C      2           + dNidy*ddsdde(2,4)*dNjdy
C      3           + dNidx*ddsdde(4,1)*dNjdx

C             stiff_p(2,2) = dNidy*ddsdde(2,2)*dNjdy
C      1                + dNidx*ddsdde(4,4)*dNjdx
C      2           + dNidy*ddsdde(2,4)*dNjdx
C      3           + dNidx*ddsdde(4,2)*dNjdy
C             
C            stiff_p(2,3) = dNidy*ddsdde(2,1)*dNjcxdx
C      1                + dNidx*ddsdde(4,4)*dNjcxdy
C      2           + dNidy*ddsdde(2,4)*dNjcxdy
C      3           + dNidx*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(2,4) = dNidy*ddsdde(2,2)*dNjcxdy
C      1                + dNidx*ddsdde(4,4)*dNjcxdx
C      2           + dNidy*ddsdde(2,4)*dNjcxdx
C      3           + dNidx*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(2,5) = dNidy*ddsdde(2,1)*dNjcydx
C      1                + dNidx*ddsdde(4,4)*dNjcydy
C      2           + dNidy*ddsdde(2,4)*dNjcydy
C      3           + dNidx*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(2,6) = dNidy*ddsdde(2,2)*dNjcydy
C      1                + dNidx*ddsdde(4,4)*dNjcydx
C      2           + dNidy*ddsdde(2,4)*dNjcydx
C      3           + dNidx*ddsdde(4,2)*dNjcydy
C c THREE ===========================================================          
C             stiff_p(3,1) = dNicxdx*ddsdde(1,1)*dNjdx
C      1              + dNicxdy*ddsdde(4,4)*dNjdy
C      2           + dNicxdx*ddsdde(1,4)*dNjdy
C      3          + dNicxdy*ddsdde(4,1)*dNjdx

C             stiff_p(3,2) = dNicxdx*ddsdde(1,2)*dNjdy
C      1           + dNicxdy*ddsdde(4,4)*dNjdx
C      2           + dNicxdx*ddsdde(1,4)*dNjdx
C      3           + dNicxdy*ddsdde(4,2)*dNjdy
C             
C            stiff_p(3,3) = dNicxdx*ddsdde(1,1)*dNjcxdx
C      1           + dNicxdy*ddsdde(4,4)*dNjcxdy
C      2           + dNicxdx*ddsdde(1,4)*dNjcxdy
C      3           + dNicxdy*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(3,4) = dNicxdx*ddsdde(1,2)*dNjcxdy
C      1           + dNicxdy*ddsdde(4,4)*dNjcxdx
C      2           + dNicxdx*ddsdde(1,4)*dNjcxdx
C      3           + dNicxdy*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(3,5) = dNicxdx*ddsdde(1,1)*dNjcydx
C      1           + dNicxdy*ddsdde(4,4)*dNjcydy
C      2           + dNicxdx*ddsdde(1,4)*dNjcydy
C      3           + dNicxdy*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(3,6) = dNicxdx*ddsdde(1,2)*dNjcydy
C      1           + dNicxdy*ddsdde(4,4)*dNjcydx
C      2           + dNicxdx*ddsdde(1,4)*dNjcydx
C      3           + dNicxdy*ddsdde(4,2)*dNjcydy
C c FOUR =================================================================          
C             stiff_p(4,1) = dNicxdy*ddsdde(2,1)*dNjdx
C      1           + dNicxdx*ddsdde(4,4)*dNjdy
C      2           + dNicxdy*ddsdde(2,4)*dNjdy
C      3           + dNicxdx*ddsdde(4,1)*dNjdx

C             stiff_p(4,2) = dNicxdy*ddsdde(2,2)*dNjdy
C      1           + dNicxdx*ddsdde(4,4)*dNjdx
C      2           + dNicxdy*ddsdde(2,4)*dNjdx
C      3           + dNicxdx*ddsdde(4,2)*dNjdy
C             
C            stiff_p(4,3) = dNicxdy*ddsdde(2,1)*dNjcxdx
C      1           + dNicxdx*ddsdde(4,4)*dNjcxdy
C      2           + dNicxdy*ddsdde(2,4)*dNjcxdy
C      3           + dNicxdx*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(4,4) = dNicxdy*ddsdde(2,2)*dNjcxdy
C      1           + dNicxdx*ddsdde(4,4)*dNjcxdx
C      2           + dNicxdy*ddsdde(2,4)*dNjcxdx
C      3           + dNicxdx*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(4,5) = dNicxdy*ddsdde(2,1)*dNjcydx
C      1           + dNicxdx*ddsdde(4,4)*dNjcydy
C      2           + dNicxdy*ddsdde(2,4)*dNjcydy
C      3           + dNicxdx*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(4,6) = dNicxdy*ddsdde(2,2)*dNjcydy
C      1           + dNicxdx*ddsdde(4,4)*dNjcydx
C      2           + dNicxdy*ddsdde(2,4)*dNjcydx
C      3           + dNicxdx*ddsdde(4,2)*dNjcydy
C            
C c 5IVE =================================================================          
C             stiff_p(5,1) = dNicydx*ddsdde(1,1)*dNjdx
C      1           + dNicydy*ddsdde(4,4)*dNjdy
C      2           + dNicydx*ddsdde(1,4)*dNjdy
C      3           + dNicydy*ddsdde(4,1)*dNjdx

C             stiff_p(5,2) = dNicydx*ddsdde(1,2)*dNjdy
C      1           + dNicydy*ddsdde(4,4)*dNjdx
C      2           + dNicydx*ddsdde(1,4)*dNjdx
C      3           + dNicydy*ddsdde(4,2)*dNjdy
C             
C            stiff_p(5,3) = dNicydx*ddsdde(1,1)*dNjcxdx
C      1           + dNicydy*ddsdde(4,4)*dNjcxdy
C      2           + dNicydx*ddsdde(1,4)*dNjcxdy
C      3           + dNicydy*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(5,4) = dNicydx*ddsdde(1,2)*dNjcxdy
C      1           + dNicydy*ddsdde(4,4)*dNjcxdx
C      2           + dNicydx*ddsdde(1,4)*dNjcxdx
C      3           + dNicydy*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(5,5) = dNicydx*ddsdde(1,1)*dNjcydx
C      1           + dNicydy*ddsdde(4,4)*dNjcydy
C      2           + dNicydx*ddsdde(1,4)*dNjcydy
C      3           + dNicydy*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(5,6) = dNicydx*ddsdde(1,2)*dNjcydy
C      1           + dNicydy*ddsdde(4,4)*dNjcydx
C      2           + dNicydx*ddsdde(1,4)*dNjcydx
C      3           + dNicydy*ddsdde(4,2)*dNjcydy
C c 6IX =================================================================         
C              stiff_p(6,1) = dNicydy*ddsdde(2,1)*dNjdx
C      1           + dNicydx*ddsdde(4,4)*dNjdy
C      2           + dNicydy*ddsdde(2,4)*dNjdy
C      3           + dNicydx*ddsdde(4,1)*dNjdx

C             stiff_p(6,2) = dNicydy*ddsdde(2,2)*dNjdy
C      1           + dNicydx*ddsdde(4,4)*dNjdx
C      2           + dNicydy*ddsdde(2,4)*dNjdx
C      3           + dNicydx*ddsdde(4,2)*dNjdy
C             
C            stiff_p(6,3) = dNicydy*ddsdde(2,1)*dNjcxdx
C      1           + dNicydx*ddsdde(4,4)*dNjcxdy
C      2           + dNicydy*ddsdde(2,4)*dNjcxdy
C      3           + dNicydx*ddsdde(4,1)*dNjcxdx
C            
C          stiff_p(6,4) = dNicydy*ddsdde(2,2)*dNjcxdy
C      1           + dNicydx*ddsdde(4,4)*dNjcxdx
C      2           + dNicydy*ddsdde(2,4)*dNjcxdx
C      3           + dNicydx*ddsdde(4,2)*dNjcxdy
C          
C           stiff_p(6,5) = dNicydy*ddsdde(2,1)*dNjcydx
C      1           + dNicydx*ddsdde(4,4)*dNjcydy
C      2           + dNicydy*ddsdde(2,4)*dNjcydy
C      3           + dNicydx*ddsdde(4,1)*dNjcydx
C           
C            stiff_p(6,6) = dNicydy*ddsdde(2,2)*dNjcydy
C      1           + dNicydx*ddsdde(4,4)*dNjcydx
C      2           + dNicydy*ddsdde(2,4)*dNjcydx
C      3           + dNicydx*ddsdde(4,2)*dNjcydy
C             !IF (KINC .EQ. 1) THEN
C               !              WRITE(*,*) "stiff_p"    
C               !WRITE(*,*)(stiff_p)
C               !ENDIF
C             do jdof = 1, ndof
C               icol = jdof + incr_col
C               do idof = 1, ndof
C                 irow = idof + incr_row
C                 stiff(irow,icol) = stiff(irow,icol) +
C      1               stiff_p(idof,jdof)*dvol
C                 
C                 
C               end do
C             end do
C             
C             
C           end do
C       end do
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------
c------------------------------------------------------------------------    
                 
 
      
          
c --------------------------------------------------------------------   
C           ALPHA = PARAMS(1)
C           BETA  = PARAMS(2)
C           GAMMA = PARAMS(3)
C           
C           DADU = ONE/(BETA*dtime**2)
C           DVDU = GAMMA/(BETA*dtime)
C           
C           AM = MATMUL(TRANSPOSE(AN),AN)

C           
C           CMAS = CMAS + (AM * Rho * dvol)

C           
C c------------------------------------------------------------------------
C c
C c       assemble rhs and lhs
C c

C         do k1=1, ndof*nnode
C          rhs(K1,1)=rhs(K1,1)- force(k1)    
C           do k2=1, ndof*nnode
C             amatrx(k1, k2) = amatrx(k1, k2) + stiff(k1,k2)
C           end do
C       end do
C       

C              !               write(*,*) 'stress'
C              !WRITE(*,*)(stress)
C         end do
C       end do       ! end loop on material integration points
      !write(*,*)('amatrx*********************')
      !       DO K1=1,24 
      !          
      !      WRITE(*,*)(amatrx(K1,K2),K2=1,24)
      !ENDDO 
! Dynamic ----------------------------------------------------        
!      IF (LFLAGS(1).EQ.11 .OR. LFLAGS(1).EQ.12) THEN
!c IF Lumped MASS--------------------------------------------    
!      do k1 = 1, ndof*nnode
!        do  k2 = 1, ndof*nnode
!          AMASS(k1,k1)= CMAS(k1,k2) + AMASS(k1,k1)
!          enddo
!      enddo
!c IF CONSISTENT MASS--------------------------------------------
!      !      do k1 = 1, ndof*nnode
!      !  do  k2 = 1, ndof*nnode
!      !    AMASS(k1,k2)= CMAS(k1,k2) 
!      !    enddo
!      !enddo
!
!      amatrx = amatrx + AMASS * DADU
!          
!          do k1=1, ndof*nnode
!              do k2=1, ndof*nnode
!       rhs(k1,1) = rhs(k1,1)- AMASS(k1,k2)*a(k2) 
!          enddo
!      end do
!!      endif
      !           WRITE(*,*) 'stress'
      !!DO K1=1,4
      ! WRITE(*,*)(stress)
      !ENDDO 
      !WRITE(*,*) 'force*******************************'
      !DO K1=1,ndofel
      ! WRITE(*,*)(force(K1))
      !ENDDO
      !

      
      !      WRITE(*,*) 'RHS (Right-Hand-Side) Vector'
      !      DO K1=1,ndofel
      !WRITE(*,*) (rhs(K1,1))
      !ENDDO
      pnewdt = pnewdtLocal
      
      ! write(*,*)('amatrx**************************************')
      !       DO K1=1,24 
      !      WRITE(*,*)(amatrx(K1,K2),K2=1,24)
      !ENDDO
c
 999  continue
c
       if (jelem.eq.(9)) then
C               write(*,*) 'u(top-left///////////\\\\\\\\\'
C             if (kinc.eq.1) then
C             write(*,*) 'kinc:'
           write(*,*)(kinc)
C                end if
            do n=1, ndofel
            write(*,*)(u(n))
            end do
               end if 
      return

      end
