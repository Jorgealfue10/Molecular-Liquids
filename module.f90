module generator
        implicit none

        contains

        !Subroutine for calculating the SW potential of a single trial move of the MC simulation.
        !This subroutine takes the number of atoms, the randomly selected atom k, the initial geometry (ri),
        !the trial geometry (rf), the neighbour lists matrix (mat), the cell size (l) and the parameters for 
        !the SW potential calculation (A,B,lambda,gamma) as input variables. 
        !It returns the change in the potential after the trial move (dv).
        subroutine potcalc(N,k,ri,rf,mat,rc,l,A,B,lambda,gamma,dv)
                    implicit none
                    !Intrinsic variables of the subroutine.
                    integer :: i,j,m
                    !Angle formed by the atoms in the 3-body coefficients.
                    real(kind=8) :: angle
                    !Vector distances and module of the distance between the pair of atoms.
                    real(kind=8),dimension(3) :: drij,drim,drjm,drji
                    real(kind=8) :: modlij,modlim,modljm,modlji
                    real(kind=8),dimension(3) :: drijf,drimf,drjmf,drjif
                    real(kind=8) :: modlijf,modlimf,modljmf,modljif
                    !Coefficients for calculating the 2 and 3-body potential terms.
                    real(kind=8) :: f2,f2f,f3,f3f
                    real(kind=8) :: v2,v2f,v3,v3f
                    !Input and output variables.
                    integer,intent(in) :: N,k
                    integer,intent(in) :: mat(N,N)
                    real(kind=8),intent(in) :: ri(3,N),rf(3,N)
                    real(kind=8),intent(in) :: A,B,lambda,gamma,rc,l
                    real(kind=8),intent(inout) :: dv

                    !Resetting potential coefficients.
                    v2=0.
                    v2f=0.
                    v3=0.
                    v3f=0.
                    f2=0.
                    f2f=0.
                    f3=0.
                    f3f=0.

                    !Iterating over all the neighbour atoms of atom k.
                    j=1
                    do while (mat(k,j).ne.0)

                      !Computing the MI distance vector and modulus of both geometries.
                      call MIdist(ri(:,k),ri(:,mat(k,j)),l,modlij,drij)
                      call MIdist(rf(:,k),rf(:,mat(k,j)),l,modlijf,drijf)

                      !If the modulus distance between the atoms is lower than the cutoff,
                      !calculate the f2 2-body potential coefficient.
                      if (modlij.lt.rc**2.) then
                        f2=A*(B*(sqrt(modlij)**(-4.))-1)*exp(1./(sqrt(modlij)-rc))
                      else
                        f2=0.
                      endif

                      if (modlijf.lt.rc**2.) then
                        f2f=A*(B*(sqrt(modlijf)**(-4.))-1)*exp(1./(sqrt(modlijf)-rc))
                      else
                        f2f=0.
                      endif

                      !Storing the coefficients in the 2-body potential variables.
                      v2=v2+f2
                      v2f=v2f+f2f

                      !Iterating over the neighbour list to calculate the 3-body potential coefficient
                      !when considering the atom k at the center of the triatomic system.
                      m=j+1
                      do while (mat(k,m).ne.0)

                        !Computing the MI distance vector and modulus between k and the 3rd atom
                        !for both geometries.
                        call MIdist(ri(:,k),ri(:,mat(k,m)),l,modlim,drim)
                        call MIdist(rf(:,k),rf(:,mat(k,m)),l,modlimf,drimf)

                        !If the distance modulus between the atoms with atom k is lower than the cut-off
                        !distance, the 3-body potential coefficient is calculated, if not it is 0.
                        if ((modlij.ge.rc**2.).or.(modlim.ge.rc**2.)) then
                          f3=0.                            
                        else
                          angle=theta(drij,drim,modlij,modlim)
                          f3=lambda*exp((gamma/(sqrt(modlij)-rc))+(gamma/(sqrt(modlim)-rc)))*&
                          (cos(angle)+1./3.)**2.               
                        endif 

                        if ((modlijf.ge.rc**2.).or.(modlimf.ge.rc**2.)) then
                          f3f=0.                            
                        else
                          angle=theta(drijf,drimf,modlijf,modlimf)
                          f3f=lambda*exp((gamma/(sqrt(modlijf)-rc))+(gamma/(sqrt(modlimf)-rc)))*&
                          (cos(angle)+1./3.)**2. 
                        endif

                        !Storing the coefficients in the 3-body potential variables.
                        v3=v3+f3
                        v3f=v3f+f3f

                        !Incrementing index for the new third atom.
                        m=m+1

                      enddo

                      !Resetting the 3-body coefficients.
                      f3=0.
                      f3f=0.

                      !Computing the MI distance vector and modulus between k and the second atom
                      !for both geometries.
                      call MIdist(ri(:,mat(k,j)),ri(:,k),l,modlji,drji)
                      call MIdist(rf(:,mat(k,j)),rf(:,k),l,modljif,drjif)

                      !do i=1,3
                      !if (drij(i).ne.drji(i)) print*,"wut",drij(i),drji(i)
                      !if (drijf(i).ne.drjif(i)) print*,"wut-f",drijf(i),drjif(i)
                      !enddo

                      !Iterating over the neighbour list of atom mat(k,j) to calculate the 3-body potential coefficient
                      !when considering the atom mat(k,j) at the center of the triatomic system.
                      m=1
                      do while (mat(mat(k,j),m).ne.0)
                        !Checking that when iterating over the neighbout list of atom mat(k,j), 
                        !the atom k is not considered again.
                        if (mat(mat(k,j),m).ne.k) then 

                          !Computing the MI distance vector and modulus between mat(k,j) and the third atom
                          !for both geometries.
                          call MIdist(ri(:,mat(k,j)),ri(:,mat(mat(k,j),m)),l,modljm,drjm)
                          call MIdist(rf(:,mat(k,j)),rf(:,mat(mat(k,j),m)),l,modljmf,drjmf)

                          !If the distance modulus between the atoms with atom mat(k,j) is lower than the cut-off
                          !distance, the 3-body potential coefficient is calculated, if not it is 0.
                          if ((modlji.ge.rc**2.).or.(modljm.ge.rc**2.)) then
                            f3=0.                            
                          else
                            angle=theta(drji,drjm,modlji,modljm)
                            f3=lambda*exp((gamma/(sqrt(modlji)-rc))+(gamma/(sqrt(modljm)-rc)))*&
                            (cos(angle)+1./3.)**2. 
                          endif       
                                                
                          if ((modljif.ge.rc**2.).or.(modljmf.ge.rc**2.)) then
                            f3f=0.     
                          else
                            angle=theta(drjif,drjmf,modljif,modljmf)
                            f3f=lambda*exp((gamma/(sqrt(modljif)-rc))+(gamma/(sqrt(modljmf)-rc)))*&
                            (cos(angle)+1./3.)**2.
                          endif  

                          !Incrementing index for the new third atom.
                          m=m+1
                        else
                          !Incrementing index for the new third atom.
                          m=m+1
                          f3=0.
                          f3f=0.
                        endif

                        !Storing the coefficients in the 3-body potential variables.
                        v3=v3+f3
                        v3f=v3f+f3f

                    enddo
                    !Incrementing index for the new second atom.
                    j=j+1
                  enddo

                  !Computing the potential difference between the trial structure and the previous one.
                  dv=(v2f-v2)+(v3f-v3)

        end subroutine potcalc

        !Subroutine for calculating the neighbour matrix. It takes the number of atoms (N), 
        !a geometry matrix (r), the distance to consider neighbours (rneigh) and the cell size (l)
        !as input variables. It returns a matrix with N lines, each line corresponds to the atom i
        !and it contains the indexes of the atoms whose distance modulus with i is lower than rneigh.
        subroutine nlists(N,r,mat,nneigh,rneigh,l)
                    implicit none
                    !Intrinsic variables.
                    integer :: i,j,m
                    real(kind=8) :: modl,dr(3)
                    !Input and output variables.
                    integer,intent(in) :: N,nneigh
                    integer,intent(out) :: mat(N,nneigh)
                    real(kind=8),intent(in) :: r(3,N),l,rneigh

                    !Setting all neighbour lists as 0.
                    mat=0

                    !The modulus distance between each atom (i) with all of the atoms (j), except itself,
                    !is calculated. If the distance is lower than rneigh, it is considered as a neighbour
                    !and the atom number is stored in the i line m column of the mat matrix.
                    do i=1,N
                      m=1
                      do j=1,N
                        call MIdist(r(:,i),r(:,j),l,modl,dr)
                        if ((modl.lt.rneigh**2.).and.(i.ne.j)) then
                          mat(i,m)=j
                          m=m+1
                        endif
                      enddo
                    enddo

        end subroutine nlists

        !Subroutine for updating the histogram containing the information for the RDF. It takes the number
        !of atoms (N), the number of distances (nr), the previous number of times Ni has been updated (nacc),
        !the previour distance histogram (Ni), the vector containing the distances for RDF (gdist), the distance
        !step (dgr) and the cell size (l) as input variables. The output is the updated Ni and the updated nacc.
        subroutine countNi(N,nr,nacc,r,Ni,gdist,dgr,l)
                    implicit none
                    !Intrinsic variables.
                    integer :: i,j,m
                    real(kind=8) :: modl,dr(3)
                    !Input and output variables.
                    integer,intent(in) :: N,nr
                    integer,intent(inout) :: Ni(nr),nacc
                    real(kind=8),intent(in) :: r(3,N),gdist(nr),dgr,l

                    !In this subroutine all unique distance modulus between atoms are calculated.
                    !If the calculated distance is inside a range between the distance gdist(m) and
                    !that distance plus dgr, the corresponding element for that distance of the 
                    !histogram Ni(m) will increase by one unit. Lastly, nacc will increment also by 1.
                    do m=1,nr
                        do i=1,N-1
                            do j=i+1,N
                                call MIdist(r(:,i),r(:,j),l,modl,dr)
                                if ((modl.ge.gdist(m)**2.).and.(modl.lt.(gdist(m)+dgr)**2.)) then
                                    Ni(m)=Ni(m)+1
                                endif
                            enddo
                        enddo
                    enddo

                    nacc=nacc+1

        end subroutine countNi

        !Subroutine for calculating the MI distance vector and squared modulus between two atoms. It takes the 
        !position of the first atom (r1), the position of the second atom (r2) and the cell size (l) as
        !input variables. The MI squared distance modulus (modl) and the MI distance vector (dr) are the output.
        subroutine MIdist(r1,r2,l,modl,dr)
                    implicit none
                    !Intrinsic variables.
                    integer :: i
                    !Input and output variables.
                    real(kind=8),intent(in) :: r1(3),r2(3),l
                    real(kind=8),intent(out) :: modl,dr(3)

                    !The vector distance between the two atoms is calculated. Then, coordinate by coordinate
                    !the values of this vector are check. If the distance coordinate is greater than the half
                    !of the box length, the atoms are really far apart in the same cell. So, depending on the
                    !sign of this distance coordinate the length of the unit cell is substracted or added.

                    dr=r1-r2

                    do i=1,3
                        if (dr(i).lt.(-l*0.5)) then
                          dr(i)=dr(i)+l
                        else if (dr(i).gt.(l*0.5)) then
                          dr(i)=dr(i)-l
                        endif
                    enddo

                    !Finally, the squared modulus distance is calculated.
                    modl=dr(1)**2.+dr(2)**2.+dr(3)**2.

        end subroutine MIdist

        !Subroutine for calculating the RDF (gr). It takes the number of atoms (N), the number of distances (nr),
        !value of pi, the reduced density (reddens), the histogram Ni, the step in distance (dgr) and the vector
        !containing the distances (gdist) as input variables. The output is the vector containing the RDF values.
        subroutine RDF(N,nr,nacc,pi,reddens,Ni,dgr,gdist,gr)
                    implicit none
                    !Intrinsic variables.
                    integer :: i
                    real(kind=8) :: vol
                    !Input and output variables.
                    integer,intent(in) :: N,nr,nacc
                    integer,intent(in) :: Ni(nr)
                    real(kind=8),intent(in) :: pi,reddens,dgr
                    real(kind=8),intent(in) :: gdist(nr)
                    real(kind=8),intent(out) :: gr(nr)

                    !For all the distances in the vector gdist, using the corresponding formula, the 
                    !RDF function is calculated.
                    do i=1,nr
                      vol=(gdist(i)+dgr)**3.-gdist(i)**3.
                      gr(i)=(3.*real(Ni(i)))/(real(nacc)*real(N)*0.5*4.*pi*vol*reddens)
                    enddo

        end subroutine RDF

        !Function for calculating the angle between three atoms. It takes the MI vector distances
        !and the MI squared distance modulus as input. It returns the value for the angle in radians.
        real(kind=8) function theta(dr1,dr2,modl1,modl2)
                    implicit none
                    real(kind=8),intent(in) :: dr1(3),dr2(3),modl1,modl2

                    theta=acos((dot_product(dr1,dr2))/(sqrt(modl1)*sqrt(modl2)))

      end function theta

end module generator