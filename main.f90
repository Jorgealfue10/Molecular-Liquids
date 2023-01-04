program mainSW
        !use modSW
        use generator

        implicit none
        !Program parameters.
        real(kind=8),parameter :: kb=8.314462618d-3,pi=3.14159265359
        !General variables.
        integer :: i,k,p
        real(kind=8) :: l !Size of unit cell
        character(len=25) :: text

        !Montecarlo and SW potential parameters and variables.
        !Parameters.
        !real(kind=8),parameter :: A=7.049556277,B=0.6022245584,lambda=21.0,gamma=1.20
        real(kind=8) :: A,B,lambda,gamma
        !Integer variables.
        integer :: N,steps,counter,countstruct,equil,prod,sweep
        !Real type variables.
        real(kind=8) :: T,w,rc,sigma,pr,at
        !Potential variables.
        real(kind=8) :: dv

        !Module distance variables.
        real(kind=8) :: modlref,drmax
        !Distance vectors and random changes.
        real(kind=8),dimension(3) :: drref,wpos,rprint
        !Position matrices.
        real(kind=8),allocatable,dimension(:,:) :: ri,rf,rref
        !Atomic symbols.
        character(len=2) :: symbol

        !SW Neighbour list not yet implemented.
        !Approx. number of neighbours.
        integer :: nneigh
        !Matrix.
        integer,allocatable :: mat(:,:)
        !Max distances.
        real(kind=8) :: rskin,rneigh

        !RDF parameters and variables.
        !Number of distances and number of times the histogram is updated.
        integer :: nacc,nr
        !Histogram.
        integer,allocatable,dimension(:) :: Ni
        !Density, reduced density values and volume.
        real(kind=8) :: reddens
        !Initial and final RDF distances and step for RDF.
        real(kind=8) :: gro,grf,dgr
        !Radial distribution function and list of distances.
        real(kind=8),allocatable,dimension(:) :: gr,gdist

        !To obtain the proportion of accepted moves in the MC simulation.
        integer :: succes,fail
        !CPU Time related values. 
        real(kind=8) :: TMC1,TMC2,T1,T2

        call cpu_time(T1)

        !Opening input parameters file.
        open(15,file="Input/parameters.inp",action="read")
        !Opening initial geometry input file.
        open(16,file="Input/init-geom.xyz",action="read")
        !Open xyz output file.
        open(17,file="Output/prod_traj.xyz",action="write")
        !Open RDF output file.
        open(18,file="Output/gfunction.dat",action="write")
        !Open general output file.
        open(19,file="Output/SWM-data.out",action="write")

        !Reading input parameters.
        read(15,*)
        read(15,*)
        read(15,*) text, equil               !Number of equilibration sweeps
        read(15,*) text, prod                !Number of production sweeps
        read(15,*) text, l                   !Box length
        read(15,*) text, sigma               !Sigma value
        read(15,*) text, T                   !Temperature
        read(15,*) text, reddens             !Density
        read(15,*) text, rc                  !Cutoff distance value
        read(15,*) text, rskin               !Skin distance
        read(15,*) text, drmax               !Maximum distance in each Montecarlo step
        read(15,*) text, A,B                 !Parameters for calculating 2-body SW potential
        read(15,*) text, lambda,gamma        !Parameters for calculating 3-body SW potential
        read(15,*) text, gro, grf, dgr       !RDF initial, final and increment distances
        
        close(15)

        !Reading the number of atoms.
        read(16,*) N
        read(16,*)

        !Defining number of distances for the RDF.
        nr=nint((grf-gro)/dgr)

        !Allocating position matrices and neighbour lists matrix.
        allocate(ri(3,N),rf(3,N),rref(3,N))!,mat(N,N))
        !Allocating histogram, RDF distances vector and RDF function.
        allocate(Ni(nr),gdist(nr),gr(nr))

        !Reading the initial atoms positions.
        do i=1,N
                read(16,*) symbol, ri(:,i)
        enddo

        close(16)

        ri=ri/sigma

        !Generating distance vector for RDF.
        do i=0,nr
                gdist(i)=gro+dgr*i
        enddo

        !Total number of steps.
        steps=(equil+prod)*N 

        !Defining the neighbour distance.
        rneigh=rskin+rc
        !Defining an approximate number of neighbours, based on the density and the volume
        !formed by the neighbour distance.
        nneigh=int(reddens*(4./3.)*pi*rneigh**3.+10.)
        !Allocating the Neighbour lists matrix.
        allocate(mat(N,nneigh))
        mat=0
        !Generating Neighbour lists.
        call nlists(N,ri,mat,nneigh,rneigh,l)
        !Storing reference geometry.
        rref=ri

        !Setting counters to initial value.
        sweep=0
        counter=0
        countstruct=0
        Ni=0
        nacc=0
        succes=0
        fail=0

        !Writing the head of the ouput.
        write(19,*) "--------------------------------------------------------------------"
        write(19,*) "        Stiliger-Weber Model for Silicon Molecular Liquid"
        write(19,*) "--------------------------------------------------------------------"
        write(19,*)
        write(19,*) "Based on the Manual provided by Jeremy Harvey and the work by "
        write(19,*) "F. H. Stillinger and T. A. Weber [1]."
        write(19,*)
        write(19,*) "[1] F. H. Stillinger and T. A. Weber, Phys.Rev. B, 1985, 31, 5626 - 5271"
        write(19,*)
        write(19,*) "--------------------------------------------------------------------"
        write(19,*) "Summary of the input"
        write(19,*) "--------------------------------------------------------------------"
        write(19,*)
        write(19,*) "----------------General parameters of the simulation----------------"
        write(19,*)
        write(19,fmt='(A35,1x,I6)') "Number of atoms in the simulation: ", N
        write(19,fmt='(A31,1x,I10)') "Number of equilibration sweeps: ", equil
        write(19,fmt='(A29,3x,I10)') "Number of production sweeps: ", prod
        write(19,fmt='(A34,1x,I20)') "Total number of Montecarlo Steps: ", steps
        write(19,fmt='(A47,1x,F16.10)') "Unit cubic cell length (Reduced Units (r.u.)): ", l
        write(19,fmt='(A22,1x,F8.6,1x,A24,F8.6)') "Temperature (T) (ru): ", T, "; Density (rho) (r.u.): ", reddens
        write(19,*)
        write(19,*) "----------------Potential and Si related parameters-----------------"
        write(19,*)
        write(19,fmt='(A32,1x,F8.6)') "Sigma value for Si (Angstroms): ", sigma
        write(19,fmt='(A23,1x,F8.6,1x,A44,1x,F6.4)') "Cut-off distance (ru): ", rc, &
                                                     "; Skin distance for Neighbour lists (ru): ", rskin
        write(19,fmt='(A50,1x,F8.6)') "Maximum distance for considering neighbour atoms: ", rneigh
        write(19,fmt='(A65,1x,F8.6)') "Maximum distance allowed for random displacement of each atom (ru): ", drmax
        write(19,fmt='(A45,1x,F16.14,1x,A5,1x,F16.14)') "Parameters for 2-body SW potential terms. A: ", A, "; B: ", B
        write(19,fmt='(A50,1x,F12.8,1x,A9,1x,F10.8)') "Parameters for 3-body SW potential terms. lambda: ", lambda, &
                                                      "; gamma: ", gamma
        write(19,*)
        write(19,*) "---------------Radial Distribution function parameters---------------"
        write(19,*)
        write(19,fmt='(A23,1x,F10.8,1x,A21,1x,F10.8)') "Initial distance (ru): ",gro, "Final distance (ru): ",grf
        write(19,fmt='(A20,1x,F6.4,1x,A24,1x,I6)') "Distance step (ru): ",dgr, "Total number of points: ", nr
        write(19,*)
        write(19,*) "--------------------------------------------------------------------"
        write(19,*)

        call cpu_time(TMC1)

        write(19,*) "--------------------------------------------------------------------"
        write(19,*) "                 STARTING MONTECARLO SIMULATION"
        write(19,*)

        do p=1,steps                

                dv=0.

                !Redifining final positions matrix in each step.
                rf=ri

                !Increasing counter.
                counter=counter+1

                !Randomly determining the atom that is going to be displaced.
                call random_number(at)
                k=nint(at*N)
                do while (k.eq.0)                 !Repeating the process until k is not 0.
                        call random_number(at)
                        k=nint(at*N)
                enddo

                !Calling random number to the random displacement.
                call random_number(wpos)
                !Modifying the position of the k atom.
                rf(:,k)=rf(:,k)+2.*drmax*(wpos(:)-0.5)

                !Checking that the displaced atom is kept inside the unit cell boundaries.
                do i=1,3
                        if (rf(i,k).le.0.) then
                                rf(i,k)=rf(i,k)+l
                        elseif (rf(i,k).gt.l) then
                                rf(i,k)=rf(i,k)-l
                        endif
                enddo

                !Calculating the SW potential.
                call potcalc(N,k,ri,rf,mat,rc,l,A,B,lambda,gamma,dv)

                !For output pursposes.
                rprint(:)=ri(:,k)

                !Checking if the trial move is going to be accepted or not, 
                !based on the potential change.
                if (dv.le.0.) then
                        ri=rf
                        succes=succes+1
                else
                        call random_number(w)
                        pr=exp(-dV/T)
                        if (pr.ge.w) then
                                ri=rf
                                succes=succes+1
                        else
                                ri=ri
                                fail=fail+1
                        endif
                endif

                !Checking if the displacement of the k atom is big enough to remake the neighbour lists.
                !modlref=distMI(rref(:,k),ri(:,k),l)
                call MIdist(rref(:,k),ri(:,k),l,modlref,drref)
                if (modlref.gt.((rskin/2.)**2.)) then
                        call nlists(N,ri,mat,nneigh,rneigh,l)
                        rref(:,k)=ri(:,k)
                endif
                        

                !For each N steps, the sweep counter is increased. 
                !If the number of sweeps is higher than the equilibration sweeps of the input
                !the structure of the molecular liquid is printed and the RDF Ni histogram is updated.
                if (counter.eq.N) then
                        sweep=sweep+1
                        if (sweep.gt.equil) then
                                if (countstruct.eq.2) then
                                        write(17,*) N
                                        write(17,*) "Structure: ", sweep
                                        do i=1,N
                                                write(17,*) symbol,ri(:,i)
                                        enddo
                                        countstruct=0
                                endif
                                countstruct=countstruct+1

                                call countNi(N,nr,nacc,ri,Ni,gdist,dgr,l)

                                write(19,*) "-----------------------------------------------------------------"
                                write(19,*) "Moved atom: ", k
                                write(19,fmt='(A24,1x,3(F8.6,1x))') "Original geometry (ru): ", rprint(:)
                                write(19,fmt='(A21,1x,3(F8.6,1x))') "Trial geometry (ru): ", rf(:,k)
                                write(19,*) "       Geometry stored in file prod_traj.xyz"
                                write(19,fmt='(A14,1x,I12)') "Current step: ", p
                                write(19,fmt='(A26,1x,I10,1x,A24,I10)') "Number of accepted moves: ", succes,&
                                                                        "Number of failed moves: ", fail
                                write(19,fmt='(A28,1x,F16.12)') "Variation of potential (ru): ", dv
                                write(19,*)

                        endif
                        counter=0
                endif

                !print*,p
                
        enddo

        call cpu_time(TMC2)

        close(17)

        !The RDF histogram is divided by the number of times it has been updated.
        !Ni=Ni/nacc

        !Calculation of the RDF.
        call RDF(N,nr,nacc,pi,reddens,Ni,dgr,gdist,gr)

        do i=1,nr
                write(18,*) gdist(i),gr(i)
        enddo

        close(18)

        call cpu_time(T2)

        write(19,*) "--------------------------------------------------------------------"
        write(19,*) "                      END OF SIMULATION"
        write(19,*) "--------------------------------------------------------------------"
        write(19,fmt='(A32,1x,I10)') "Total number of accepted moves: ", succes
        write(19,fmt='(A38,1x,F5.2)') "Proportion of accepted trial moves (%)", (real(succes)/real(steps))*100.
        write(19,*) "--------------------------------------------------------------------"
        write(19,*)
        write(19,fmt='(A20,1x,F12.6)') "Total CPU time (s): ", T2-T1
        write(19,fmt='(A28,1x,F12.6)') "MC simulation CPU time (s): ", TMC2-TMC1
        write(19,*)
        write(19,*)
        write(19,*)
        write(19,*) "--------------------------------------------------------------------"
        write(19,*) "Program done by: Jorge Alonso de la Fuente"

        close(19)

stop
end program mainSW