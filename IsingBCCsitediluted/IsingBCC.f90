program IsingBCC
    implicit none
    save
!teste
    !decalrando variaveis
    intrinsic random_seed, random_number

    integer ::MCx, MCc
    integer ::passo
    integer ::magnetizacao, eneJ1, eneJ2
    double precision energia
    double precision  somaMag, somaMag2, somaMag4
    double precision  somaEne, somaEne2, somaEne4
    double precision  mediaMag,  mediaMag2, mediaMag4
    double precision  mediaEne,mediaEne2, mediaEne4
    double precision  susceptibilidade, cumuM
    double precision  calor
    real :: tempoInicial, tempoFinal, tempoCPU
    integer, dimension(:,:,:,:), allocatable ::sigma
    byte, dimension(:), allocatable ::ant, suc
    byte, dimension(:,:,:,:), allocatable ::bond_i,bond_j,bond_k
    byte :: histograma
    byte ::L
    double precision::t0, p   !p=0 sistema puro
    double precision::tin,tfi,dt
    real::J2,J1
    real::rando
    integer::spin
    integer::amostras
    double precision ::w(-8:8,-6:6)
    double precision ::PJ1,PJ2


    !Abrindo arquivos
    open(22,file='hist.dat')
    open(12,file='dados.dat')

    !inicilizando das variaveis
    !    call leDados
    call leargumentos
    J1=1-p
    allocate(ant(1:L))
    allocate(suc(1:L))
    allocate(bond_i(1:2,1:L,1:L,1:L))
    allocate(bond_j(1:2,1:L,1:L,1:L))
    allocate(bond_k(1:2,1:L,1:L,1:L))
    allocate(sigma(1:2,1:L,1:L,1:L))
    call cpu_time(tempoInicial)
    call initRandomSeed()
    call iniciaContorno
    call iniciaSigma
    call iniciaBond
    call dilui
    call geraLigacao
    if (histograma==1) then
        t0 = tin
        call atualiza
        do amostras = 1 , 50
        call iniciaSigma
        call iniciaBond
        call dilui
        call geraLigacao
        call iniciaVariaveis
        do passo = 1 , MCx
            call metropolis
            call wolff
        end do
        do passo = 1 , MCc
            call metropolis
            call wolff
            call calcularMagEng
            write(*,*) eneJ1, eneJ2,  magnetizacao
        end do
        end do
    else
        write(*,*) '#iniciando simulação com :'
        write(*,*) '# L =', L
        write(*,*) '# p =', p
        write(*,*) '# Tin =', Tin
        write(*,*) '# Tfi =', Tfi
        write(*,*) '# dt =', dt
        write(*,*) '# J1 =', J1
        write(*,*) '# J2 =', J2
        write(*,*) '# MCx =', MCx
        write(*,*) '# MCc =', MCc
        write(*,*) '# histograma =', histograma
        call iniciaSigma
        call iniciaBond
        call dilui
        call geraLigacao
        t0 = tin
        do while (t0<=tfi)
            call iniciaVariaveis
            call atualiza
            !repetiÃ§Ãµes termalizaÃ§Ã£o
            do passo = 1 , MCx
                call metropolis
                call wolff
            end do
            do passo = 1 , MCc
                call metropolis
                call wolff
                call calcularMagEng
                call calcularSoma
            end do
            call cpu_time(tempoFinal)
            TempoCPU=tempoFinal-tempoInicial
            call calcularMedia
            t0=t0+dt
        end do
    end if
    !fechando arquivos
    close(10)
    close(20)
    deallocate(ant)
    deallocate(suc)
    deallocate(bond_i)
    deallocate(bond_j)
    deallocate(bond_k)
    deallocate(sigma)
  !  write(*,*) '# fim do programa'
   !fim do programa

CONTAINS
    !-----------------------------------------------------------------------------
    subroutine calcularSoma
        double precision ::mag, mag2, mag4, ene,ene2, ene4
        mag=float(magnetizacao)/(2*L**3)
        mag2=mag*mag
        mag4=mag2*mag2
        ene =  energia
        ene2 = ene*ene
        ene4 = ene2*ene2
        somaMag=somaMag+mag
        somaMag2=somaMag2+mag2
        somaMag4=somaMag4+mag4
        somaEne=somaEne+ene
        somaEne2=somaEne2+ene2
        somaEne4=somaEne4+ene4

    end subroutine calcularSoma
    !-----------------------------------------------------------------------------
    subroutine iniciaVariaveis

        somaMag=0
        somaMag2=0
        somaMag4=0
        somaEne=0
        somaEne2=0
        somaEne4=0

    end subroutine iniciaVariaveis

     !-----------------------------------------------------------------------------
    subroutine calcularMedia
        integer  ::numeroSitios

        numeroSitios = 2*L**3
        mediaMag = somaMag/MCc
        mediaMag2 = somaMag2/MCc
        mediaMag4 = somaMag4/MCc
        mediaEne=somaEne/MCc
        mediaEne2=somaEne2/MCc
        mediaEne4=somaEne4/MCc
        calor= (MediaEne2-mediaEne*mediaEne)/NumeroSitios/t0/t0
        susceptibilidade= (mediaMag2 - MediaMag*MediaMag)*NumeroSitios/t0
        cumuM=1-mediaMag4/(3*mediaMag2*mediaMag2)  ! rever essa equaÃ§Ã£o soma ou media
        write(*,*)t0, susceptibilidade, calor, cumuM, mediaMag, mediaEne, TempoCPU
    end subroutine calcularMedia
     !-----------------------------------------------------------------------------
    subroutine leDados
        read(12,*) L
        read(12,*) p
        read(12,*) tin
        read(12,*) tfi
        read(12,*) dt
        read(12,*) J1
        read(12,*) J2
        read(12,*) MCx
        read(12,*) MCc
        read(12,*) histograma
    end subroutine leDados
    !-----------------------------------------------------------------------------
    subroutine leargumentos
        character(len=32) :: arg
        Integer, Parameter :: wp = Selected_real_kind( 12, 70 )
        Real( wp ) :: a
        call get_command_argument(1, arg)
        read (arg,*) L
        call get_command_argument(2, arg)
        read(arg,*) a
        p=a;
        call get_command_argument(3, arg)
        read(arg,*) a
        !read(arg,*) tin
        tin=a
        call get_command_argument(4, arg)
        !read(arg,*) tfi
        read(arg,*) a
        tfi=a
        call get_command_argument(5, arg)
        !read(arg,*) dt
        read(arg,*) a
        dt=a
        call get_command_argument(6, arg)
        !read(arg,*) J1
        read(arg,*) a
        J1=a
        call get_command_argument(7, arg)
        read(arg,*) a
        J2=a
        !read(arg,*) J2
        call get_command_argument(8, arg)
        read(arg,*) MCx
        call get_command_argument(9, arg)
        read(arg,*) MCc
        call get_command_argument(10, arg)
        read(arg,*) histograma
    end subroutine leargumentos
    !-----------------------------------------------------------------------------
    subroutine iniciaContorno
        byte :: i
        do i = 1, L
            ant(i) = i - 1
            suc(i) = i + 1
        end do
        ant(1) = L
        suc(L) = 1
    end subroutine iniciaContorno

    !-----------------------------------------------------------------------------
    subroutine iniciaSigma
        sigma=1
    end subroutine iniciaSigma
    !-----------------------------------------------------------------------------
    subroutine iniciaBond
        bond_i=0
        bond_j=0
        bond_k=0
    end subroutine iniciaBond
    !-----------------------------------------------------------------------------
    subroutine  dilui
        integer :: n
        byte    ::  i,j,k,subrede
        real    ::aux
        n=0
        do while(n<int(2*p*L*L*L))
            call random_number(rando)
            i=int(L*rando)+1
            call random_number(rando)
            j=int(L*rando)+1
            call random_number(rando)
            k=int(L*rando)+1
            call random_number(rando)
            subrede=int(2*rando)+1
            if (abs(sigma(subrede, i,j,k))==1) then
                sigma(subrede, i,j,k)=0
                n=n+1
            end if
        end do
    end subroutine dilui

    !-----------------------------------------------------------------------------

    subroutine wolff
        byte :: i,j,k, subrede
        !sorteia sitio semente
        do
            call random_number(rando)
            i=int(L*rando)+1
            call random_number(rando)
            j=int(L*rando)+1
            call random_number(rando)
            k=int(L*rando)+1
            call random_number(rando)
            subrede=int(2*rando)+1
            if (abs(sigma(subrede, i,j,k))==1) exit
        end do
        spin= sigma(subrede, i,j,k)
        sigma(subrede, i,j,k)=-sigma(subrede, i,j,k)
        call perimetro(subrede, i,j,k)
    end subroutine wolff
       !-----------------------------------------------------------------------------
    recursive subroutine perimetro(subrede, i,j,k)
        byte, intent(in) :: i,j,k, subrede
        byte :: sub
        if(subrede==1) then
            sub=2
            if (sigma(sub,i,j,k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,j,k)=-spin
                    call perimetro(sub,i,j,k)
                end if
            end if
            !
            if (sigma(sub,suc(i),j,k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,suc(i),j,k)=-spin
                    call perimetro(sub,suc(i),j,k)
                end if
            end if
            !
            if (sigma(sub,suc(i),suc(j),k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,suc(i),suc(j),k)=-spin
                    call perimetro(sub,suc(i),suc(j),k)
                end if
            end if
            !
            if (sigma(sub,i,suc(j),k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,suc(j),k)=-spin
                    call perimetro(sub,i,suc(j),k)
                end if
            end if
            !
            if (sigma(sub,i,j,suc(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,j,suc(k))=-spin
                    call perimetro(sub,i,j,suc(k))
                end if
            end if
            !
            if (sigma(sub,suc(i),j,suc(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,suc(i),j,suc(k))=-spin
                    call perimetro(sub,suc(i),j,suc(k))
                end if
            end if
            !
            if (sigma(sub,suc(i),suc(j),suc(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,suc(i),suc(j),suc(k))=-spin
                    call perimetro(sub,suc(i),suc(j),suc(k))
                end if
            end if
            !
            if (sigma(sub,i,suc(j),suc(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,suc(j),suc(k))=-spin
                    call perimetro(sub,i,suc(j),suc(k))
                end if
            end if
        else
            !--------------------------------------------------------------------------
            sub=1
            if (sigma(sub,ant(i),ant(j),k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,ant(i),ant(j),k)=-spin
                    call perimetro(sub,ant(i),ant(j),k)
                end if
            end if
            !
            if (sigma(sub,i,ant(j),k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,ant(j),k)=-spin
                    call perimetro(sub,i,ant(j),k)
                end if
            end if
            !
            if (sigma(sub,i,j,k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,j,k)=-spin
                    call perimetro(sub,i,j,k)
                end if
            end if
            !
            if (sigma(sub,ant(i),j,k)==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,ant(i),j,k)=-spin
                    call perimetro(sub,ant(i),j,k)
                end if
            end if
            !
            if (sigma(sub,ant(i),ant(j),ant(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,ant(i),ant(j),ant(k))=-spin
                    call perimetro(sub,ant(i),ant(j),ant(k))
                end if
            end if
            !
            if (sigma(sub,i,ant(j),ant(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,ant(j),ant(k))=-spin
                    call perimetro(sub,i,ant(j),ant(k))
                end if
            end if
            !
            if (sigma(sub,i,j,ant(k))==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,i,j,ant(k))=-spin
                    call perimetro(sub,i,j,ant(k))
                end if
            end if
            !
            if (sigma(sub,ant(i),j,ant(k)) ==spin)then
                call random_number(rando)
                if (rando<PJ1) then
                    sigma(sub,ant(i),j,ant(k))=-spin
                    call perimetro(sub,ant(i),j,ant(k))
                end if
            end if
        end if
        !--------------------------------------------------------------------------
        sub=subrede;
        if (sigma(sub,ant(i),j,k)*Bond_i(sub,ant(i),j,k)==spin)then
            call random_number(rando)
            if (rando<PJ2) then
                sigma(sub,ant(i),j,k)=-spin
                call perimetro(sub,ant(i),j,k)
            end if
        end if
        !
        if (sigma(sub,suc(i),j,k)*Bond_i(sub,i,j,k)==spin)then
            call random_number(rando)
            if (rando<PJ2) then
                sigma(sub,suc(i),j,k)=-spin
                call perimetro(sub,suc(i),j,k)
            end if
        end if
        !
        if (sigma(sub,i,ant(j),k)*Bond_j(sub,i,ant(j),k)==spin)then
            call random_number(rando)
            if (rando<PJ2) then
                sigma(sub,i,ant(j),k)=-spin
                call perimetro(sub,i,ant(j),k)
            end if
        end if
        !
        if (sigma(sub,i,suc(j),k)*Bond_j(sub,i,j,k)==spin)then
            call random_number(rando)
            if (rando<PJ2) then
                sigma(sub,i,suc(j),k)=-spin
                call perimetro(sub,i,suc(j),k)
            end if
        end if
        !
        if (sigma(sub,i,j,ant(k))*Bond_k(sub,i,j,ant(k))==spin)then
            call random_number(rando)
            if (rando<PJ2) then
                sigma(sub,i,j,ant(k))=-spin
                call perimetro(sub,i,j,ant(k))
            end if
        end if
        if (sigma(sub,i,j,suc(k))*Bond_k(sub,i,j,k)==spin)then
            call random_number(rando)
            if (rando<PJ1) then
                sigma(sub,i,j,suc(k))=-spin
                call perimetro(sub,i,j,suc(k))
            end if
        end if
        return
    end subroutine perimetro

    !-----------------------------------------------------------------------------

    subroutine geraLigacao
        byte::  i,j,k

        do i=1,L
            do j=1,L
                do k=1,L
                    if( isolada(1,i,j,k) )then
                        bond_i(2,i,j,k)=1
                        bond_i(2,i,suc(j),k)=1
                        bond_i(2,i,j,suc(k))=1
                        bond_i(2,i,suc(j),suc(k))=1

                        bond_j(2,i,j,k)=1
                        bond_j(2,suc(i),j,k)=1
                        bond_j(2,i,j,suc(k))=1
                        bond_j(2,suc(i),j,suc(k))=1

                        bond_k(2,i,j,k)= 1
                        bond_k(2,suc(i),j,k)= 1
                        bond_k(2,i,suc(j),k)= 1
                        bond_k(2,suc(i),suc(j),k)= 1
                    end if
                end do
            end do
        end do

        do i=1,L
            do j=1,L
                do k=1,L
                    if( isolada(2,i,j,k) )then

                        Bond_i(1,ant(i),ant(j),ant(k))=1
                        Bond_i(1,ant(i),ant(j), k)=1
                        Bond_i(1,ant(i), j, k)=1
                        Bond_i(1,ant(i),j, ant(k))=1

                        Bond_j(1,ant(i), ant(j),ant(k))=1
                        Bond_j(1,ant(i), ant(j), k)=1
                        Bond_j(1, i,  ant(j), k)=1
                        Bond_j(1, i,  ant(j),ant(k))=1

                        Bond_k(1,ant(i),ant(j),ant(k))=1
                        Bond_k(1,i,ant(j),ant(k))=1
                        Bond_k(1,i,j,ant(k))=1
                        Bond_k(1,ant(i),j,ant(k))=1
                    end if
                end do
            end do
        end do
    end subroutine geraLigacao

    !-----------------------------------------------------------------------------

    !    function isolada(subrede,i,j,k) ! verifica primeiros e segudos vizinhos
    !        !Variaveis mudas
    !        logical :: isolada
    !        byte::  i,j,k
    !        integer ::subrede
    !
    !        if(subrede==1)then
    !            isolada= (ABS(&
    !                sigma(2,i,j,k)* &
    !                sigma(2,suc(i),j,k)* &
    !                sigma(2,suc(i),suc(j),k)* &
    !                sigma(2,i,suc(j),k)* &
    !                sigma(2,i,j,suc(k))* &
    !                sigma(2,suc(i),j,suc(k))* &
    !                sigma(2,suc(i),suc(j),suc(k))* &
    !                sigma(2,i,suc(j),suc(k))&
    !                )==1) .AND. &
    !                (sigma(1,i,j,k)==0).and. &
    !                (sigma(1,ant(i),j,k)* &
    !                sigma(1,suc(i),j,k)* &
    !                sigma(1,i,ant(j),k)* &
    !                sigma(1,i,suc(j),k)* &
    !                sigma(1,i,j,ant(k))* &
    !                sigma(1,i,j,suc(k))==1)
    !        else
    !            isolada= (ABS(&
    !                sigma(1,ant(i),ant(j),k)* &
    !                sigma(1,i,ant(j),k)* &
    !                sigma(1,i,j,k)* &
    !                sigma(1,ant(i),j,k)* &
    !                sigma(1,ant(i),ant(j),ant(k))*&
    !                sigma(1,i,ant(j),ant(k))* &
    !                sigma(1,i,j,ant(k))* &
    !                sigma(1,ant(i),j,ant(k)) &
    !                )==1) .AND. (sigma(2,i,j,k)==0).and. &
    !                (sigma(2,ant(i),j,k)* &
    !                sigma(2,suc(i),j,k)* &
    !                sigma(2,i,ant(j),k)* &
    !                sigma(2,i,suc(j),k)* &
    !                sigma(2,i,j,ant(k))* &
    !                sigma(2,i,j,suc(k))==1)
    !        end if
    !    end function  isolada

    !-----------------------------------------------------------------------------
    function isolada(subrede,i,j,k)! verifica apenas primeiros vizinhos
        !Variaveis mudas
        logical :: isolada
        byte::  i,j,k
        integer ::subrede

        if(subrede==1)then
            isolada= (ABS(&
                sigma(2,i,j,k)* &
                sigma(2,suc(i),j,k)* &
                sigma(2,suc(i),suc(j),k)* &
                sigma(2,i,suc(j),k)* &

                sigma(2,i,j,suc(k))* &
                sigma(2,suc(i),j,suc(k))* &
                sigma(2,suc(i),suc(j),suc(k))* &
                sigma(2,i,suc(j),suc(k))&
                )==1) .AND. &
                (sigma(1,i,j,k)==0)
        else
            isolada= (ABS(&
                sigma(1,ant(i),ant(j),k)* &
                sigma(1,i,ant(j),k)* &
                sigma(1,i,j,k)* &
                sigma(1,ant(i),j,k)* &
                sigma(1,ant(i),ant(j),ant(k))*&
                sigma(1,i,ant(j),ant(k))* &
                sigma(1,i,j,ant(k))* &
                sigma(1,ant(i),j,ant(k)) &
                )==1) .AND. (sigma(2,i,j,k)==0)
        end if
    end function  isolada
    !-----------------------------------------------------------------------------

    subroutine atualiza
        !Variaveis Locais
        double precision:: aux
        byte::  i,j

        do i=-8,8
            do j=-6,6
                aux=(i*J1+j*J2)/t0
                w(i,j)=1
                if (aux >0) w(i,j)=exp(-2*aux)
            end do
        end do
        PJ1= 1- exp(-2*J1/t0)
        PJ2= 1- exp(-2*J2/t0)

    end subroutine atualiza


    !-----------------------------------------------------------------
    subroutine CalcularMagEng
        integer :: i,j,k,NJ1,NJ2
        integer :: soma_nj1
        integer :: soma_nj2
        integer :: somaSigma

        soma_nj1=0
        soma_nj2=0
        somaSigma=0
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(1,i,j,k)*(&
                        sigma(2,i,j,k)+ &
                        sigma(2,suc(i),j,k)+ &
                        sigma(2,suc(i),suc(j),k)+ &
                        sigma(2,i,suc(j),k)+ &
                        sigma(2,i,j,suc(k))+ &
                        sigma(2,suc(i),j,suc(k))+ &
                        sigma(2,suc(i),suc(j),suc(k))+ &
                        sigma(2,i,suc(j),suc(k))   )

                    NJ2= sigma(1,i,j,k)*(&
                        sigma(1,suc(i),j,k)*Bond_i(1,i,j,k)+ &
                        sigma(1,i,suc(j),k)*Bond_j(1,i,j,k)+ &
                        sigma(1,i,j,suc(k))*Bond_k(1,i,j,k) ) +&
                        sigma(2,i,j,k)*( &
                        sigma(2,suc(i),j,k)*Bond_i(2,i,j,k)+ &
                        sigma(2,i,suc(j),k)*Bond_j(2,i,j,k)+ &
                        sigma(2,i,j,suc(k))*Bond_k(2,i,j,k))

                    soma_nj1=NJ1+ soma_nj1
                    soma_nj2=NJ2+ soma_nj2
                end do
            end do
        end do
        energia = -(J1*soma_nj1+ J2*soma_nj2)
        eneJ1=soma_nj1
        eneJ2=soma_nj2
        magnetizacao = abs(sum(sigma(:,:,:,:)))
    end subroutine CalcularMagEng
    !-----------------------------------------------------------------
    subroutine metropolis
        integer i, j, k, NJ1, NJ2
        ! varre subrede=1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(1,i,j,k)*(&
                        sigma(2,i,j,k)+ &
                        sigma(2,suc(i),j,k)+ &
                        sigma(2,suc(i),suc(j),k)+ &
                        sigma(2,i,suc(j),k)+ &
                        sigma(2,i,j,suc(k))+ &
                        sigma(2,suc(i),j,suc(k))+ &
                        sigma(2,suc(i),suc(j),suc(k))+ &
                        sigma(2,i,suc(j),suc(k))   )

                    NJ2= sigma(1,i,j,k)*(&
                        sigma(1,ant(i),j,k)*Bond_i(1,ant(i),j,k)+ &
                        sigma(1,suc(i),j,k)*Bond_i(1,i,j,k)+ &
                        sigma(1,i,ant(j),k)*Bond_j(1,i,ant(j),k)+ &
                        sigma(1,i,suc(j),k)*Bond_j(1,i,j,k)+ &
                        sigma(1,i,j,ant(k))*Bond_k(1,i,j,ant(k))+ &
                        sigma(1,i,j,suc(k))*Bond_k(1,i,j,k) )
                    call random_number(rando)
                    if (rando<W(NJ1,NJ2)) then
                        sigma(1,i,j,k)=-sigma(1,i,j,k)

                    end if
                end do
            end do
        end do

        ! varre subrede=2
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(2,i,j,k)*( &
                        sigma(1,ant(i),ant(j),k)+ &
                        sigma(1,i,ant(j),k)+ &
                        sigma(1,i,j,k)+ &
                        sigma(1,ant(i),j,k)+ &
                        sigma(1,ant(i),ant(j),ant(k))+ &
                        sigma(1,i,ant(j),ant(k))+ &
                        sigma(1,i,j,ant(k))+ &
                        sigma(1,ant(i),j,ant(k)) )

                    NJ2= sigma(2,i,j,k)*(  &
                        sigma(2,ant(i),j,k)*Bond_i(2,ant(i),j,k)+ &
                        sigma(2,suc(i),j,k)*Bond_i(2,i,j,k)+ &
                        sigma(2,i,ant(j),k)*Bond_j(2,i,ant(j),k)+ &
                        sigma(2,i,suc(j),k)*Bond_j(2,i,j,k)+ &
                        sigma(2,i,j,ant(k))*Bond_k(2,i,j,ant(k))+ &
                        sigma(2,i,j,suc(k))*Bond_k(2,i,j,k))
                    call random_number(rando)
                    if(rando<W(NJ1,NJ2)) sigma(2,i,j,k)=-sigma(2,i,j,k)
                end do
            end do
        end do
    end subroutine metropolis
    !-----------------------------------------------------------------
    subroutine initRandomSeed()
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call RANDOM_SEED(size = n)
        allocate(seed(n))
        call SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
    end subroutine initRandomSeed

end program IsingBCC
