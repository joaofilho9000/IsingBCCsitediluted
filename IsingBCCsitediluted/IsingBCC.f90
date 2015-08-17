program IsingBCC
    implicit none
    save

    !decalração de variaveis
    intrinsic random_seed, random_number
    byte, parameter ::L=30
    integer ::MCx, MCc
    integer ::passo
    integer ::magnetizacao
    double precision energia
    double precision  somaMag, somaMag2, somaMag4
    double precision  somaEne, somaEne2, somaEne4

    double precision  mediaMag,  mediaMag2, mediaMag4
    double precision  mediaEne,mediaEne2, mediaEne4
    double precision  susceptibilidade, cumuM
    double precision  calor, cumuE


    integer ::sigma(-1:1,1:L,1:L,1:L)
    byte ::ant(1:L), suc(1:L)
    byte ::bond_i(-1:1,1:L,1:L,1:L)
    byte ::bond_j(-1:1,1:L,1:L,1:L)
    byte ::bond_k(-1:1,1:L,1:L,1:L)
    real::t0, p   !p=0 sistema puro
    real::J2
    real::tin,tfi,dt
    double precision ::w(-8:8,-6:6)
    real    ::rando
    byte :: histograma
    !Abrindo arquivos
    open(22,file='hist.dat')
    open(12,file='dados.dat')

    !inicilização das variaveis

    call leDados

    write(*,*) 'iniciando simulação com :'
    write(*,*) 'p =', p
    write(*,*) 't0 =', t0
    write(*,*) 'J2 =', J2
    write(*,*) 'MCx =', MCx
    write(*,*) 'MCc =', MCc
    call initRandomSeed()
    call iniciaContorno
    call iniciaSigma
    call iniciaBond
    call dilui
    call geraLigacao
    do t0 = tin , tfi, dt
        call iniciaVariaveis
        call atualiza(t0,W)
        !repetições termalização
        do passo = 1 , MCx
            call metropolis
        end do
        if (histograma==1) then
            write(*,*) 'ok'
            do passo = 1 , MCc
                call metropolis
                call salvar
            end do
        else
            do passo = 1 , MCc
                call metropolis
                call calcularMagEng
                call calcularSoma
            end do
            call calcularMedia
        end if
    end do
    !fechando arquivos
    close(10)
    close(20)
    write(*,*) 'fim do programa'
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
        cumuM=1-mediaMag4/(3*mediaMag2*mediaMag2)  ! rever essa equação soma ou media
        write(*, *)t0, susceptibilidade, calor, cumuM, mediaMag, mediaEne
    end subroutine calcularMedia
     !-----------------------------------------------------------------------------
    subroutine leDados
        read(12,*) p
        read(12,*) tin
        read(12,*) tfi
        read(12,*) dt
        read(12,*) j2
        read(12,*) MCx
        read(12,*) MCc
        read(12,*) histograma
    end subroutine leDados

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
        byte::  i,j,k
        integer ::sinal
        !        do sinal =-1,1,2
        !            do i = 1 , L
        !                do j = 1 , L
        !                    do k = 1 , L
        !                        sigma(sinal,i,j,k)=1
        !                    end do
        !                end do
        !            end do
        !        end do
        sigma=1
    end subroutine iniciaSigma
    !-----------------------------------------------------------------------------
    subroutine iniciaBond
        byte::  i,j,k
        integer ::sinal
        !        do sinal =-1,1,2
        !            do i = 1 , L
        !                do j = 1 , L
        !                    do k = 1 , L
        !                        bond_i(sinal,i,j,k)=0
        !                        bond_j(sinal,i,j,k)=0
        !                        bond_k(sinal,i,j,k)=0
        !                    end do
        !                end do
        !            end do
        !        end do
        bond_i=0
        bond_j=0
        bond_k=0
    end subroutine iniciaBond
    !-----------------------------------------------------------------------------
    subroutine  dilui
        integer :: sinal,n
        byte::  i,j,k

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
            aux=rando+0.5
            sinal=int(aux/abs(aux))
            if (abs(sigma(sinal, i,j,k))==1) then
                sigma(sinal, i,j,k)=0
                n=n+1
            end if
        end do
    end subroutine dilui

    !-----------------------------------------------------------------------------
    subroutine geraLigacao
        byte::  i,j,k
        integer ::sinal
        integer :: n
        n=0
        sinal=1
        do i=1,L
            do j=1,L
                do k=1,L
                    if( isolada(sinal,i,j,k) )then
                        bond_i(-sinal,i,j,k)=1
                        bond_i(-sinal,i,suc(j),k)=1
                        bond_i(-sinal,i,j,suc(k))=1
                        bond_i(-sinal,i,suc(j),suc(k))=1

                        bond_j(-sinal,i,j,k)=1
                        bond_j(-sinal,suc(i),j,k)=1
                        bond_j(-sinal,i,j,suc(k))=1
                        bond_j(-sinal,suc(i),j,suc(k))=1

                        bond_k(-sinal,i,j,k)= 1
                        bond_k(-sinal,suc(i),j,k)= 1
                        bond_k(-sinal,i,suc(j),k)= 1
                        bond_k(-sinal,suc(i),suc(j),k)= 1
                    end if
                end do
            end do
        end do
        sinal=-1
        do i=1,L
            do j=1,L
                do k=1,L
                    if( isolada(sinal,i,j,k) )then
                        Bond_i(-sinal,ant(i),ant(j),ant(k))=1
                        Bond_i(-sinal,ant(i),    j, ant(k))=1
                        Bond_i(-sinal,ant(i),ant(j),    k)=1
                        Bond_i(-sinal,ant(i),    j,     k)=1

                        Bond_j(-sinal,ant(i), ant(j),ant(k))=1
                        Bond_j(-sinal,ant(i), ant(j),    k)=1
                        Bond_j(-sinal,    i,  ant(j),ant(k))=1
                        Bond_j(-sinal,    i,  ant(j),    k)=1

                        Bond_k(-sinal,ant(i),ant(j),ant(k))=1
                        Bond_k(-sinal,i,ant(j),ant(k))=1
                        Bond_k(-sinal,i,j,ant(k))=1
                        Bond_k(-sinal,ant(i),j,ant(k))=1
                    end if
                end do
            end do
        end do
    end subroutine geraLigacao

    !-----------------------------------------------------------------------------

    function isolada(sinal,i,j,k)
        !Variaveis mudas
        logical :: isolada
        byte::  i,j,k
        integer ::sinal

        if(sinal==1)then
            isolada= (ABS(&
                sigma(-sinal,i,j,k)* &
                sigma(-sinal,suc(i),j,k)* &
                sigma(-sinal,suc(i),suc(j),k)* &
                sigma(-sinal,i,suc(j),k)* &

                sigma(-sinal,i,j,suc(k))* &
                sigma(-sinal,suc(i),j,suc(k))* &
                sigma(-sinal,suc(i),suc(j),suc(k))* &
                sigma(-sinal,i,suc(j),suc(k))&
                )==1) .AND. &
                (sigma(sinal,i,j,k)==0)
        end if
        if(sinal==-1)then
            isolada= (ABS(&
                sigma(-sinal,ant(i),ant(j),k)* &
                sigma(-sinal,i,ant(j),k)* &
                sigma(-sinal,i,j,k)* &
                sigma(-sinal,ant(i),j,k)* &
                sigma(-sinal,ant(i),ant(j),ant(k))*&
                sigma(-sinal,i,ant(j),ant(k))* &
                sigma(-sinal,i,j,ant(k))* &
                sigma(-sinal,ant(i),j,ant(k)) &
                )==1) .AND. (sigma(sinal,i,j,k)==0)
        end if
    end function  isolada

    !-----------------------------------------------------------------------------

    subroutine atualiza(t0,w)
        !Variaveis mudas
        real :: t0
        double precision :: w(-8:8,-6:6)
        !Variaveis Locais
        double precision:: aux
        byte::  i,j

        do i=-8,8
            do j=-6,6
                aux=(i+j*J2)/t0
                w(i,j)=1
                if (aux >0) w(i,j)=exp(-2*aux)
            end do
        end do
    end subroutine atualiza

    !-----------------------------------------------------------------------------
    subroutine salvar
        integer :: i,j,k,sinal, NJ1,NJ2
        integer :: soma_nj1
        integer :: soma_nj2 , magnetizacao
        integer :: somaSigma

        soma_nj1=0
        soma_nj2=0
        somaSigma=0

        sinal=1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(sinal,i,j,k)*(&
                        sigma(-sinal,i,j,k)+ &
                        sigma(-sinal,suc(i),j,k)+ &
                        sigma(-sinal,suc(i),suc(j),k)+ &
                        sigma(-sinal,i,suc(j),k)+ &
                        sigma(-sinal,i,j,suc(k))+ &
                        sigma(-sinal,suc(i),j,suc(k))+ &
                        sigma(-sinal,suc(i),suc(j),suc(k))+ &
                        sigma(-sinal,i,suc(j),suc(k))   )

                    NJ2= sigma(sinal,i,j,k)*(&
                        sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                        sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                        sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                        sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                        sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                        sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k) )
                    soma_nj1=NJ1+ soma_nj1
                    soma_nj2=NJ2+ soma_nj2
                end do
            end do
        end do

        sinal=-1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(sinal,i,j,k)*( &
                        sigma(-sinal,ant(i),ant(j),k)+ &
                        sigma(-sinal,i,ant(j),k)+ &
                        sigma(-sinal,i,j,k)+ &
                        sigma(-sinal,ant(i),j,k)+ &
                        sigma(-sinal,ant(i),ant(j),ant(k))+ &
                        sigma(-sinal,i,ant(j),ant(k))+ &
                        sigma(-sinal,i,j,ant(k))+ &
                        sigma(-sinal,ant(i),j,ant(k)) )

                    NJ2= sigma(sinal,i,j,k)*(  &
                        sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                        sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                        sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                        sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                        sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                        sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k))

                    soma_nj1=NJ1+ soma_nj1
                    soma_nj2=NJ2+ soma_nj2

                end do
            end do
        end do

        !        do sinal=-1,1,2
        !            do i=1,L
        !                do j=1,L
        !                    do k=1,L
        !                        somaSigma= somaSigma + sigma(sinal,i,j,k)
        !                    end do
        !                end do
        !            end do
        !        end do
        magnetizacao = sum(sigma(:,:,:,:))
        write(22,*) soma_nj1, soma_nj2,  magnetizacao
    end subroutine salvar
    !-----------------------------------------------------------------
    subroutine CalcularMagEng
        integer :: i,j,k,sinal, NJ1,NJ2
        integer :: soma_nj1
        integer :: soma_nj2
        integer :: somaSigma

        soma_nj1=0
        soma_nj2=0
        somaSigma=0

        sinal=1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(sinal,i,j,k)*(&
                        sigma(-sinal,i,j,k)+ &
                        sigma(-sinal,suc(i),j,k)+ &
                        sigma(-sinal,suc(i),suc(j),k)+ &
                        sigma(-sinal,i,suc(j),k)+ &
                        sigma(-sinal,i,j,suc(k))+ &
                        sigma(-sinal,suc(i),j,suc(k))+ &
                        sigma(-sinal,suc(i),suc(j),suc(k))+ &
                        sigma(-sinal,i,suc(j),suc(k))   )

                    NJ2= sigma(sinal,i,j,k)*(&
                         sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                         sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                         sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k) )
                    soma_nj1=NJ1+ soma_nj1
                    soma_nj2=NJ2+ soma_nj2
                end do
            end do
        end do

        sinal=-1
        do i=1,L
            do j=1,L
                do k=1,L

                    NJ2= sigma(sinal,i,j,k)*(  &
                        sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                        sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                        sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k))


                    soma_nj2=NJ2+ soma_nj2

                end do
            end do
        end do

                do sinal=-1,1,2
                    do i=1,L
                        do j=1,L
                            do k=1,L
                                somaSigma= somaSigma + sigma(sinal,i,j,k)
                            end do
                        end do
                    end do
                end do


        energia = -(soma_nj1+ J2*soma_nj2)
        magnetizacao = somaSigma!sum(sigma(:,:,:,:))
    end subroutine CalcularMagEng
    !-----------------------------------------------------------------
    subroutine metropolis
        integer i, j, k, sinal, NJ1, NJ2
        sinal=1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(sinal,i,j,k)*(&
                        sigma(-sinal,i,j,k)+ &
                        sigma(-sinal,suc(i),j,k)+ &
                        sigma(-sinal,suc(i),suc(j),k)+ &
                        sigma(-sinal,i,suc(j),k)+ &
                        sigma(-sinal,i,j,suc(k))+ &
                        sigma(-sinal,suc(i),j,suc(k))+ &
                        sigma(-sinal,suc(i),suc(j),suc(k))+ &
                        sigma(-sinal,i,suc(j),suc(k))   )

                    NJ2= sigma(sinal,i,j,k)*(&
                        sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                        sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                        sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                        sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                        sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                        sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k) )
                    call random_number(rando)
                    if (rando<W(NJ1,NJ2)) then
                        sigma(sinal,i,j,k)=-sigma(sinal,i,j,k)
                    end if
                end do
            end do
        end do

        sinal=-1
        do i=1,L
            do j=1,L
                do k=1,L
                    NJ1= sigma(sinal,i,j,k)*( &
                        sigma(-sinal,ant(i),ant(j),k)+ &
                        sigma(-sinal,i,ant(j),k)+ &
                        sigma(-sinal,i,j,k)+ &
                        sigma(-sinal,ant(i),j,k)+ &
                        sigma(-sinal,ant(i),ant(j),ant(k))+ &
                        sigma(-sinal,i,ant(j),ant(k))+ &
                        sigma(-sinal,i,j,ant(k))+ &
                        sigma(-sinal,ant(i),j,ant(k)) )

                    NJ2= sigma(sinal,i,j,k)*(  &
                        sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                        sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                        sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                        sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                        sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                        sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k))
                    call random_number(rando)
                    if(rando<W(NJ1,NJ2)) sigma(sinal,i,j,k)=-sigma(sinal,i,j,k)
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
