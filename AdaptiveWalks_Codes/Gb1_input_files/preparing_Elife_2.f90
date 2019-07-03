! --- Programa que le dados e calcula a distancia de hamming --- !

program lendo_dados
implicit none 

! --- Declarando variaveis importantes 
integer     M, nViz, Viz
dimension   M(0 : 159999, 4), nViz(0 : 159999), Viz(0 : 159999, 76)

! --- Declarando variaveis adicionais 
integer     i, j, k, soma

! --- Abrindo arquivo 
open(10, file = "elife_numero_sequencia_numerica.txt")

! --- Percorrendo arquivo e lendo os dados 
do i = 0, 159999
    ! --- Lendo dados 
    read(10, *) j, M(i, 1), M(i, 2), M(i, 3), M(i, 4)
end do 

! --- Fechando arquivo 
close(10)

! --- Zerando o vetor que armazena o numero de vizinhos 
nViz = 0

! --- Percorrendo os pares de sequencia 
do i = 0, 159998
    do j = (i + 1), 159999
        ! --- Zerando a variavel que conta os elementos iguais
        soma = 0
        ! --- Percorrendo as sequencias e comparando 
        do k = 1, 4 
            if (M(i, k) .eq. M(j, k)) soma = soma + 1
        end do 
     
        ! --- Verificando se devemos afirmar as sequencias i, j como vizinhas 
        if (soma .eq. 3) then
            ! --- Adicionando o vizinho na sequencia i 
            nViz(i) = nViz(i) + 1 
            Viz(i, nViz(i)) = j 

            ! --- Adicionando o vizinho na sequencia j 
            nViz(j) = nViz(j) +  1 
            Viz(j, nViz(j)) = i
        end if

    end do
end do 

! --- Imprimindo resultado 
open(10, file = "elife_sequence_degree.txt")
open(11, file = "elife_sequence_neighbors_correta.txt")
do i = 0, 159999
    write(10, *) i, nViz(i)
    do j = 1, nViz(i)
        write(11, *) i, (j - 1), Viz(i, j)
    end do 
end do 
close(10)
close(11)

! --- Encerrando o programa 
stop 
end