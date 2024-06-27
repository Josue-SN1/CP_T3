#include <stdio.h>
#include "mpi.h"
#include <time.h>
#include <stdlib.h>

#define ROOT 0
#define TRUE 1
#define FALSE 0
#define LEFT_CHILD(my_rank) 2 * my_rank + 1
#define RIGHT_CHILD(my_rank) 2 * my_rank + 2
#define FATHER(my_rank) (my_rank - 1)/2

void BubbleSort(int *vector, int size, int rank)
{

    for (int i = 0; i < size - 1; i++)
    {
        for (int j = 0; j < size - 1 - i; j++)
        {
            if (vector[j] > vector[j + 1])
            {
                int temp = vector[j];
                vector[j] = vector[j + 1];
                vector[j + 1] = temp;
            }
        }
    }
}

int *initializeVector(int n)
{
    // Alocar memória para o vetor no heap
    int *vetor = (int *)malloc(n * sizeof(int));
    if (vetor == NULL)
    {
        printf("Erro de alocação de memória.\n");
        return NULL;
    }

    // Preencher o vetor com valores ao contrário
    for (int i = 0; i < n; i++)
    {
        vetor[i] = n - i;
    }
    return vetor;
}

void interleaveArray(int **arr, int size)
{
    int mid = size / 2;
    int *result = (int *)malloc(size * sizeof(int));

    if (result == NULL)
    {
        fprintf(stderr, "Erro ao alocar memória\n");
        return;
    }

    int i = 0, j = mid, k = 0;

    while (i < mid && j < size)
    {
        if ((*arr)[i] <= (*arr)[j])
        {
            result[k++] = (*arr)[i++];
        }
        else
        {
            result[k++] = (*arr)[j++];
        }
    }

    while (i < mid)
    {
        result[k++] = (*arr)[i++];
    }

    while (j < size)
    {
        result[k++] = (*arr)[j++];
    }

    // Libera o espaço do array original
    free(*arr);

    // Faz o ponteiro arr apontar para o novo array result
    *arr = result;
}

void SHOW(int *vetor, int SIZE)
{
    for (size_t i = 0; i < SIZE; i++)
    {
        printf("%d  ", vetor[i]);
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int my_rank;
    MPI_Status status;
    int world_size;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int vector_size = 1000000;
    int DELTA = vector_size/((world_size+1)/2);
    int *vector = (int *)malloc(vector_size * sizeof(int));
    
    double start;
    double end;

    switch (my_rank)
    {
    case ROOT:
        start = MPI_Wtime();
	vector = initializeVector(vector_size);
        printf("Rodando com vetor de tamanaho %d,  DELTA  %d e quantidade de processos %d\n\n", vector_size, DELTA, world_size);
	printf("Tempo inicial: %f\n\n",(double)start);
        break;

    default:
        MPI_Recv(vector, vector_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // não sou a raiz, tenho pai
        MPI_Get_count(&status, MPI_INT, &vector_size);                                                 // descubro tamanho da mensagem recebida
        break;
    }

    // dividir ou conquistar?

    switch (vector_size <= DELTA)
    {
    case TRUE:
        BubbleSort(vector, vector_size, my_rank); // conquisto
        break;

    case FALSE:
        MPI_Send(&vector[0], vector_size / 2, MPI_INT, LEFT_CHILD(my_rank), 1, MPI_COMM_WORLD);
        MPI_Send(&vector[vector_size / 2], (vector_size / 2) + (vector_size % 2), MPI_INT, RIGHT_CHILD(my_rank), 1, MPI_COMM_WORLD);

        MPI_Recv(&vector[0], vector_size, MPI_INT, LEFT_CHILD(my_rank), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&vector[vector_size / 2], vector_size, MPI_INT, RIGHT_CHILD(my_rank), MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        interleaveArray(&vector, vector_size);
        break;
    }

    switch (my_rank)
    {
    case ROOT:
	end = MPI_Wtime();
	printf("Tempo final: %f\n\n",end);
	double finalTime = end - start;
//	SHOW(vector,vector_size);
	printf("Final time: %f\n\n",finalTime);
	printf("Time taken to sort the vectors: %f seconds\n", finalTime);
	
        break;
    default:
        MPI_Send(&vector[0], vector_size, MPI_INT, FATHER(my_rank), 1, MPI_COMM_WORLD);
        break;
    }

    MPI_Finalize();

    free(vector);
}

