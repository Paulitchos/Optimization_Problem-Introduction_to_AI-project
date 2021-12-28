#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "algoritmo.h"
#include "utils.h"

#define DEFAULT_RUNS 10

/*
Ex:
nome do ficheiro e runs (!= iterações)
gcc *.c -o ./dist/main && ./dist/main grafo_20.txt 100
*/

int main(int argc, char *argv[])
{
	printf("trepa colinas probabilistico\n");
    char    nome_fich[100];
    int     vert, num_iter, k, runs, custo, best_custo=0;
    int     *grafo, *sol, *best;
	float   mbf = 0.0;

	int prints;
	printf("Print cada sol? (0/1)\n");
	scanf("%d",&prints);
	int c;
	while ((c = getchar()) != '\n' && c != EOF) { }

	if(argc == 3)
	{
		runs = atoi(argv[2]);
		strcpy(nome_fich, argv[1]);
	}
	else
        if(argc == 2)
        {
            runs = DEFAULT_RUNS;
            strcpy(nome_fich, argv[1]);
        }
        else
        {
            runs = DEFAULT_RUNS;
            printf("Nome do Ficheiro: ");
            gets(nome_fich);
        }
	if(runs <= 0)
		return 0;
	init_rand();
    // Preenche matriz de adjacencias
    grafo = init_dados(nome_fich, &vert, &num_iter);
	sol = malloc(sizeof(int)*vert);
	best = malloc(sizeof(int)*vert);
	if(sol == NULL || best == NULL)
	{
		printf("Erro na alocacao de memoria");
		exit(1);
	}
	for(k=0; k<runs; k++){
		// Gerar solucao inicial
		gera_sol_inicial(sol, vert);
		if (prints) escreve_sol(sol, vert); // ***
		// Trepa colinas
		//custo = trepa_colinas(sol, grafo, vert, num_iter);
		custo = esferiamento(sol, grafo, vert, num_iter); //esfriamento
		// Escreve resultados da repeticao k
		if (prints)  printf("\nRepeticao %d:", k); // ***
		if (prints) escreve_sol(sol, vert); // ***
		if (prints) printf("Custo final: %2d\n", custo); // ***
		mbf += custo;
		if(k==0 || best_custo > custo)
		{
			best_custo = custo;
			substitui(best, sol, vert);
		}
    }
	// Escreve eresultados globais
	printf("\n====================\n");
	printf("\nMBF(mean best feat): %f\n", mbf/k); //médio
	printf("\nMelhor solucao encontrada\n");
	if (prints) escreve_sol(best, vert); // ***
	printf("Custo final (melhor): %2d\n", best_custo);
	printf("\n====================\n");
	free(grafo);
    free(sol);
	free(best);
    return 0;
}
