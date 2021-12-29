#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "algoritmo.h"
#include "funcao.h"
#include "utils.h"

#define PROB 0.01

// old one: "c": "cd $dir && gcc $fileName -o /tmp/$fileNameWithoutExt -lm && /tmp/$fileNameWithoutExt",
// windows: "c": "cd $dir && gcc *.c -o C:/Users/yeshe/AppData/Local/Temp/$fileNameWithoutExt -lm && C:/Users/yeshe/AppData/Local/Temp/$fileNameWithoutExt",

// Gera um vizinho
// Parametros: solucao actual, vizinho, numero de vertices
//swap two vertices
void gera_vizinho(int a[], int b[], int n)
{
    int i, p1, p2;

    for(i=0; i<n; i++)
        b[i]=a[i];
	// Encontra posicao com valor 0
    do
        p1=random_l_h(0, n-1);
    while(b[p1] != 0);
	// Encontra posicao com valor 0
    do
        p2=random_l_h(0, n-1);
    while(b[p2] != 1);
	// Troca
    b[p1] = 1;
    b[p2] = 0;
}

/*if erro > 0 then current <- $next 
    else current <- next with probibility e^(erro/t)

    */

int esferiamento(int sol[], int *mat, int vert, int num_iter) { // trepa_colinas alterado

    

    int *nova_sol, custo, custo_viz, i; 
    float temp, erro, prob_aceitar, max = 100, min = 5; // Temperaturas inicial e final

	nova_sol = malloc(sizeof(int)*vert);
    if(nova_sol == NULL){ printf("Erro na alocacao de memoria"); exit(1); }
	// Avalia solucao inicial
    custo = calcula_fit(sol, mat, vert); //fico a saber qual é o custo da minha solução atual
    temp = max;
    
    while (temp > min) {
        temp -= (max-min) / (float)num_iter; // Descer temperatura
        // Gera vizinho
		gera_vizinho(sol, nova_sol, vert);
		// Avalia vizinho
		custo_viz = calcula_fit(nova_sol, mat, vert);
        if(custo_viz<=custo) {
            substitui(sol, nova_sol, vert);
            custo = custo_viz;
        } else {
            erro = custo_viz - custo; //tem que ser modulo
            // Aceitar com determinada probabiçlidade - calcular probabilidade usando funcao exp()
            /*  if erro > 0 then current <- $next 
                else current <- next with probibility e^(erro/temp) */

            erro = custo - custo_viz; prob_aceitar = exp(erro/temp);
            if (prob_aceitar>rand_01()){
                substitui(sol, nova_sol, vert); custo = custo_viz;
            }
            /*    
            if (erro > 0) {
                custo = custo_viz;
            } else if ( rand_01() < exp(erro/temp) ) {
                custo = custo_viz;
            }
            */
        }
    }
    free(nova_sol);
    return custo;
}

// Trepa colinas first-choice
// Parametros: solucao, matriz de adjacencias, numero de vertices e numero de iteracoes
// Devolve o custo da melhor solucao encontrada
int trepa_colinas(int sol[], int *mat, int vert, int num_iter)
{
    int *nova_sol, custo, custo_viz, i; 
    int *nova_sol2, custo_viz2;  // Vizinhança 2

	nova_sol = malloc(sizeof(int)*vert);
    nova_sol2 = malloc(sizeof(int)*vert); // Vizinhança 2
    if(nova_sol == NULL)
    {
        printf("Erro na alocacao de memoria");
        exit(1);
    }
	// Avalia solucao inicial
    custo = calcula_fit(sol, mat, vert); //fico a saber qual é o custo da minha solução atual
    // continuar, aceitar a melhor entre sol, nova sol e nova sol 2
    for(i=0; i<num_iter; i++)
    {
		// Gera vizinho
		gera_vizinho(sol, nova_sol, vert);
        gera_vizinho(sol, nova_sol2, vert);
		// Avalia vizinho
		custo_viz = calcula_fit(nova_sol, mat, vert);
        //custo_viz2 = calcula_fit(nova_sol2, mat, vert); // Vizinhança 2
		// Aceita vizinho se o custo diminuir (problema de minimizacao)
        // se o custo for menor, ent, aceito, substituo a solução atual pela nova

        if (custo_viz <= custo) {
            substitui(sol, nova_sol, vert);
            custo=custo_viz;
        } else { //solução pior tmb pode ser aceite
            if (rand_01() < PROB) { //isto ajuda a fugir de máximos locais
                substitui(sol, nova_sol, vert);
                custo=custo_viz;
            }
        }

        /* // Vizinhança 2
        if(custo_viz <= custo_viz2) // mudado de < para <=
        {
            if (custo_viz < custo){
                substitui(sol, nova_sol, vert);
                custo=custo_viz;
            }
        } else { // viz 2 é melhor
            if (custo_viz2 < custo){
                substitui(sol, nova_sol2, vert);
                custo=custo_viz2;
            }
        }
        */
    }
    free(nova_sol);
    free(nova_sol2);
    return custo;
}
