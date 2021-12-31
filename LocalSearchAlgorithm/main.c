#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define DEFAULT_RUNS 10
#define PROB 0.01

/*

c FILE: myciel3.col
c SOURCE: Michael Trick (trick@cmu.edu)
c DESCRIPTION: Graph based on Mycielski transformation. 
c              Triangle free (clique number 2) but increasing
c              coloring number
p edge 11 20 <vertices> <arestas = populacao>
e 1 2 // Existe uma aresta entre vertices 1 e 2
e 1 4 // Existe uma aresta entre vertices 1 e 4
e 1 7 // Existe uma aresta entre vertices 1 e 7
e 1 9 // ...
e 2 3

Tamanho da população - arestas
Probabilidades do operador de recombinação / mutação - pm e pr
*/

// algoritmo
int trepa_colinas(int sol[], int *mat, int vert, int num_iter);
int esferiamento(int sol[], int *mat, int vert, int num_iter);

// utils
int* read_file(char *nome, int *n, int *iter);
void gera_sol_inicial(int *sol, int v);
void escreve_sol(int *sol, int vert);
void substitui(int a[], int b[], int n);
void init_rand(void);
int random_l_h(int min, int max);
float rand_01(void);

//funcao
int calcula_fit(int a[], int *mat, int vert);

/*
Ex:
nome do ficheiro e runs (!= iterações)
gcc *.c -o ./dist/main && ./dist/main grafo_20.txt 100
*/

int main(int argc, char *argv[]){
	printf("trepa colinas probabilistico\n");
    char    nome_fich[100] = "",tmp[100] = "";
    int     vert, num_iter, k, runs, custo, best_custo=0;
    int     *grafo, *sol, *best;
	float   mbf = 0.0;

	int prints;
	printf("Print cada solucao? (0/1)\n");
	scanf("%d",&prints);
	int c;
	while ((c = getchar()) != '\n' && c != EOF) { }

	if(argc == 3){
		runs = atoi(argv[2]);
		strcpy(nome_fich, argv[1]);
	}else if(argc == 2){
            runs = DEFAULT_RUNS;
            strcpy(nome_fich, argv[1]);
    }else{
        runs = DEFAULT_RUNS;
        printf("Nome do Ficheiro: ");
        fgets(tmp,100,stdin);
        strtok(tmp, "\n");
    }
    strcat(nome_fich,"./Instancias/");
    strcat(nome_fich,tmp);
    //printf("%s",nome_fich);
	if(runs <= 0)
		return 0;
	init_rand();
    // Preenche matriz de adjacencias
    grafo = read_file(nome_fich, &vert, &num_iter);
	sol = malloc(sizeof(int)*vert);
	best = malloc(sizeof(int)*vert);
	if(sol == NULL || best == NULL){
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


/*  /$$   /$$ /$$$$$$$$ /$$$$$$ /$$        /$$$$$$ 
   | $$  | $$|__  $$__/|_  $$_/| $$       /$$__  $$
   | $$  | $$   | $$     | $$  | $$      | $$  \__/
   | $$  | $$   | $$     | $$  | $$      |  $$$$$$ 
   | $$  | $$   | $$     | $$  | $$       \____  $$
   | $$  | $$   | $$     | $$  | $$       /$$  \ $$
   |  $$$$$$/   | $$    /$$$$$$| $$$$$$$$|  $$$$$$/
    \______/    |__/   |______/|________/ \______/  */


// Leitura do ficheiro de input
// Recebe: nome do ficheiro, numero de vertices (ptr), numero de iteracoes (ptr)
// Devolve a matriz de adjacencias
int* read_file(char *nome, int *n, int *iter){
	FILE *f;
	int *p, *q;
	int i, j;
    int c;

	// Numero de iteracoes
	printf("Escolha o numero de iteracoes: ");
    scanf("%d",iter);
	while ((c = getchar()) != '\n' && c != EOF) { }

	f=fopen(nome, "rt");
	if(!f) { printf("Erro no acesso ao ficheiro dos dados\n"); exit(1); }

    // salta as linhas de comentario
	char ded;
	fscanf(f," %c", &ded);
	while (ded == 'c') {
		fscanf(f," %*[^\n]");
		fscanf(f," %c", &ded);
	}

    if (ded == 'p'){  //  <vertices> <arestas>
		fscanf(f," edge %d %*d",n );
	}

	// Alocacao dinamica da matriz
	p = malloc(sizeof(int)*(*n)*(*n));
	if(!p){ printf("Erro na alocacao de memoria\n");exit(1);}

    memset(p, 0, sizeof(int)*(*n)*(*n));
	q=p;
    // Preenchimento da matriz unidimensional, tratada como bidimensional, cada linha representa um vertice e tem uma coluna para cada outro verticer, onde 0 é que há uma aresta entre eles, e 1 quer dizer que não há
    int vert1, vert2;
    for (i=0; i< *n; i++){
		fscanf(f, " e %d %d", &vert1, &vert2);
        q[(vert1-1)*(*n)+(vert2-1)]=1;
	}

	fclose(f);

    // print da matriz
    int k=0;
    for (i=0; i<*n; i++) {
        for(j=0; j<*n; j++){
            printf("%d ",q[k++]);
        }
        putchar('\n');
    }

	return p;
}

// Gera a solucao inicial
// Parametros: solucao, numero de vertices
void gera_sol_inicial(int *sol, int v)
{
	int i, x;

	for(i=0; i<v; i++)
        sol[i]=0;
	for(i=0; i<v/2; i++)
    {
        do
			x = random_l_h(0, v-1);
        while(sol[x] != 0);
        sol[x]=1;
    }
}

// Escreve solucao
// Parametros: solucao e numero de vertices
void escreve_sol(int *sol, int vert)
{
	int i;

	printf("\nConjunto A: ");
	for(i=0; i<vert; i++)
		if(sol[i]==0)
			printf("%2d  ", i);
	printf("\nConjunto B: ");
	for(i=0; i<vert; i++)
		if(sol[i]==1)
			printf("%2d  ", i);
	printf("\n");
}

// copia vector b para a (tamanho n)
void substitui(int a[], int b[], int n)
{
    int i;
    for(i=0; i<n; i++)
        a[i]=b[i];
}

// Inicializa o gerador de numeros aleatorios
void init_rand()
{
	srand((unsigned)time(NULL));
}

// Devolve valor inteiro aleatorio entre min e max
int random_l_h(int min, int max)
{
	return min + rand() % (max-min+1);
}

// Devolve um valor real aleatorio do intervalo [0, 1]
float rand_01()
{
	return ((float)rand())/(float)RAND_MAX;
}


/*   /$$$$$$  /$$        /$$$$$$   /$$$$$$  /$$$$$$$  /$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$ 
    /$$__  $$| $$       /$$__  $$ /$$__  $$| $$__  $$|_  $$_/|__  $$__/| $$$    /$$$ /$$__  $$
   | $$  \ $$| $$      | $$  \__/| $$  \ $$| $$  \ $$  | $$     | $$   | $$$$  /$$$$| $$  \ $$
   | $$$$$$$$| $$      | $$ /$$$$| $$  | $$| $$$$$$$/  | $$     | $$   | $$ $$/$$ $$| $$  | $$
   | $$__  $$| $$      | $$|_  $$| $$  | $$| $$__  $$  | $$     | $$   | $$  $$$| $$| $$  | $$
   | $$  | $$| $$      | $$  \ $$| $$  | $$| $$  \ $$  | $$     | $$   | $$\  $ | $$| $$  | $$
   | $$  | $$| $$$$$$$$|  $$$$$$/|  $$$$$$/| $$  | $$ /$$$$$$   | $$   | $$ \/  | $$|  $$$$$$/
   |__/  |__/|________/ \______/  \______/ |__/  |__/|______/   |__/   |__/     |__/ \______/  */


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

/*  /$$$$$$$$ /$$   /$$ /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$ 
   | $$_____/| $$  | $$| $$$ | $$ /$$__  $$ /$$__  $$ /$$__  $$
   | $$      | $$  | $$| $$$$| $$| $$  \__/| $$  \ $$| $$  \ $$
   | $$$$$   | $$  | $$| $$ $$ $$| $$      | $$$$$$$$| $$  | $$
   | $$__/   | $$  | $$| $$  $$$$| $$      | $$__  $$| $$  | $$
   | $$      | $$  | $$| $$\  $$$| $$    $$| $$  | $$| $$  | $$
   | $$      |  $$$$$$/| $$ \  $$|  $$$$$$/| $$  | $$|  $$$$$$/
   |__/       \______/ |__/  \__/ \______/ |__/  |__/ \______/  */

int calcula_fit(int a[], int *mat, int vert)
{
	int total=0;
	int i, j;

	for(i=0; i<vert; i++)
		if(a[i]==0)
		{
			for(j=0; j<vert;j++)
				if(a[j]==1 && *(mat+i*vert+j)==1)
				    total++;
		}
	return total;
}  