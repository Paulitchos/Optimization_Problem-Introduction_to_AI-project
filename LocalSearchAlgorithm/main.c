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
void escreve_sol2(int *sol, int vert);
// algoritmo
int trepa_colinas(int sol[], int *mat, int vert, int num_iter);
int esferiamento(int sol[], int *mat, int vert, int num_iter, int print);
int flip();

// utils
int* read_file(char *nome, int *n, int *iter);
void gera_sol_inicial(int *sol, int v, int * mat, int prints);
void escreve_sol(int *sol, int vert);
void substitui(int a[], int b[], int n);
void init_rand(void);
int random_l_h(int min, int max);
float rand_01(void);
void printsol(int * sol, int vert);

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
    int     *mat, *sol, *best;
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

    mat = read_file(nome_fich, &vert, &num_iter);
    // Preenche matriz de adjacencias
    sol = malloc(sizeof(int)*vert);
    best = malloc(sizeof(int)*vert);
    if(sol == NULL || best == NULL){
        printf("Erro na alocacao de memoria");
        exit(1);
    }
    for(k=0; k<runs; k++){
        // Gerar solucao inicial
        gera_sol_inicial(sol, vert, mat, prints);
        if (prints) escreve_sol(sol, vert); // ***
        // Trepa colinas
        //custo = trepa_colinas(sol, grafo, vert, num_iter);
        custo = trepa_colinas(sol, mat, vert, num_iter);
        //custo = esferiamento(sol, mat, vert, num_iter, prints); //esfriamento // sol fica com a melhor soulução e retorna o custo dela
        // Escreve resultados da repeticao k
        printf("\nRepeticao %d:", k); // ***
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
    printf("\n==========RESULTADOS==========\n");
    printf("\nMelhor solucao encontrada\n");
    if (prints) escreve_sol(best, vert); // ***

    printf("\nMBF(mean best feat): %f\n", mbf/k); //médio
    printf("Quantos verts em media das runs: %f\n", vert-(mbf/k));
    printf("Custo final (melhor): %2d \n", best_custo); // (quanto menor melhor, = número de verts que nao pertencem a sol)
    printf("Número de vertices da melhor sol: %2d ", vert-best_custo);
    printf("\n====================\n");
    free(mat);
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
int* read_file(char *nome, int *verts, int *num_iter){
	FILE *f;
	int *p, *q;
	int i, j;
    int c;
    int arestas;

	// Numero de iteracoes
	printf("Escolha o numero de iteracoes: ");
    scanf("%d",num_iter);
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
		fscanf(f," edge %d %d",verts, &arestas );
	}

	// Alocacao dinamica da matriz
	p = malloc(sizeof(int)*(*verts)*(*verts));
	if(!p){ printf("Erro na alocacao de memoria\n");exit(1);}

    memset(p, 0, sizeof(int)*(*verts)*(*verts));
	q=p;
    // Preenchimento da matriz unidimensional, tratada como bidimensional, cada linha representa um vertice e tem uma coluna para cada outro verticer, onde 0 é que há uma aresta entre eles, e 1 quer dizer que não há
    int vert1, vert2;
    for (i=0; i< arestas; i++){
		fscanf(f, " e %d %d", &vert1, &vert2);
        q[(vert1-1)*(*verts)+(vert2-1)]=1;
        q[(vert2-1)*(*verts)+(vert1-1)]=1;
	}

	fclose(f);

    // print da matriz
    int k=0;
    for (i=0; i<*verts; i++) {
        for(j=0; j<*verts; j++){
            printf("%d ",q[k++]);
        }
        putchar('\n');
    }

	return p;
}

int solucaovalida(int *sol, int v, int * mat,int prints);

// Gera a solucao inicial
// Parametros: solucao, numero de vertices
void gera_sol_inicial(int *sol, int v, int * mat, int prints){

    int i, x;

    do {

        for(i=0; i<v; i++) //meter tudo a zeros no array sol de tamanho v(vértices)
            sol[i]=0;

        int num_verts_sol = random_l_h(0, v-1);
        for(i=0; i<num_verts_sol; i++){
            // mete aleatoriamente um dos indexs que estão a zero a um
            do
                x = random_l_h(0, v-1);
            while(sol[x] != 0);
            sol[x]=1;
        }
    
    //escreve_sol2(sol, v);

    } while (solucaovalida(sol, v, mat, prints) == 0);


    /*
	int i, x;

	for(i=0; i<v; i++) //meter tudo a zeros no array sol de tamanho v(vértices)
        sol[i]=0;
	for(i=0; i<v/2; i++){
        // mete aleatoriamente um dos indexs que estão a zero a um
        do
			x = random_l_h(0, v-1);
        while(sol[x] != 0);
        sol[x]=1;
    }
    */
}

int solucaovalida(int *sol, int v, int * mat, int prints){
    int i=0, j=0;
    for(i=0; i<v; i++)
		if(sol[i]==1){
			for(j=0; j<v;j++)
				if(sol[j]==1 && mat[i*v+j]==1){
				    if (prints) printf(" sol inval ");
                    return 0; 
                }
		}
    return 1;
}

// Escreve solucao
// Parametros: solucao e numero de vertices
void escreve_sol(int *sol, int vert)
{
	int i;

	printf("\nVertices que nao pertecem a sol: ");
	for(i=0; i<vert; i++)
		if(sol[i]==0)
			printf("%2d  ", i+1);
	printf("\nVertices que pertecem a sol: ");
	for(i=0; i<vert; i++)
		if(sol[i]==1)
			printf("%2d  ", i+1);
	printf("\n");
}

void escreve_sol2(int *sol, int vert)
{
    int i;

    for(i=0;i<vert;i++) {
        printf(" %d |", sol[i]);
    }
    printf("\n\n\n");
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


int tudoAUns(int *sol, int v);
int tudoAZeros(int *sol, int v);

// Gera um vizinho
// Parametros: solucao actual, vizinho, numero de vertices
//swap two vertices
void gera_vizinho(int sol[], int nova_sol[], int vert)
{
    int i, p1, p2;

    for(i=0; i<vert; i++)
        nova_sol[i]=sol[i];

	// Encontra posicao com valor 0
    //if (random_l_h(0,1) == 1 ){
    if (flip()==1){ 
        if (tudoAUns(sol, vert) == 0){
            do
                p1=random_l_h(0, vert-1);
            while(nova_sol[p1] != 0);
        }
        nova_sol[p1] = 1;
    }
   
        
    //}

	// Encontra posicao com valor 1
    if (flip()==1){ 
        if (tudoAZeros(sol,vert) == 0){
            do
                p2=random_l_h(0, vert-1);
            while(nova_sol[p2] != 1);
        }
        nova_sol[p2] = 0;
    }
	// Troca
    
    
}

int tudoAUns(int *sol, int v){
    for(int i=0; i<v; i++)
        if (sol[i]!=1) 
            return 0; 
    return 1;
}

int tudoAZeros(int *sol, int v){
    for(int i=0; i<v; i++)
        if (sol[i]!=0) 
            return 0; 
    return 1;
}

void printsol(int * sol, int vert){
    for(int i=0; i<vert; i++)
        printf("%d ", sol[i]);
    putchar('|');
}
/*if erro > 0 then current <- $next 
    else current <- next with probibility e^(erro/t)

    */

   /*
   sol - [0, 0 ,1, 0 , 1, 0 , 1, 0...
   mat - tabela com as ligações de todos os vertices
   */

int esferiamento(int sol[], int *mat, int vert, int num_iter, int prints) { // trepa_colinas alterado

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
        do {
            gera_vizinho(sol, nova_sol, vert);
        } while (solucaovalida(nova_sol, vert, mat, prints) == 0);
		
		// Avalia vizinho
		custo_viz = calcula_fit(nova_sol, mat, vert);
        if(custo_viz<=custo) {
            substitui(sol, nova_sol, vert); // sol fica = nova_sol
            custo = custo_viz;
        } else {
            erro = custo_viz - custo; //tem que ser modulo
            // Aceitar com determinada probabiçlidade - calcular probabilidade usando funcao exp()
            /*  if erro > 0 then current <- $next 
                else current <- next with probibility e^(erro/temp) */

            erro = custo - custo_viz; 
            prob_aceitar = exp(erro/temp); // mutação
            if (prob_aceitar>rand_01()){
                substitui(sol, nova_sol, vert); 
                custo = custo_viz;
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

int verifica_validade(int *sol, int *mat,int vert) {
    int i = 0, j = 0;
    for (i = 0; i < vert; i++)
        if (sol[i] == 1) {
            for (j = 0; j < vert; j++)
                if (sol[j] == 1 && mat[i * vert + j] == 1) {
                    return 0;
                }
        }
    return 1;
}

int trepa_colinas(int sol[], int *mat, int vert, int num_iter)
{
    int *nova_sol, custo, custo_viz, i;
    int verifica=0;
    //matriz para a nova solucao
    nova_sol = malloc(sizeof(int)*vert);
    if(nova_sol == NULL){ printf("Erro na alocacao de memoria"); exit(1);}
    //------------------------------------------//

    // Avalia solucao inicial
    custo= calcula_fit(sol,mat,vert);
    //-------------------------------------//

    for(i=0; i<num_iter; i++)
    {
        do {
            gera_vizinho(sol, nova_sol, vert);
            verifica = verifica_validade(nova_sol, mat, vert);
        }while(verifica==0);
        // Avalia vizinho
        custo_viz= calcula_fit(nova_sol,mat,vert);

        if(custo_viz >= custo) //trocar isto pela linha de baixo
        {
            substitui(sol, nova_sol, vert);
            custo = custo_viz;
        }
    }
    free(nova_sol);
    return custo;
}

/*
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
        
    }
    free(nova_sol);
    free(nova_sol2);
    return custo;
}
*/
/*  /$$$$$$$$ /$$   /$$ /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$ 
   | $$_____/| $$  | $$| $$$ | $$ /$$__  $$ /$$__  $$ /$$__  $$
   | $$      | $$  | $$| $$$$| $$| $$  \__/| $$  \ $$| $$  \ $$
   | $$$$$   | $$  | $$| $$ $$ $$| $$      | $$$$$$$$| $$  | $$
   | $$__/   | $$  | $$| $$  $$$$| $$      | $$__  $$| $$  | $$
   | $$      | $$  | $$| $$\  $$$| $$    $$| $$  | $$| $$  | $$
   | $$      |  $$$$$$/| $$ \  $$|  $$$$$$/| $$  | $$|  $$$$$$/
   |__/       \______/ |__/  \__/ \______/ |__/  |__/ \______/  */

int calcula_fit(int sol[], int *mat, int vert) // mat é uma tabela de dimensões vert*vert
{
    // calcula o custo: 
    /*
    O custo da bissecção é dado pelo número de arcos que efetuam ligações entre
    vértices que pertencem a conjuntos diferentes (ou seja, arcos que ligam um vértice
    de V1 a um vértice de V2);
    */
	int custo=0;
	int i, j;

    //syntax for mat[l][c] will be mat[l*sizeY+c]

    for(i=0; i<vert; i++)
		if(sol[i]==1)
			custo++;

	return custo;
}  

int flip()
{
	if ((((float)rand()) / (float)RAND_MAX) < 0.5)
		return 0;
	else
		return 1;
}