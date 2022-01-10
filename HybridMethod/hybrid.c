#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define DEFAULT_RUNS	10
#define MAX_OBJ 1000
#define PROB 0.01

#define PRINT_0_1 0

int global = 0;

// EStrutura para armazenar parametros
struct info
{
    // Tamanho da popula��o
    int     popsize; // número de soluções
    // Probabilidade de muta��o
    float   pm;
    // Probabilidade de recombina��o
    float   pr;
    // Tamanho do torneio para sele��o do pai da pr�xima gera��o
	int     tsize;
	// Constante para avalia��o com penaliza��o
	float   ro; //
	// N�mero de objetos que se podem colocar na mochila
    int     numGenes; // obj // numero máximo de objetos das nossas sol // = número de vertices
	// Capacidade da mochila
	int     capacity; // cap
	// N�mero de gera��es
    int     numGenerations;  //max_gen
};

// Individuo (solu��o)
typedef struct individual chrom, *pchrom;

struct individual
{
    // Solu��o (objetos que est�o dentro da mochila)
    int     p[MAX_OBJ]; // tem 1 se o objeto pertence à solução, zero se não
    // Valor da qualidade da solu��o
	float   fitness;
    // 1 se for uma solu��o v�lida e 0 se n�o for
	int     valido;
};

// algoritmo
void tournament(pchrom pop, struct info d, pchrom parents);
void genetic_operators(pchrom parents, struct info d, pchrom offspring);
void crossover(pchrom parents, struct info d, pchrom offspring);
void mutation(pchrom offspring,struct info d);
void printsol(int * sol, int vert);

// utils
int * read_file(char * nome_fich,struct info * EA_param);
pchrom init_pop(struct info d);
void print_pop(pchrom pop, struct info d);
chrom get_best(pchrom pop, struct info d, chrom best);
void write_best(chrom x, struct info d);
void init_rand();
int random_l_h(int min, int max);
float rand_01();
int flip();

int solucaovalida(int *sol, int v, int * mat);
int calcula_fit(int sol[], int *mat, int vert);
int verifica_validade(int *sol, int *mat,struct info d);
int verifica_validade_l(int *sol, int *mat,int vert);




// ============= Local search functions ============= 
// funcao
void evaluate(pchrom pop, struct info d, int * mat);
void gera_vizinho(int sol[], int nova_sol[], int vert,int * mat);

void escreve_sol2(int *sol, int vert);
// algoritmo
int trepa_colinas(int sol[], int *mat, int vert, int num_iter);
int esferiamento(int sol[], int *mat, int vert, int num_iter, int print);

// utils
void gera_sol_inicial(int *sol, int v, int * mat, int prints);
void escreve_sol(int *sol, int vert);
void substitui(int a[], int b[], int n);
void init_rand(void);
int random_l_h(int min, int max);
float rand_01(void);
void printsol(int * sol, int vert);
int penalizacao(int *sol, int *mat,int vert, int fit_antigo);



/*
• Tamanho da população;
• Probabilidade de mutação;
• Probabilidade de recombinação;
• Tamanho do torneio;
• Número de gerações;
• Número de objetos do problema a otimizar;
• Capacidade da mochila;
*/

/*
pop: 100 população inicial
pm: 0.01 probabilidade de mutação
pr: 0.3 probabilidade de cruzamento
tsize: 2 tamanho do torneio
max_gen: 2500 máximo de gerações
obj: 100 número de objetos
cap: 250 capacidade da mochila (peso)
Weight Profit peso e lucro de cada objeto
2 8
5 1
10 5
*/

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


/* $$\      $$\  $$$$$$\  $$$$$$\ $$\   $$\ 
  $$$\    $$$ |$$  __$$\ \_$$  _|$$$\  $$ |
  $$$$\  $$$$ |$$ /  $$ |  $$ |  $$$$\ $$ |
  $$\$$\$$ $$ |$$$$$$$$ |  $$ |  $$ $$\$$ |
  $$ \$$$  $$ |$$  __$$ |  $$ |  $$ \$$$$ |
  $$ |\$  /$$ |$$ |  $$ |  $$ |  $$ |\$$$ |
  $$ | \_/ $$ |$$ |  $$ |$$$$$$\ $$ | \$$ |
  \__|     \__|\__|  \__|\______|\__|  \__| */
                                         


int main(int argc, char *argv[])
{
	char        nome_fich[100] = "", tmp[100] = "";
	struct info EA_param;
	pchrom      pop = NULL, parents = NULL;
	chrom       best_run, best_ever;
	int         gen_actual, r, runs, i, inv, *mat;
	float       mbf = 0.0;
    int custo, best_custo = 0 ;

    // L� os argumentos de entrada
	if (argc == 3)
	{
		runs = atoi(argv[2]);
		strcpy(nome_fich, argv[1]);
	}
	else
        // Se o n�mero de execu��es do processo n�o for colocado nos argumentos de entrada, define-o com um valor por defeito
        if (argc == 2)
        {
            runs = DEFAULT_RUNS;
            strcpy(nome_fich, argv[1]);
        }
        // Se o nome do ficheiro de informa��es n�o for colocado nos argumentos de entrada, pede-o novamente
        else
        {
            runs = DEFAULT_RUNS;
            printf("Nome do Ficheiro: ");
            fgets(tmp,100,stdin);
			strtok(tmp, "\n"); // para retirar o '\n' do final
        }
	// adicionar caminho ao nome do ficheiro
	strcat(nome_fich,"./Instancias/");
	//strcat(nome_fich,"./EvolutionaryAlgorithm/refe/");
	strcat(nome_fich,tmp);
	printf("%s",nome_fich);

    // Se o n�mero de execu��es do processo for menor ou igual a 0, termina o programa
	if (runs <= 0)
		return 0;
    //Inicializa a gera��o dos n�meros aleat�rios
	init_rand();
    // Preenche a matriz com dados dos objectos (peso e valor) e a estrutura EA_param que foram definidos no ficheiro de input
	mat = read_file(nome_fich,&EA_param);
    int half_the_iterations = EA_param.numGenerations;
    int prints = PRINT_0_1;
	// Faz um ciclo com o n�mero de execu��es definidas

    // local algorthm
    int vert;
    int num_iter;
    int * sol;
    int * best;

	for (r=0; r<runs; r++)
	{
        // Gera��o da popula��o inicial

		pop = init_pop(EA_param); // pop ficam com array de soluções (structs chrom), que dentro delas têm arrays de 1s e zeros aleatorios

        // Avalia a popula��o inicial
		evaluate(pop, EA_param, mat); // mat tem os dados direitinhos do ficheiro de texto

		// Como ainda n�o existe, escolhe-se como melhor solu��o a primeira da popula��o (poderia ser outra qualquer)
		best_run = pop[0];
        // Encontra-se a melhor solu��o dentro de toda a popula��o
		best_run = get_best(pop, EA_param, best_run);
        // Reserva espa�o para os pais da popula��o seguinte
		parents = malloc(sizeof(chrom)*EA_param.popsize);
		if (parents==NULL){printf("Erro na alocacao de memoria\n");exit(1);}

		// Ciclo de optimiza��o
		gen_actual = 1;
        /*
        while (gen_actual <= EA_param.numGenerations)
		{

            // Torneio bin�rio para encontrar os progenitores (ficam armazenados no vector parents)
			tournament(pop, EA_param, parents); // parents fica com uma série de soluções um pouco melhores que as anteriores, vai haver soluções repetidas, vão se perder algumas das soluções piores
            // Aplica os operadores gen�ticos aos pais (os descendentes ficam armazenados na estrutura pop)
			genetic_operators(parents, EA_param, pop); // pop fica com o crossover dos parents e da propria pop, e é aplicada mutação à pop
            // Avalia a nova popula��o (a dos filhos)
			evaluate(pop, EA_param, mat); // calcula se é válida e o fitness de cada solução
            // Actualiza a melhor solu��o encontrada
			best_run = get_best(pop, EA_param, best_run);
			gen_actual++;

			//write_best(*pop, EA_param);
			//printf("validade:%d\n",pop->valido);

		}
        */


        vert = EA_param.numGenes;
        num_iter = EA_param.numGenerations;
        sol = malloc(sizeof(int)*EA_param.numGenes);
        best = malloc(sizeof(int)*EA_param.numGenes);
        for(int z=0; z<EA_param.numGenes; z++){ //meter tudo a zeros no array sol de tamanho v(vértices)
            sol[z]=0;
            best[z]=0;
            pop->p[z] = 0;
            best_run.p[z] = 0;
        }
        if(sol == NULL || best == NULL){ printf("Erro na alocacao de memoria"); exit(1);}
        for (int z = 0; z<EA_param.numGenes; z++) {
            sol[z] = pop->p[z];
        }

        for (int z = 0; z<EA_param.numGenes; z++) {
            best[z] = best_run.p[z];
        }

        // Gerar solucao inicial
        
        gera_sol_inicial(sol, vert, mat, prints);
        if (prints) escreve_sol(sol, vert); // ***
        // Trepa colinas
        //custo = trepa_colinas(sol, grafo, vert, num_iter);

        // ============================ TREPA COLINAS ========================== //
        //custo = trepa_colinas(sol, mat, vert, num_iter);
        //int trepa_colinas(int sol_t[], int *mat_t, int vert_t, int num_iter_t)

        //int *sol_t = sol;
        //int *mat_t = mat;
        //int vert_t = vert;
        //int num_iter_t = num_iter;
        int *nova_sol_t, custo_t=0, fit_viz_t, i_t;
        int verifica_t=0;
        //matriz para a nova solucao
        nova_sol_t = malloc(sizeof(int)*vert);
        if(nova_sol_t == NULL){ printf("Erro na alocacao de memoria"); exit(1);}
        //------------------------------------------//

        // Avalia solucao inicial
        custo_t= calcula_fit(sol,mat,vert);
        //-------------------------------------//

        for(i_t=0; i_t<num_iter; i_t++)
        {

            if (i_t%2 == 0){

                //printf("making local\n");
                do {
                    gera_vizinho(sol, nova_sol_t, vert, mat);
                    verifica_t = verifica_validade_l(nova_sol_t, mat, vert);
                    //verifica = 1; // PENALIZAÇAO
                }while(verifica_t==0);
                // Avalia vizinho
                fit_viz_t= calcula_fit(nova_sol_t,mat,vert);

                if(fit_viz_t >= custo_t){ // PARA ACEITAR SOLUÇÕES DE CUSTO IGUAL OU NAO
                    substitui(sol, nova_sol_t, vert);
                    custo_t = fit_viz_t;
                }

                // ========= tradução =========== //
                for (int z = 0; z<EA_param.numGenes; z++) {
                    pop->p[z] = sol[z];
                }

                for (int z = 0; z<EA_param.numGenes; z++) {
                    best_run.p[z] = best[z];
                }
                // ========= tradução =========== //

            } else {

                //printf("making evolutionary\n");
                // Torneio bin�rio para encontrar os progenitores (ficam armazenados no vector parents)
                tournament(pop, EA_param, parents); // parents fica com uma série de soluções um pouco melhores que as anteriores, vai haver soluções repetidas, vão se perder algumas das soluções piores
                // Aplica os operadores gen�ticos aos pais (os descendentes ficam armazenados na estrutura pop)
                genetic_operators(parents, EA_param, pop); // pop fica com o crossover dos parents e da propria pop, e é aplicada mutação à pop
                // Avalia a nova popula��o (a dos filhos)
                evaluate(pop, EA_param, mat); // calcula se é válida e o fitness de cada solução
                // Actualiza a melhor solu��o encontrada
                best_run = get_best(pop, EA_param, best_run);
                

                // ========= tradução =========== //
                for (int z = 0; z<EA_param.numGenes; z++) {
                    sol[z] = pop->p[z];
                }

                for (int z = 0; z<EA_param.numGenes; z++) {
                    best[z] = best_run.p[z];
                }
                // ========= tradução =========== //
            }

        }
        free(nova_sol_t);

        custo = custo_t;

        // ============================ TREPA COLINAS ========================== //


/*
		while (gen_actual <= EA_param.numGenerations)
		{

            // Torneio bin�rio para encontrar os progenitores (ficam armazenados no vector parents)
			tournament(pop, EA_param, parents); // parents fica com uma série de soluções um pouco melhores que as anteriores, vai haver soluções repetidas, vão se perder algumas das soluções piores
            // Aplica os operadores gen�ticos aos pais (os descendentes ficam armazenados na estrutura pop)
			genetic_operators(parents, EA_param, pop); // pop fica com o crossover dos parents e da propria pop, e é aplicada mutação à pop
            // Avalia a nova popula��o (a dos filhos)
			evaluate(pop, EA_param, mat); // calcula se é válida e o fitness de cada solução
            // Actualiza a melhor solu��o encontrada
			best_run = get_best(pop, EA_param, best_run);
			gen_actual++;

			//write_best(*pop, EA_param);
			//printf("validade:%d\n",pop->valido);

		}
*/

        //custo = trepa_colinas(sol, mat, vert, num_iter);
        //custo = esferiamento(sol, mat, vert, num_iter, prints); //esfriamento // sol fica com a melhor soulução e retorna o custo dela
        // Escreve resultados da repeticao k
        printf("\nRepeticao %d:", r); // ***
        if (prints) escreve_sol(sol, vert); // ***
        if (prints) printf("Custo final: %2d\n", custo); // ***
        mbf += custo;
        if(r==0 || best_custo > custo)
        {
            best_custo = custo;
            substitui(best, sol, vert);
        }

		// Contagem das solu��es inv�lidas
		for (inv=0, i=0; i<EA_param.popsize; i++)
			if (pop[i].valido == 0){
				inv++;
				//printf("sol %d invalida\n", i);
			} else {
				//printf("sol %d VALIDA\n", i);
			}
		// Escreve resultados da repeti��o que terminou
		//printf("\nRepeticao %d:", r);
		//write_best(best_run, EA_param);
		//printf("\nPercentagem Invalidos: %f\n", 100*(float)inv/EA_param.popsize);
		mbf += best_run.fitness;
		if (r==0 || best_run.fitness > best_ever.fitness)
			best_ever = best_run;
        // Liberta a mem�ria usada
		free(parents);
		free(pop);
	}
    free(mat);
    free(sol);
    free(best);
	// Escreve resultados globais
	printf("\n\nMBF(mean best fitness): %f\n", mbf/r); 
	printf("\nMelhor solucao encontrada");
	printf("\nvalida?: %d\n",best_ever.valido);
	write_best(best_ever, EA_param);

/*
    printf("\n========LOCAL resultados globais=====\n"); 
    printf("\nMBF(mean best feat): %f\n", vert-(mbf/r)); //médio
    printf("Quantos verts em media das runs: %f\n",mbf/r );
    printf("Custo final (melhor): %2d \n", vert-best_custo); // (quanto menor melhor, = número de verts que nao pertencem a sol)
    printf("Número de vertices da melhor sol: %2d ", best_custo);
    printf("|||%d||", vert);
    printf("\n====================\n");
    */
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

// Inicializa��o do gerador de n�meros aleat�rios
void init_rand()
{
	srand((unsigned)time(NULL));
}

// Leitura dos par�metros e dos dados do problema
// Par�metros de entrada: Nome do ficheiro e matriz a preencher com os dados dos objectos (peso e valor)
// Par�metros de sa�da: Devolve a estrutura com os par�metros
int * read_file(char *filename, struct info * pEA_param )
{
	
	pEA_param->popsize = 200; //fscanf(f, " pop: %d", &x.popsize);
	pEA_param->pm = 0.01; //fscanf(f, " pm: %f", &x.pm);
	pEA_param->pr = 0.7; //fscanf(f, " pr: %f", &x.pr);
	pEA_param->tsize = 2; //fscanf(f, " tsize: %d", &x.tsize);
	pEA_param->numGenerations = 700;//fscanf(f, " max_gen: %d", &x.numGenerations); //max_gen
	//x.capacity = 250;//fscanf(f, " cap: %d", &x.capacity);
	pEA_param->ro = 0.0;
	if (pEA_param->numGenes > MAX_OBJ){printf("Number of itens is superior to MAX_OBJ\n");exit(1);}

	FILE *f;
	int *p, *q;
	int i, j;
    int c;
	int arestas;

	// Numero de iteracoes
	//printf("Escolha o numero de iteracoes: ");
    //scanf("%d",num_iter);
	//while ((c = getchar()) != '\n' && c != EOF) { }

	f=fopen(filename, "rt");
	if(!f) { printf("Erro no acesso ao ficheiro dos dados\n"); exit(1); }

    // salta as linhas de comentario
	char ded;
	fscanf(f," %c", &ded);
	while (ded == 'c') {
		fscanf(f," %*[^\n]");
		fscanf(f," %c", &ded);
	}

    if (ded == 'p'){  //  <vertices> <arestas>
		fscanf(f," edge %d %d", &pEA_param->numGenes, &arestas);
	}

	// Alocacao dinamica da matriz
	p = malloc(sizeof(int)*(pEA_param->numGenes)*(pEA_param->numGenes));
	if(!p){ printf("Erro na alocacao de memoria\n");exit(1);}

    memset(p, 0, sizeof(int)*(pEA_param->numGenes)*(pEA_param->numGenes));

	q=p;
    // Preenchimento da matriz unidimensional, tratada como bidimensional, cada linha representa um vertice e tem uma coluna para cada outro verticer, onde 0 é que há uma aresta entre eles, e 1 quer dizer que não há
    int vert1, vert2;
    for (i=0; i< arestas; i++){
		fscanf(f, " e %d %d", &vert1, &vert2);
        p[(vert1-1)*(pEA_param->numGenes)+(vert2-1)]=1;
		p[(vert2-1)*(pEA_param->numGenes)+(vert1-1)]=1;
	}

	fclose(f);

    // print da matriz
    int k=0;
    for (i=0; i<pEA_param->numGenes; i++) {
        for(j=0; j<pEA_param->numGenes; j++) {
            printf("%d ",q[k++]);
        }
        putchar('\n');
    }

	return p;

}

// Simula o lan�amento de uma moeda, retornando o valor 0 ou 1
int flip()
{
	if ((((float)rand()) / (float)RAND_MAX) < 0.5)
		return 0;
	else
		return 1;
}


// Criacao da populacao inicial. O vector e alocado dinamicamente
// Par�metro de entrada: Estrutura com par�metros do problema
// Par�metro de sa�da: Preenche da estrutura da popula��o apenas o vector bin�rio com os elementos que est�o dentro ou fora da mochila
pchrom init_pop(struct info d)
{

	int     i, j,x;
	pchrom  indiv;

	indiv = malloc(sizeof(chrom)*d.popsize);
	if (indiv==NULL){printf("Erro na alocacao de memoria\n");exit(1);}



        for(i=0; i<d.popsize; i++) //meter tudo a zeros no array sol de tamanho v(vértices)
            for(j=0; j<d.numGenes; j++)
				indiv[i].p[j] = 0;

		i = j = 0;

        for(i=0; i<d.popsize; i++){
			int num_verts_sol = random_l_h(0, d.numGenes-1);
			for(j=0; j<num_verts_sol; j++){
				do{
					x = random_l_h(0, d.numGenes-1);
				}while(indiv[i].p[x] != 0);
				indiv[i].p[x]=1;
			}
		}

	return indiv;
}

void printsol(int * sol, int vert){
    for(int i=0; i<vert; i++)
        printf("%d ", sol[i]);
    putchar('|');
}

// Actualiza a melhor solu��o encontrada
// Par�metro de entrada: populacao actual (pop), estrutura com par�metros (d) e a melhor solucao encontrada at� a gera��oo imediatamente anterior (best)
// Par�metro de sa�da: a melhor solucao encontrada at� a gera��o actual
chrom get_best(pchrom pop, struct info d, chrom best)
{
	int i;

	for (i=0; i<d.popsize; i++) //popsize é o tamanho do array pop
	{
		if (best.fitness < pop[i].fitness && pop[i].valido == 1)
			best=pop[i];
	}
	return best;
}

// Devolve um valor inteiro distribuido uniformemente entre min e max
int random_l_h(int min, int max)
{
	return min + rand() % (max-min+1);
}

// Devolve um valor real distribuido uniformemente entre 0 e 1
float rand_01()
{
	return ((float)rand())/(float)RAND_MAX;
}

// Escreve uma solu��o na consola
// Par�metro de entrada: populacao actual (pop) e estrutura com par�metros (d)
void write_best(chrom x, struct info d)
{
	int i;

	printf("Best individual: %4.1f\n", x.fitness);
	for (i=0; i<d.numGenes; i++)
		printf("%d", x.p[i]);
	putchar('\n');
}
                                   

/*   /$$$$$$  /$$        /$$$$$$   /$$$$$$  /$$$$$$$  /$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$ 
    /$$__  $$| $$       /$$__  $$ /$$__  $$| $$__  $$|_  $$_/|__  $$__/| $$$    /$$$ /$$__  $$
   | $$  \ $$| $$      | $$  \__/| $$  \ $$| $$  \ $$  | $$     | $$   | $$$$  /$$$$| $$  \ $$
   | $$$$$$$$| $$      | $$ /$$$$| $$  | $$| $$$$$$$/  | $$     | $$   | $$ $$/$$ $$| $$  | $$
   | $$__  $$| $$      | $$|_  $$| $$  | $$| $$__  $$  | $$     | $$   | $$  $$$| $$| $$  | $$
   | $$  | $$| $$      | $$  \ $$| $$  | $$| $$  \ $$  | $$     | $$   | $$\  $ | $$| $$  | $$
   | $$  | $$| $$$$$$$$|  $$$$$$/|  $$$$$$/| $$  | $$ /$$$$$$   | $$   | $$ \/  | $$|  $$$$$$/
   |__/  |__/|________/ \______/  \______/ |__/  |__/|______/   |__/   |__/     |__/ \______/  */


// Preenche uma estrutura com os progenitores da pr�xima gera��o, de acordo com o resultados do torneio binario (tamanho de torneio: 2)
// Par�metros de entrada: popula��o actual (pop), estrutura com par�metros (d) e popula��o de pais a encher
void tournament(pchrom pop, struct info d, pchrom parents){
	int i, x1, x2;

	// Realiza popsize torneios
	for(i=0; i<d.popsize;i++) // tamanho do array pop com as soluções
	{
		x1 = random_l_h(0, d.popsize-1);
		do
			x2 = random_l_h(0, d.popsize-1);
		while(x1==x2);
		if(pop[x1].fitness > pop[x2].fitness)		// Problema de maximizacao
			parents[i]=pop[x1];
		else
			parents[i]=pop[x2];
	}
}

// Operadores geneticos a usar na gera��o dos filhos
// Par�metros de entrada: estrutura com os pais (parents), estrutura com par�metros (d), estrutura que guardar� os descendentes (offspring)
void genetic_operators(pchrom parents, struct info d, pchrom pop)
{
/*			printf("before crossover\n");
			for (int l=0; l<d.popsize; l++) {
				printsol(pop[l].p, d.numGenes);
				//putchar('\n');
			}
			printf("\n=========================\n");
*/    // Recombina��o com um ponto de corte
	crossover(parents, d, pop); // quando sucesso na probabilidade de recombinação, pop fica com metade dos bits do parent e mantém metade dos seus
/*			printf("after crossover\n");
			for (int l=0; l<d.popsize; l++) {
				printsol(pop[l].p, d.numGenes);
				//putchar('\n');
			}
			printf("\n=========================\n");
*/	// Muta��o bin�ria
	mutation(pop, d); // para cada bit de cada solução vê atravéz da probabilidade de mutação se lhe dá flip
	/*
			printf("after mutation\n");
			for (int l=0; l<d.popsize; l++) {
				printsol(pop[l].p, d.numGenes);
				//putchar('\n');
			}
			printf("\n=========================\n");
*/
			//static int i;
			//i++;
			//if (i==10) exit(0);
}

// Preenche o vector descendentes com o resultado das opera��es de recombina��o
// Par�metros de entrada: estrutura com os pais (parents), estrutura com par�metros (d), estrutura que guardar� os descendentes (offspring)
void crossover(pchrom parents, struct info d, pchrom offspring)
{
	int i, j, point;

	for (i=0; i<d.popsize; i+=2) // tam de parents
	{
		if (rand_01() < d.pr) // rand_01 retorna um valor entre 0 e 1, pr = probabilidade de recombinação
		{
			point = random_l_h(0, d.numGenes-1); // num genes é o tamanho de obj para a mochila
			for (j=0; j<point; j++)
			{
				offspring[i].p[j] = parents[i].p[j];
				offspring[i+1].p[j] = parents[i+1].p[j];
			}
			for (j=point; j<d.numGenes; j++)
			{
				offspring[i].p[j]= parents[i+1].p[j];
				offspring[i+1].p[j] = parents[i].p[j];
			}
		}else{
			offspring[i] = parents[i];
			offspring[i+1] = parents[i+1];
		}
	}
}

// Muta��o bin�ria com v�rios pontos de muta��o
// Par�metros de entrada: estrutura com os descendentes (offspring) e estrutura com par�metros (d)
void mutation(pchrom offspring, struct info d)
{
	int i, j;

	for (i=0; i<d.popsize; i++)
		for (j=0; j<d.numGenes; j++)
			if (rand_01() < d.pm)
				offspring[i].p[j] = !(offspring[i].p[j]);
}


/*  /$$$$$$$$ /$$   /$$ /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$ 
   | $$_____/| $$  | $$| $$$ | $$ /$$__  $$ /$$__  $$ /$$__  $$
   | $$      | $$  | $$| $$$$| $$| $$  \__/| $$  \ $$| $$  \ $$
   | $$$$$   | $$  | $$| $$ $$ $$| $$      | $$$$$$$$| $$  | $$
   | $$__/   | $$  | $$| $$  $$$$| $$      | $$__  $$| $$  | $$
   | $$      | $$  | $$| $$\  $$$| $$    $$| $$  | $$| $$  | $$
   | $$      |  $$$$$$/| $$ \  $$|  $$$$$$/| $$  | $$|  $$$$$$/
   |__/       \______/ |__/  \__/ \______/ |__/  |__/ \______/  */


/* CERQUEIRA
float eval_individual(int sol[], struct info d, int *mat, int *v){
	int total=0;
    int i ,verifica;
    verifica=verifica_validade(sol, mat,d);
    if(verifica==0) {
        *v=0;
        return 0;
    }else{
    	for(i=0;i<d.numGenes;i++){
        	if(sol[i]==1)
            	total++;
    		}
    	*v=1;
    	return total;
	}
}

int verifica_validade(int *sol, int *mat,struct info d) {
    int i = 0, j = 0;
    for (i = 0; i < d.numGenes; i++)
        if (sol[i] == 1) {
            for (j = 0; j < d.numGenes; j++)
                if (sol[j] == 1 && mat[i * d.numGenes + j] == 1) {
                    return 0;
                }
        }
    return 1;
}
*/

// Calcula a qualidade de uma solu��o
// Par�metros de entrada: solu��o (sol), capacidade da mochila (d), matriz com dados do problema (mat) e numero de objectos (v)
// Par�metros de sa�da: qualidade da solu��o (se a capacidade for excedida devolve 0)

float eval_individual(int sol[], struct info d, int *mat, int *v){
	int     i, min;
	float   ro;

		if (solucaovalida(sol, d.numGenes, mat) == 0) { // solução inválida
			//*v = 0;
			/*
			for (i=0; sol[i] == 0 && i<d.numGenes; i++) // Procurar primeiro objeto na mochila
				; 
			min = i;
			for (; i<d.numGenes; i++) // Procura objeto menos valioso na mochila
				if(sol[i]==1 && mat[i][1] < mat[min][1])
					min = i;
			sol[min] = 0; // retirá-lo
			*/
			
				int x=0;
				int num_verts_sol = random_l_h(0, d.numGenes-1);
				for(int j=0; j<d.numGenes; j++){
					do{
						x = random_l_h(0, d.numGenes-1);
					}while(sol[x] == 0);
					sol[x]=1;
				}
				*v =0;
			
		} else {
			*v = 1;
		}

	int fit = calcula_fit(sol,mat,d.numGenes);
	if (*v == 0 ) fit = 0;
	// v é a validade da solução (0 - inválido, 1 - válido)
	return fit;
}

int calcula_fit(int sol[], int *mat, int vert) // mat é uma tabela de dimensões vert*vert
{
    // calcula o custo: 
    /*
    O custo da bissecção é dado pelo número de arcos que efetuam ligações entre
    vértices que pertencem a conjuntos diferentes (ou seja, arcos que ligam um vértice
    de V1 a um vértice de V2);
    */
	int qualidade=0;
	int i, j;

    //syntax for mat[l][c] will be mat[l*sizeY+c]

    for(i=0; i<vert; i++)
		if(sol[i]==1)
			qualidade++;

	return qualidade;
}  

int solucaovalida(int *sol, int v, int * mat){
    int i=0, j=0;
    for(i=0; i<v; i++)
		if(sol[i]==1)
		{
			for(j=0; j<v;j++)
				if(sol[j]==1 && mat[i*v+j]==1){
				    //printf(" sol inval ");
                    return 0; 
                }
		}
    return 1;
}


// Avaliacao da popula��o
// Par�metros de entrada: populacao (pop), estrutura com parametros (d) e matriz com dados do problema (mat)
// Par�metros de sa�da: Preenche pop com os valores de fitness e de validade para cada solu��o
void evaluate(pchrom pop, struct info d, int * mat)
{
	int i;

	for (i=0; i<d.popsize; i++)
		pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido);
}

 //============================ ====================================================================================
/*
 $$        /$$$$$$   /$$$$$$   /$$$$$$  /$$      
| $$       /$$__  $$ /$$__  $$ /$$__  $$| $$      
| $$      | $$  \ $$| $$  \__/| $$  \ $$| $$      
| $$      | $$  | $$| $$      | $$$$$$$$| $$      
| $$      | $$  | $$| $$      | $$__  $$| $$      
| $$      | $$  | $$| $$    $$| $$  | $$| $$      
| $$$$$$$$|  $$$$$$/|  $$$$$$/| $$  | $$| $$$$$$$$
|________/ \______/  \______/ |__/  |__/|________/
*/      
//============================ ====================================================================================
                                            

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

    } while (solucaovalida(sol, v, mat) == 0);


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
void gera_vizinho(int sol[], int nova_sol[], int vert, int * mat)
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
    if (tudoAZeros(sol,vert) == 0){
    do {
        if (flip()==1){ 
            if (tudoAZeros(sol,vert) == 0){
                do{
                    if (tudoAZeros(nova_sol,vert) == 1) break;
                    p2=random_l_h(0, vert-1);
                }while(nova_sol[p2] != 1);
            }
            nova_sol[p2] = 0;
        }

    }while(verifica_validade_l(nova_sol, mat, vert) == 0);
	// Troca
    }
    
    
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
            gera_vizinho(sol, nova_sol, vert, mat);
        } while (solucaovalida(nova_sol, vert, mat) == 0);
		
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


//int verifica_validade(int *sol, int *mat,struct info d);
//local algorythm:
int verifica_validade_l(int *sol, int *mat,int vert) {
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

int verifica_validade(int *sol, int *mat,struct info d) {
    int i = 0, j = 0;
    for (i = 0; i < d.numGenes; i++)
        if (sol[i] == 1) {
            for (j = 0; j < d.numGenes; j++)
                if (sol[j] == 1 && mat[i * d.numGenes + j] == 1) {
                    return 0;
                }
        }
    return 1;
}

int trepa_colinas(int sol_t[], int *mat_t, int vert_t, int num_iter_t)
{
    int *nova_sol_t, custo_t=0, fit_viz_t, i_t;
    int verifica_t=0;
    //matriz para a nova solucao
    nova_sol_t = malloc(sizeof(int)*vert_t);
    if(nova_sol_t == NULL){ printf("Erro na alocacao de memoria"); exit(1);}
    //------------------------------------------//

    // Avalia solucao inicial
    custo_t= calcula_fit(sol_t,mat_t,vert_t);
    //-------------------------------------//

    for(i_t=0; i_t<num_iter_t; i_t++)
    {
        do {
            gera_vizinho(sol_t, nova_sol_t, vert_t, mat_t);
            verifica_t = verifica_validade_l(nova_sol_t, mat_t, vert_t);
            //verifica = 1; // PENALIZAÇAO
        }while(verifica_t==0);
        // Avalia vizinho
        fit_viz_t= calcula_fit(nova_sol_t,mat_t,vert_t);

        if(fit_viz_t >= custo_t){ // PARA ACEITAR SOLUÇÕES DE CUSTO IGUAL OU NAO
            substitui(sol_t, nova_sol_t, vert_t);
            custo_t = fit_viz_t;
        }
    }
    free(nova_sol_t);
    return custo_t;
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
