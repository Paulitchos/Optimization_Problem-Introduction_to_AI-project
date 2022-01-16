#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define DEFAULT_RUNS	50
#define MAX_OBJ 1000

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
float eval_individual(int sol[], struct info d, int *mat, int *v, float best_fitness_found);


// funcao
void evaluate(pchrom pop, struct info d, int * mat, float best_fitness_found);


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
	// Faz um ciclo com o n�mero de execu��es definidas

	float best_fitness_found = 0;

	for (r=0; r<runs; r++)
	{
        // Gera��o da popula��o inicial
		pop = init_pop(EA_param); // pop ficam com array de soluções (structs chrom), que dentro delas têm arrays de 1s e zeros aleatorios

        // Avalia a popula��o inicial
		evaluate(pop, EA_param, mat,  best_fitness_found); // mat tem os dados direitinhos do ficheiro de texto

		// Como ainda n�o existe, escolhe-se como melhor solu��o a primeira da popula��o (poderia ser outra qualquer)
		best_run = pop[0];
        // Encontra-se a melhor solu��o dentro de toda a popula��o
		best_run = get_best(pop, EA_param, best_run);
        // Reserva espa�o para os pais da popula��o seguinte
		parents = malloc(sizeof(chrom)*EA_param.popsize);
		if (parents==NULL){printf("Erro na alocacao de memoria\n");exit(1);}

		
		// Ciclo de optimiza��o
		gen_actual = 1;
		while (gen_actual <= EA_param.numGenerations)
		{
            // Torneio bin�rio para encontrar os progenitores (ficam armazenados no vector parents)
			tournament(pop, EA_param, parents); // parents fica com uma série de soluções um pouco melhores que as anteriores, vai haver soluções repetidas, vão se perder algumas das soluções piores
            // Aplica os operadores gen�ticos aos pais (os descendentes ficam armazenados na estrutura pop)
			genetic_operators(parents, EA_param, pop); // pop fica com o crossover dos parents e da propria pop, e é aplicada mutação à pop
            // Avalia a nova popula��o (a dos filhos)
			evaluate(pop, EA_param, mat, best_run.fitness); // calcula se é válida e o fitness de cada solução
            // Actualiza a melhor solu��o encontrada
			best_run = get_best(pop, EA_param, best_run);
			gen_actual++;

		// write_best(*pop, EA_param);
		//	printf("validade:%d\n",pop->valido);

			//global ++;
			//if (global == 20) exit (0);
			//for (int l=0; l<EA_param.popsize; l++) {
			//	printsol(pop[l].p, EA_param.numGenes / 2);
			//	putchar('\n');
			//}
			//printf("\n=========================\n");

			//printsol(best_run.p, EA_param.numGenes);
			//printf("\n=========================\n");

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
		printf("\nRepeticao %d:", r);
		write_best(best_run, EA_param);
		printf("\nPercentagem Invalidos: %f\n", 100*(float)inv/EA_param.popsize);
		mbf += best_run.fitness;
		if (r==0 || best_run.fitness > best_ever.fitness)
			best_ever = best_run;
        // Liberta a mem�ria usada
		free(parents);
		free(pop);
	}
	// Escreve resultados globais
	printf("\n\nMBF(mean best fitness): %f\n", mbf/r); 
	printf("\nMelhor solucao encontrada");
	printf("\nvalida?: %d\n",best_ever.valido);
	write_best(best_ever, EA_param);
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
	
	pEA_param->popsize = 10; //fscanf(f, " pop: %d", &x.popsize);
	pEA_param->pm = 0.1; //fscanf(f, " pm: %f", &x.pm);
	pEA_param->pr = 0.3; //fscanf(f, " pr: %f", &x.pr);
	pEA_param->tsize = 2; //fscanf(f, " tsize: %d", &x.tsize);
	pEA_param->numGenerations = 1000;//fscanf(f, " max_gen: %d", &x.numGenerations); //max_gen
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
int flip(){
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
            
/*
        int num_verts_sol = random_l_h(0, v-1);
        for(i=0; i<num_verts_sol; i++){
            // mete aleatoriamente um dos indexs que estão a zero a um
            do
                x = random_l_h(0, v-1);
            while(sol[x] != 0);
            sol[x]=1;
        }
*/
/*
	for (i=0; i<d.popsize; i++)
	{
		for (j=0; j<d.numGenes; j++){ // numGenes = Número de vertices máximos da sol
			indiv[i].p[j] = flip();
		}
	}
	*/
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
	//for (i=0; i<d.numGenes; i++)
	//	printf("%d", x.p[i]);
	//putchar('\n');
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
	/*int i, j,num;

	for (i=0; i<d.popsize; i++)
		for (j=0; j<d.numGenes; j++)
			if (rand_01() < d.pm)
				/*do{
                    num = random_l_h(0,(d.numGenes-1));
                }while(num == offspring[i].p[j]);*/
				//offspring[i].p[j] = !(offspring[i].p[j]); 
	
	int i, j, g1, g2, p1, count;

	for (i = 0; i < d.popsize; i++) {
		if (rand_01() < d.pm) { // Mutar individuo
			
			do { // Garantir que o grupo escolhido tem pelo menos um elemento
				count = 0; 
				g1 = random_l_h(0, d.numGenes-1);

				for (j = 0; j < d.numGenes; j++) {
					if (offspring[i].p[j] == g1) {
						count++;
						break;
					}
				}
			} while (count == 0);

            do
                g2 = random_l_h(0, d.numGenes-1);
            while (g2 == g1);

            do
				p1 = random_l_h(0, d.numGenes-1);
			while (offspring[i].p[p1] != g1);

            offspring[i].p[p1] = g2;
        }
	}
}


/*  /$$$$$$$$ /$$   /$$ /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$ 
   | $$_____/| $$  | $$| $$$ | $$ /$$__  $$ /$$__  $$ /$$__  $$
   | $$      | $$  | $$| $$$$| $$| $$  \__/| $$  \ $$| $$  \ $$
   | $$$$$   | $$  | $$| $$ $$ $$| $$      | $$$$$$$$| $$  | $$
   | $$__/   | $$  | $$| $$  $$$$| $$      | $$__  $$| $$  | $$
   | $$      | $$  | $$| $$\  $$$| $$    $$| $$  | $$| $$  | $$
   | $$      |  $$$$$$/| $$ \  $$|  $$$$$$/| $$  | $$|  $$$$$$/
   |__/       \______/ |__/  \__/ \______/ |__/  |__/ \______/  */

int solucaovalida(int *sol, int v, int * mat);
int calcula_fit(int sol[], int *mat, int vert);
int verifica_validade(int *sol, int *mat,struct info d);

// Calcula a qualidade de uma solu��o
// Par�metros de entrada: solu��o (sol), capacidade da mochila (d), matriz com dados do problema (mat) e numero de objectos (v)
// Par�metros de sa�da: qualidade da solu��o (se a capacidade for excedida devolve 0)

int num_arestas_sol(int *sol, int v, int * mat){
	int num_arestas;
	int i=0, j=0;
    for(i=0; i<v; i++)
		if(sol[i]==1)
		{
			for(j=0; j<v;j++)
				if(sol[j]==1 && mat[i*v+j]==1){
				    //printf(" sol inval ");
                    num_arestas++; 
                }
		}
	return num_arestas;
}

float eval_individual(int sol[], struct info d, int *mat, int *v, float best_fitness_found){
	int     i, min;
	float   ro;

		if (solucaovalida(sol, d.numGenes, mat) == 0) { // solução inválida			
				int fit = calcula_fit(sol,mat,d.numGenes);
				int num_arestas_mal = num_arestas_sol(sol,d.numGenes,mat);
				*v =0;
				if ( fit <= num_arestas_mal) {
					fit = 0;
					return fit;
				}
				return fit  - num_arestas_mal;
			
		} else {
			*v = 1;
		}

	int fit = calcula_fit(sol,mat,d.numGenes);
	
	//if (*v == 0 ) fit = 0;
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
void evaluate(pchrom pop, struct info d, int * mat, float best_fitness_found)
{
	int i;

	for (i=0; i<d.popsize; i++)
		pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido, best_fitness_found);
}