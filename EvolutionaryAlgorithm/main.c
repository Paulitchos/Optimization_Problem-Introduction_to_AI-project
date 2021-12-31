#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DEFAULT_RUNS	10
#define MAX_OBJ 1000

// Individuo (solucao)
typedef struct individual chrom, *pchrom;
struct individual{
    // Solucao (objetos que estao dentro da mochila)
    int     p[MAX_OBJ];
    // Valor da qualidade da solucao
	float   fitness;
    // 1 se for uma solucao valida e 0 se nao for
	int     valido;
};

// EStrutura para armazenar parametros
struct info{
    // Tamanho da populacao
    int     popsize;
    // Probabilidade de mutacao
    float   pm;
    // Probabilidade de recombinacao
    float   pr;
    // Tamanho do torneio para selecao do pai da proxima geracao
	int     tsize;
	// Constante para avaliacao com penalizacao
	float   ro;
	// Numero de objetos que se podem colocar na mochila
    int     numGenes;
	// Capacidade da mochila
	int     capacity;
	// Numero de geracoes
    int     numGenerations;
};

// algoritmo
void tournament(pchrom pop, struct info d, pchrom parents);
void genetic_operators(pchrom parents, struct info d, pchrom offspring);
void crossover(pchrom parents, struct info d, pchrom offspring);
void mutation(pchrom offspring,struct info d);

// utils
struct info read_file(char *s, int mat[][2]);
pchrom init_pop(struct info d);
void print_pop(pchrom pop, struct info d);
chrom get_best(pchrom pop, struct info d, chrom best);
void write_best(chrom x, struct info d);
void init_rand();
int random_l_h(int min, int max);
float rand_01();
int flip();


// funcao
void evaluate(pchrom pop, struct info d, int mat[][2]);

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

int main(int argc, char *argv[]){
	char        nome_fich[100] = "",tmp[100] = "";
	struct info EA_param;
	pchrom      pop = NULL, parents = NULL;
	chrom       best_run, best_ever;
	int         gen_actual, r, runs, i, inv, mat[MAX_OBJ][2];
	float       mbf = 0.0;

    // Le os argumentos de entrada
	if (argc == 3){
		runs = atoi(argv[2]);
		strcpy(tmp, argv[1]);
	} else if (argc == 2) { // Se o numero de execucoes do processo nao for colocado nos argumentos de entrada, define-o com um valor por defeito
            runs = DEFAULT_RUNS;
            strcpy(nome_fich, argv[1]);
    } else { // Se o nome do ficheiro de informacoes nao for colocado nos argumentos de entrada, pede-o novamente
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
    // Se o numero de execucoes do processo for menor ou igual a 0, termina o programa
	if (runs <= 0)
		return 0;
    //Inicializa a geracao dos numeros aleatorios
	init_rand();
    // Preenche a matriz com dados dos objectos (peso e valor) e a estrutura EA_param que foram definidos no ficheiro de input
	EA_param = read_file(nome_fich, mat); // mat é matriz bidimensional com MAX_OBJ, onde para cada tens o peso e o lucro
	// Faz um ciclo com o numero de execucoes definidas
	for (r=0; r<runs; r++){
		printf("Repeticao %d\n",r+1);
        // gracao da populacao inicial
		pop = init_pop(EA_param); // pop = array de chroms (struct individual) one pop.p[x] têm 0 ou 1 aleatoriamente
        // Avalia a populacao inicial
		evaluate(pop, EA_param, mat); // coloca no fitness de cada pop[i].fitness se é zero (inválido) ou 1 (válido)
		gen_actual = 1; // geracao atual
        // Encontra-se a melhor solucao dentro de toda a populacao
		best_run = get_best(pop, EA_param, best_run);
        // Reserva espaco para os pais da populacao seguinte
		parents = malloc(sizeof(chrom)*EA_param.popsize);
		if (parents==NULL){ printf("Erro na alocacao de memoria\n"); exit(1); }
		// Ciclo de optimizacao
		while (gen_actual <= EA_param.numGenerations){
            // Torneio binario para encontrar os progenitores (ficam armazenados no vector parents)
			tournament(pop, EA_param, parents);
            // Aplica os operadores geneticos aos pais (os descendentes ficam armazenados na estrutura pop)
			genetic_operators(parents, EA_param, pop);
            // Avalia a nova populacao (a dos filhos)
			evaluate(pop, EA_param, mat);
            // Actualiza a melhor solucao encontrada
			best_run = get_best(pop, EA_param, best_run);
			gen_actual++;
		}
		// Contagem das solucoes invalidas
		for (inv=0, i=0; i<EA_param.popsize; i++)
			if (pop[i].valido == 0)
				inv++;
		// Escreve resultados da repeticao que terminou
		printf("\nRepeticao %d:", r);
		write_best(best_run, EA_param);
		printf("\nPercentagem Invalidos: %f\n", 100*(float)inv/EA_param.popsize);
		mbf += best_run.fitness;
		if (r==0 || best_run.fitness > best_ever.fitness)
			best_ever = best_run;
        // Liberta a memoria usada
		free(parents);
		free(pop);
	}
	// Escreve resultados globais
	printf("\n\nMBF: %f\n", mbf/r);  // Média das fitnesses das best runs
	printf("\nMelhor solucao encontrada");
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
void init_rand(){
	srand((unsigned)time(NULL));
}

// Leitura dos par�metros e dos dados do problema
// Par�metros de entrada: Nome do ficheiro e matriz a preencher com os dados dos objectos (peso e valor)
// Par�metros de sa�da: Devolve a estrutura com os par�metros
struct info read_file(char *filename, int mat[][2]){
	struct  info x;
	FILE    *f;
	int     i;

	f = fopen(filename, "rt");
	if (!f){ printf("File not found\n"); exit(1); }
	
	// Leitura dos parametros do problema
	
	char ded;
	fscanf(f," %c", &ded);
	while (ded == 'c') {
		fscanf(f," %*[^\n]");
		fscanf(f," %c", &ded);
	}
	if (ded == 'p'){  //  <> <obj/arestas>
		fscanf(f," edge %d %d",&x.numGenes , &x.capacity );
	}

	x.popsize = 100; //fscanf(f, " pop: %d", &x.popsize);
	x.pm = 0.01; //fscanf(f, " pm: %f", &x.pm);
	x.pr = 0.7; //fscanf(f, " pr: %f", &x.pr);
	x.tsize = 2; //fscanf(f, " tsize: %d", &x.tsize);
	x.numGenerations = 2500;//fscanf(f, " max_gen: %d", &x.numGenerations);
	//x.numGenes = 100;//fscanf(f, " obj: %d", &x.numGenes);
	//x.capacity = 250;//fscanf(f, " cap: %d", &x.capacity);
	if (x.numGenes > MAX_OBJ){
		printf("Number of itens is superior to MAX_OBJ\n");
		exit(1);
	}
	x.ro = 0.0;
	// Leitura dos dados do KSP (peso e lucro)
	for (i=0; i<x.numGenes; i++){
		fscanf(f, " e %d %d", &mat[i][0], &mat[i][1]);
		printf("======================x.numGenes = %d %d\n",mat[i][0], mat[i][1]);
	}
	fclose(f);
	// Devolve a estrutura com os par�metros
	return x;
}

// Simula o lan�amento de uma moeda, retornando o valor 0 ou 1
int flip(){
	if ((((float)rand()) / (float) RAND_MAX) < 0.5)
		return 0;
	else
		return 1;
}

// Criacao da populacao inicial. O vector e alocado dinamicamente
// Par�metro de entrada: Estrutura com par�metros do problema
// Par�metro de sa�da: Preenche da estrutura da popula��o apenas o vector bin�rio com os elementos que est�o dentro ou fora da mochila
pchrom init_pop(struct info d){
	int     i, j;
	pchrom  indiv;

	indiv = malloc(sizeof(chrom)*d.popsize);
	if (indiv==NULL){ printf("Erro na alocacao de memoria\n"); exit(1);	}
	for (i=0; i<d.popsize; i++){
		for (j=0; j<d.numGenes; j++)
			indiv[i].p[j] = flip(); // 0 ou 1 aleatoriamente
	}
	return indiv;
}

// Actualiza a melhor solu��o encontrada
// Par�metro de entrada: populacao actual (pop), estrutura com par�metros (d) e a melhor solucao encontrada at� a gera��oo imediatamente anterior (best)
// Par�metro de sa�da: a melhor solucao encontrada at� a gera��o actual
chrom get_best(pchrom pop, struct info d, chrom best){
	int i;

	best = pop[0]; // Como ainda nao existe, escolhe-se como melhor solucao a primeira da populacao (poderia ser outra qualquer)
	for (i=0; i<d.popsize; i++){
		if (best.fitness < pop[i].fitness)
			best=pop[i];
	}
	return best;
}

// Devolve um valor inteiro distribuido uniformemente entre min e max
int random_l_h(int min, int max){
	return min + rand() % (max-min+1);
}

// Devolve um valor real distribuido uniformemente entre 0 e 1
float rand_01(){
	return ((float)rand())/(float)RAND_MAX;
}

// Escreve uma solu��o na consola
// Par�metro de entrada: populacao actual (pop) e estrutura com par�metros (d)
void write_best(chrom x, struct info d){
	int i;

	printf("\nBest individual: %4.1f\n", x.fitness);
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
	for(i=0; i<d.popsize;i++){
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
void genetic_operators(pchrom parents, struct info d, pchrom offspring){
    // Recombina��o com um ponto de corte
	crossover(parents, d, offspring);
	// Mutacao binaria // aleatoriamente muda bits para sair de máximos locais
	mutation(offspring, d);
}

// Preenche o vector descendentes com o resultado das opera��es de recombina��o
// Par�metros de entrada: estrutura com os pais (parents), estrutura com par�metros (d), estrutura que guardar� os descendentes (offspring)
void crossover(pchrom parents, struct info d, pchrom offspring){
	int i, j, point;

	for (i=0; i<d.popsize; i+=2){
		if (rand_01() < d.pr){
			point = random_l_h(0, d.numGenes-1);
			for (j=0; j<point; j++){
				offspring[i].p[j] = parents[i].p[j];
				offspring[i+1].p[j] = parents[i+1].p[j];
			}
			for (j=point; j<d.numGenes; j++){
				offspring[i].p[j]= parents[i+1].p[j];
				offspring[i+1].p[j] = parents[i].p[j];
			}
		}
		else{
			offspring[i] = parents[i];
			offspring[i+1] = parents[i+1];
		}
	}
}

// Muta��o bin�ria com v�rios pontos de muta��o
// Par�metros de entrada: estrutura com os descendentes (offspring) e estrutura com par�metros (d)
void mutation(pchrom offspring, struct info d){
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


// Calcula a qualidade de uma solu��o
// Par�metros de entrada: solu��o (sol), capacidade da mochila (d), matriz com dados do problema (mat) e numero de objectos (v)
// Par�metros de sa�da: qualidade da solu��o (se a capacidade for excedida devolve 0)
float eval_individual(int sol[], struct info d, int mat[][2], int *v){
	int     i;
	float   sum_weight, sum_profit;

	sum_weight = sum_profit = 0;
	// Percorre todos os objectos
	for (i=0; i < d.numGenes; i++){
        // Verifica se o objecto i esta na mochila
		if (sol[i] == 1){
            // Actualiza o peso total
			sum_weight += mat[i][0];
            // Actualiza o lucro total
			sum_profit += mat[i][1];
		}
	}
	if (sum_weight > d.capacity){
        // Solucao inv�lida
		*v = 0;
		return 0;
	}
	else{
        // Solucao v�lida
		*v = 1;
		return sum_profit;
	}
}

// Avaliacao da popula��o
// Par�metros de entrada: populacao (pop), estrutura com parametros (d) e matriz com dados do problema (mat)
// Par�metros de sa�da: Preenche pop com os valores de fitness e de validade para cada solu��o
void evaluate(pchrom pop, struct info d, int mat[][2]){
	int i;

	for (i=0; i<d.popsize; i++)
		pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido);
}