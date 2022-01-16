#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/wait.h>

int		algorithm_to_tester[2];
FILE	*fp;

int	main(void)
{
	int		pid;
	int		status;
	char	algorithm_output[50];
	int		bytes;
	int		cycle = -1;
	char	*files[] =
	{"tests/test.txt",
		"tests/file1.txt",
		"tests/file2.txt",
		"tests/file3.txt",
		"tests/file4.txt",
		"tests/file5.txt",
		"tests/file6.txt",
		NULL};
	char	*population[] =
	{"30",
		"100",
		NULL};
	char	*mutation_prob[] =
	{"0",
		"0.01",
		"0.1",
		"0.5",
		"1",
		NULL};
	char	*cross_prob[] =
	{"0",
		"0.01",
		"0.1",
		"0.5",
		"1",
		NULL};
	char	*generations[] =
	{
		"100",
		"1000",
		"10000",
		NULL};
	char	*runs[] =
	{"10",
		NULL};
	char	*repair[] =
	{"0",
		"1",
		NULL};
	char	*cross[] =
	{"1",
		"2",
		NULL};
	char	*mutation[] =
	{"1",
		"2",
		NULL};
	char	*multiplier[] =
	{"0.9",
		"0.99",
		NULL};
	char	*t_max[] =
	{"100.0",
		"1.0",
		NULL};
	char	*t_min[] =
	{"0.1",
		"0.001",
		"0.00001",
		NULL};
	char	*equalcost[] =
	{"1",
		"0",
		NULL};
	char	*neighbour[] =
	{"1",
		"2",
		NULL};
	if (pipe(algorithm_to_tester) == -1)
	{
		fprintf(stderr, "Couldn't create pipe\n");
		exit(EXIT_FAILURE);
	}
	if ((fp = fopen("mutant_fix", "w")) == NULL)
	{
		fprintf(stderr, "Couldn't open genetic algorithm test output file\n");
		close(algorithm_to_tester[0]);
		close(algorithm_to_tester[1]);
		exit(EXIT_FAILURE);
	}
	/*
	for (int i = 0; files[i]; i++)
	{
		for (int j = 0; population[j]; j++)
		{
			for (int k = 0; mutation_prob[k]; k++)
			{
				for (int l = 0; cross_prob[l]; l++)
				{
					for (int m = 0; generations[m]; m++)
					{
						for (int n = 0; runs[n]; n++)
						{
							for (int o = 0; repair[o]; o++)
							{
								for (int p = 0; cross[p]; p++)
								{
									for (int q = 0; mutation[q]; q++)
									{
										pid = fork();
										if (pid == 0)
										{
											close(1);
											if (dup(algorithm_to_tester[1]) == -1)
											{
												fprintf(stderr, "Couldn't dup pipe[1]\n");
												exit(0);
											}
											execl("genetic_algorithm",
													"genetic_algorithm",
													files[i],
													population[j],
													mutation_prob[k],
													cross_prob[l],
													generations[m],
													runs[n],
													repair[o],
													cross[p],
													mutation[q],
													NULL);
										}
										if ((bytes = read(algorithm_to_tester[0], algorithm_output, sizeof(algorithm_output) - 1)) == -1)
										{
											fprintf(stderr, "Couldn't read from pipe!\n");
											exit(EXIT_FAILURE);
										}
										algorithm_output[bytes] = '\0';
										fprintf(fp, "Genetic Algorithm,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
												files[i],
												population[j],
												mutation_prob[k],
												cross_prob[l],
												generations[m],
												runs[n],
												repair[o],
												cross[p],
												mutation[q],
												algorithm_output);
										printf("Cycle: %d\n", ++cycle);
										wait(&status);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	/*
	for (int i = 0; files[i]; i++)
	{
		for (int j = 0; multiplier[j]; j++)
		{
			for (int k = 0; t_max[k]; k++)
			{
				for (int l = 0; t_min[l]; l++)
				{
					for (int m = 0; runs[m]; m++)
					{
						for (int n = 0; repair[n]; n++)
						{
							for (int o = 0; equalcost[o]; o++)
							{
								for (int p = 0; neighbour[p]; p++)
								{
									pid = fork();
									if (pid == 0)
									{
										close(1);
										if (dup(algorithm_to_tester[1]) == -1)
										{
											fprintf(stderr, "Couldn't dup pipe[1]\n");
											exit(0);
										}
										execl("simulated_annealing",
												"simulated_annealing",
												files[i],
												t_min[l],
												t_max[k],
												multiplier[j],
												runs[m],
												repair[n],
												equalcost[o],
												neighbour[p],
												NULL);
										fprintf(stderr, "Error occured while trying to execute simulated_annealing\n");
									}
									if ((bytes = read(algorithm_to_tester[0], algorithm_output, sizeof(algorithm_output) - 1)) == -1)
									{
										fprintf(stderr, "Couldn't read from pipe!\n");
										exit(EXIT_FAILURE);
									}
									algorithm_output[bytes] = '\0';
									wait(&status);
									fprintf(fp, "Simulated Annealing,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
											files[i] + 6,
											multiplier[j],
											t_max[k],
											t_min[l],
											runs[m],
											repair[n],
											equalcost[o],
											neighbour[p],
											algorithm_output);
								}
							}
						}
					}
				}
			}
		}
	}
	*/
	for (int i = 0; files[i]; i++)
	{
		for (int j = 0; multiplier[j]; j++)
		{
			for (int k = 0; t_max[k]; k++)
			{
				for (int l = 0; t_min[l]; l++)
				{
					for (int m = 0; runs[m]; m++)
					{
						for (int n = 0; repair[n]; n++)
						{
							for (int o = 0; equalcost[o]; o++)
							{
								for (int p = 0; neighbour[p]; p++)
								{
									pid = fork();
									if (pid == 0)
									{
										close(1);
										if (dup(algorithm_to_tester[1]) == -1)
										{
											fprintf(stderr, "Couldn't dup pipe[1]\n");
											exit(0);
										}
										execl("mutating_simulated_annealing",
												"mutating_simulated_annealing",
												files[i],
												t_min[l],
												t_max[k],
												multiplier[j],
												runs[m],
												repair[n],
												equalcost[o],
												neighbour[p],
												NULL);
										fprintf(stderr, "Error occured while trying to execute mutating_simulated_annealing\n");
									}
									if ((bytes = read(algorithm_to_tester[0], algorithm_output, sizeof(algorithm_output) - 1)) == -1)
									{
										fprintf(stderr, "Couldn't read from pipe!\n");
										exit(EXIT_FAILURE);
									}
									algorithm_output[bytes] = '\0';
									wait(&status);
									fprintf(fp, "Mutant Simulated Annealing,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
											files[i] + 6,
											multiplier[j],
											t_max[k],
											t_min[l],
											runs[m],
											repair[n],
											equalcost[o],
											neighbour[p],
											algorithm_output);
								}
							}
						}
					}
				}
			}
		}
	}
	/*
	for (int i = 0; files[i]; i++)
	{
		for (int j = 0; multiplier[j]; j++)
		{
			for (int k = 0; t_max[k]; k++)
			{
				for (int l = 0; t_min[l]; l++)
				{
					for (int m = 0; runs[m]; m++)
					{
						for (int n = 0; repair[n]; n++)
						{
							for (int o = 0; equalcost[o]; o++)
							{
								for (int p = 0; neighbour[p]; p++)
								{
									pid = fork();
									if (pid == 0)
									{
										close(1);
										if (dup(algorithm_to_tester[1]) == -1)
										{
											fprintf(stderr, "Couldn't dup pipe[1]\n");
											exit(0);
										}
										execl("smarter_simulated_annealing",
												"smarter_simulated_annealing",
												files[i],
												t_min[l],
												t_max[k],
												multiplier[j],
												runs[m],
												repair[n],
												equalcost[o],
												neighbour[p],
												NULL);
										fprintf(stderr, "Error occured while trying to execute smarter_simulated_annealing\n");
									}
									if ((bytes = read(algorithm_to_tester[0], algorithm_output, sizeof(algorithm_output) - 1)) == -1)
									{
										fprintf(stderr, "Couldn't read from pipe!\n");
										exit(EXIT_FAILURE);
									}
									algorithm_output[bytes] = '\0';
									wait(&status);
									fprintf(fp, "Smarter Simulated Annealing,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
											files[i] + 6,
											multiplier[j],
											t_max[k],
											t_min[l],
											runs[m],
											repair[n],
											equalcost[o],
											neighbour[p],
											algorithm_output);
								}
							}
						}
					}
				}
			}
		}
	}
	*/
	fclose(fp);
	close(algorithm_output[0]);
	close(algorithm_output[1]);
	return (0);
}
