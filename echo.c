#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdio.h>


#define MATRIX_SIZE 100
#define LINE_SIZE 200
#define INITIATOR 0

#define SONDA_MESSAGE 10
#define ECHO_MESSAGE 20
#define CONTROL_MESSAGE 30
#define DATA_MESSAGE 40

//tema 3
#define EFFECT_MESSAGE 50
#define TOP_MESSAGE 60
#define BOTTOM_MESSAGE 70
#define SIZE_MESSAGE 80
#define CHUNK_MESSAGE 90

//tema 3
#define BLUR 1010
#define SMOOTH 2020
#define SHARPEN 3030
#define MEAN_REMOVAL 4040

#define TRUE 1 
#define FALSE 0

#define SEND 1000
#define RECEIVE 1001

	int ** matrix; 	
	int * sudokuMatrix;	
	int size;
	int rank;
	int i , j;
	int read;
	int topoSize;
	int sqrtTopoSize;
	int parent;
	int * top_nou;
	char * line ;
	size_t len;
	
	int numarSolutii;
	int numarPrimite;
	int numarDeTrimis;
	int * solutii ;
	int * primite ;
	int * aux;
	int * deTrimis;

	int const smooth_matrix[3][3] = {{1,1,1},{1,1,1},{1,1,1}};
	int const blur_matrix[3][3] = {{1,2,1},{2,4,2},{1,2,1}};
	int const sharpened_matrix[3][3] = {{0,-2,0},{-2,11,-2},{0,-2,0}};
	int const mean_removal_matrix[3][3] = {{-1,-1.-1},{-1,9,-1},{-1,-1,-1}};

	int const smooth_factor = 9;
	int const blur_factor = 16;
	int const sharpen_factor = 3;
	int const mean_removal_factor = 1;

void initTopology(int ** topology , int size ,int value);
int * parseInputAsArray(char * topologyName, int size , char * mode , int rank);
int * getSudokuFragment(char * filename , int rank);
int ** createTopologyUsingMessages(int size , int rank , int * parent , int * adiacenta , int topology[size][size] , int emptyMatrix[size][size]);
int isEmptyMessage(int * receivedMessage, int size );
int isEmptyMatrix(int size, int matrix[size][size]);
int getNumberOfNodes(char * filename , char * mode );
int getNumberOfNeighbors(int size , int rank , int parent , int topology[size][size]);
int * computeLocalSolutii(int size ,int rank , int matrix[size][size]);
void primesteSudoku (int * from  , int * to , int size, int startPoint);
void combineMatrixAdiacenta(int size ,int rank , int matrix[size][size] , int adiacenta[size]);
void logicalORMatrix(int size, int from[size][size], int to[size][size]);
void sendMatrixToAll(int size , int matrix[size][size]);
void createRoutingVector(int size , int rank , int parent , int matrix[size][size], int * vector);
void addResult(int size, int solutionCount , int * result , int * solutions);
int sudoku (int size , int line , int col , int * matrix , int * solutii);
int isValid (int size , int line , int col , int value , int * matrix);
void receiveSolution(int size , int rank , int parent , int * primite , int topology[size][size]);
int * extractSolution(int size , int rank , int * solutii);
void transportSolutieToAux(int size , int rank , int parent , int * solutii , int * aux);
void receiveMessagesFromChildren(int topologySize, int matrixSize , int rank , int parent , int neighborCount , int * primite , int topology[size][size]);
void sendMessagesWithSolutions(int topoSize , int matrixSize , int numarSolutii , int rank , int parent , int * aux , int topology[size][size]);
void generateValidSolutions (int topoSize , int sqrtSize , int rank , int * primite , int * aux , int * deTrimis);
int validateSolution(int size , int sqrtSize , int * matrix);
int * combineMatrixToMatrix (int size , int * from , int * to);
void printMatrix(int size , int matrix[size][size]);
void printArray(int size , int array[size]);
void printVectorMatrix(int size , int * matrix);
void printMessage(int source , int destination , int * array , int size , int messageTYPE , int direction);
void printMessageMatrix(int a , int b , int size , int messageType, int direction, int matrix[size][size]);
void parseImages(char *filename, int size, int topology[size][size]);
void startProcessing(int parent, int rank, int size, int topoloy[size][size]);
void send_chunks(int size,int topology[size][size],int x,int y,unsigned char *pixels,int rank,int parent);
void put_pixels(int chunk_size,unsigned char *pixels,int x,int y, unsigned char **bordered_matrix);

int main(int argc , char ** argv){
	
	int value; //DEBUG
	int initialized = FALSE ;
	numarSolutii = 0;
	
	MPI_Init(& argc ,& argv);
	MPI_Comm_size(MPI_COMM_WORLD , & size);
	MPI_Comm_rank(MPI_COMM_WORLD , & rank);
	
	//number of nodes
	//not square
	//sqrtTopoSize = getNumberOfNodes(argv[2], "r+");
	//solutii = (int *) calloc (10000 * sqrtTopoSize * sqrtTopoSize , sizeof(int));
	//primite = (int *) calloc (10000 * sqrtTopoSize * sqrtTopoSize, sizeof(int));
	
	
	// topoSize = sqrtTopoSize * sqrtTopoSize;
	topoSize = size;
	//square
	int emptyMatrix[topoSize][topoSize];
	int topology[topoSize][topoSize];
	int routingVector[topoSize];

	aux = (int * )calloc (10000 * topoSize * topoSize , sizeof(int));
	deTrimis = (int *) calloc (10000 * topoSize * topoSize , sizeof(int));
		
	for(i = 0 ; i < size ; ++i){
		routingVector[i] = -1;
	}
	
	for(i = 0 ; i < size ; ++i){
		for(j = 0 ; j < size ; ++j){
			emptyMatrix[i][j] = 0;
			topology[i][j] = 0;
		}
	}
	
	//arrays of neighbors corresponding to the line in the topology
	//matrix
	top_nou = parseInputAsArray(argv[1], size, "r+", rank);
	
	//determine the topology using message passing with wave algorithm
	matrix = createTopologyUsingMessages(topoSize , rank , &parent , top_nou , topology , emptyMatrix);
	
	//all nodes have the topology matrix after this line
	MPI_Bcast(&topology , topoSize * topoSize , MPI_INT , 0 , MPI_COMM_WORLD);
	
	//parseImages
	unsigned char *image;
	if (rank == 0){
		parseImages(argv[2], size, topology);
	}else{
		//receive image from parent
		//send to neighbors
		startProcessing(parent, rank, size, topology);
	}

	printf("Rank %d has\n", rank);
	printMatrix(topoSize, topology);

	/*
	 * Apply filters and count statistics
	 */


	/*
	createRoutingVector(topoSize , rank , parent , topology, &routingVector[0]);
	sudokuMatrix = getSudokuFragment(argv[2] , rank);
	sudoku (sqrtTopoSize , 0 , 0 , sudokuMatrix , solutii);
	int numarVecini = getNumberOfNeighbors(topoSize , rank , parent , topology);
	// }
	MPI_Barrier(MPI_COMM_WORLD);
	// printf("SQRT = %d , solutii = %d\n" , sqrtTopoSize , numarSolutii);
	transportSolutieToAux(sqrtTopoSize , rank , parent , solutii , aux);
	if(numarVecini != 0){
		receiveMessagesFromChildren(topoSize, topoSize , rank , parent , numarVecini , primite , topology);
		generateValidSolutions (topoSize , sqrtTopoSize , rank , primite , aux , deTrimis);
		if(rank != 0){
			sendMessagesWithSolutions(topoSize , topoSize , numarDeTrimis , rank , parent , deTrimis , topology);
		}
	}
	else{
		if(rank != 0)
			sendMessagesWithSolutions(topoSize , topoSize  , numarSolutii , rank , parent , aux , topology);
	}*/
	
	int k = 0 ;
	MPI_Finalize();
	return 0;
}

void printMatrix(int size , int matrix[size][size]){
	
	int i = 0 ;
	int j = 0 ;
	for (i = 0 ; i < size ; ++i){
		for(j = 0 ; j < size ; ++j){
			printf("%d " , matrix [i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void initTopology(int ** topology , int size , int value){
	
	int i  =0 ;
	int j  =0 ;
	
	for(i = 0 ; i < size ; ++i){
		topology[i] = (int *) calloc (size , sizeof(int));
	}
	
	for(i = 0; i < size ; ++i){
		for(j = 0 ; j < size ; ++j){
			topology[i][j] = value;
		}
	}

}

int ** createTopologyUsingMessages(int size , int rank , int * parent , int * adiacenta , int topology[size][size] , int emptyMatrix[size][size]){
		
	int top_nou[size][size];
	for(i = 0 ; i < size ; ++i){
		for(j = 0 ; j < size ; ++j){
			top_nou[i][j] = 0;
		}
	}
	
	int ** topology_pointer = (int **) calloc (size , sizeof(int *));
	for(i = 0 ; i < size ; ++i){
		topology_pointer[i] = (int *) calloc (size , sizeof(int ));
	}
	combineMatrixAdiacenta(size ,rank , topology , adiacenta );

	MPI_Status status;	
	
	if(rank == 0){ //INITIATOR
		for(i = 0 ; i < size ; ++i){
			if(adiacenta[i] == 1){
				MPI_Send(emptyMatrix , size * size , MPI_INT , i , SONDA_MESSAGE , MPI_COMM_WORLD);
			}
		}
	}
	
	else{
		
		//receive sonda
		MPI_Recv(top_nou , size * size , MPI_INT , MPI_ANY_SOURCE , SONDA_MESSAGE , MPI_COMM_WORLD , &status);
		*parent = status.MPI_SOURCE;
			
		//send sonde to neighbors
		for(i = 0 ; i < size ; ++i){
			if(adiacenta[i] == 1 && i != *parent){
				MPI_Send(emptyMatrix , size * size , MPI_INT , i , SONDA_MESSAGE , MPI_COMM_WORLD);
			}
		}
		
	}
		//get the number of echos I need to receive
		int numberOfEcho = 0;
		for(i = 0 ; i < size ; ++i){
			if(adiacenta[i] == 1 && i != *parent){
				numberOfEcho ++;
			}
		}
		
		while(numberOfEcho > 0){
			MPI_Recv(top_nou , size * size , MPI_INT , MPI_ANY_SOURCE , MPI_ANY_TAG , MPI_COMM_WORLD , &status);
			int source = status.MPI_SOURCE;
			
			//received echo
			if(status.MPI_TAG == ECHO_MESSAGE){
				if(top_nou != NULL){
					int empty = isEmptyMatrix(size , top_nou);
						if(!empty){
							logicalORMatrix(size , top_nou ,topology);
						}
					numberOfEcho --;
				}
			}
			//received sonda when listening for echo
			else if(status.MPI_TAG == SONDA_MESSAGE){
				MPI_Send(emptyMatrix , size * size , MPI_INT , source , ECHO_MESSAGE , MPI_COMM_WORLD);
				//delete connection
				if(topology[rank][source] == 1){
					topology[rank][source] =0;
				}
			}
			
		}

		//send echo to parent
		if(rank != 0){
			MPI_Send(topology , size * size , MPI_INT , *parent , ECHO_MESSAGE , MPI_COMM_WORLD);
		}else{
					
			for(i = 0 ; i < size ; ++i){
				for(j = 0 ; j < size ; ++j){
					topology_pointer[i][j] = topology[i][j];
				}
			}
			
			return topology_pointer;
		}
	
}

int isEmptyMessage(int * receivedMessage , int size){
	int i;

	for(i = 0 ; i < size ; ++i){
		if(receivedMessage[i] != 0)
			return FALSE;
			
		return TRUE;
	}
}

void printArray(int size , int array[size]){
	int i = 0;
	for(i = 0 ; i < size ; ++i){
		printf("%d ", array[i]);
	}
}


int * parseInputAsArray(char * topologyName, int size, char * mode , int grad){
	int i = 0 ;
	int currentPosition =0 ;
	int * outArray; 
	 
	outArray = (int *) calloc (size * size , sizeof(int));
	FILE * topologyFile = fopen(topologyName , mode );

	/* Read from file */
	while ((read = getline(&line, &len, topologyFile)) != -1) {
						
		int parinte;
		int copil;
		char * tok = strtok(line , ":- ");
		
		sscanf(tok, "%d" , &parinte);
		//get adjancency nodes
		if(parinte == grad){
			while(tok != NULL){
				tok = strtok(NULL , ":- ");
				if(tok != NULL){
					sscanf(tok, "%d" , &copil);
					outArray[copil] = 1;
				}
			}
		}				
	}	
	fclose(topologyFile);
	return outArray;
}


int getNumberOfNodes(char * filename , char * mode ){
	
	int size ;
	
	FILE * sudokuFile = fopen(filename, mode);
	fscanf(sudokuFile , "%d" , &size);
	fclose(sudokuFile);
	
	return size;
	
}

void primesteSudoku (int * from  , int * to , int size, int startPoint){
	int i = 0; 
	
	for(i = startPoint ; i < startPoint + size ; ++i){
		to[i] = from [i % size];
	}
	
	numarPrimite ++;
}

void copy(int * from , int * to , int size){
	int i ;
	for (i = 0; i < size ; ++i){
		to[i] = from[i];
	}
}

void logicalOR(int * from , int * to , int size){
	int i ;
	
	for (i = 0 ; i < size ; ++i){
		to[i] |= from[i];
	}
}

void printMessage(int a , int b , int * array , int size , int messageTYPE , int direction){
	int i ;
	if(direction == SEND)
		printf("%d trimite ",a);
	else 
		printf("%d primeste ",a);
	printArray(size, array);
	if(direction == SEND)
		printf(" catre %d ", b);
	else
		printf(" de la %d ", b);
		
	if(messageTYPE == ECHO_MESSAGE)
		printf("echo");
	else
		printf("sonda");
	printf("\n");		
}

int isEmptyMatrix(int size, int matrix[size][size]){
	int i ;
	int j ;
	for (i = 0; i < size ; ++i ){
		for (j = 0 ; j < size ; ++j){
			if(matrix[i][j] != 0)
				return FALSE;
		}
	}
	
	return TRUE;
}

void printMessageMatrix(int a , int b , int size , int messageType, int direction, int matrix[size][size]){
	int i ;
	if(direction == SEND)
		printf("%d trimite ",a);
	else 
		printf("%d primeste ",a);

	if(direction == SEND)
		printf("catre %d ", b);
	else
		printf("de la %d ", b);
		
	
	if (messageType == ECHO_MESSAGE && isEmptyMatrix(size , matrix))
		printf("echo empty");
	else if(messageType == ECHO_MESSAGE)
		printf("echo");
	else
		printf("sonda");
	printf("\n");	
}

void logicalORMatrix(int size, int from[size][size], int to[size][size]){
	int i ;
	int j ;
		
	int ** aux = (int **) calloc (size , sizeof(int *));
	for(i = 0 ; i < size ; ++i){
		aux[i] = (int *) calloc (size , sizeof(int));
	}
	
	for(i = 0 ; i < size ; ++i){
		for (j = 0 ; j < size ; ++j){
			aux[i][j] = to[i][j];
		}
	}
	
	for (i = 0 ; i < size ; ++i){
		for (j = 0 ; j < size ; ++j){
			if(from[i][j] == 1){
				aux[i][j] = 1;
			}
		}
	}
		
	for(i = 0 ; i < size ; ++i){
		for (j = 0 ; j < size ; ++j){
			to[i][j] = aux[i][j];
		}
	}
}

void combineMatrixAdiacenta(int size ,int rank , int matrix[size][size] , int adiacenta[size]){
	
	int i ;
	
	for (i = 0; i < size ; ++i){
		matrix[rank][i] |= adiacenta[i];
	}
}

void createRoutingVector(int size , int rank , int parent , int matrix[size][size], int * vector){
	
	int i = 0 ;
	int j = 0;
		
	for(i = 0 ; i < size ; ++i)
		vector[i] = -1;
		
	for(i = 0 ; i < size ; ++i){
		if(matrix[rank][i] == 1){
			vector[i] = i;
		}
	}
	
	for(i = 0 ; i < size ; ++i){
		if(vector[i] == -1 && i != rank){
			for (j = rank ; j < size ; ++j){
				if(matrix[j][i] == 1 && matrix[rank][j] == 1){
					vector[i] = j;
				}
			}
		}
	}
	
	for(i = 0; i < size ; ++i){
		if(vector[i] == -1 && i != rank){
			vector[i] = parent;
		}
	}
}

int * getSudokuFragment(char * filename , int rank){
	
	int i ;
	int j ;
	int size;
	int squareSize ;
	
	FILE * inFile = fopen (filename , "r+");
	
	fscanf(inFile , "%d\n" , &size);
	squareSize = size * size ;
	
	int * sudokuMap = (int *) calloc (squareSize * squareSize , sizeof(int));
	int * sudokuPart = (int *) calloc (squareSize , sizeof(int));
	
	for(i = 0 ; i < squareSize ; ++i){
		for (j = 0 ; j < squareSize ; ++j){
			fscanf(inFile , "%d", &sudokuMap[i * squareSize + j]);
		}
	}
	
	int linStart = (rank / size) * size ;
	int colStart =  (rank % size ) * size ;
	
	for(i = linStart ; i < linStart + size ; ++i){
		for (j = colStart ; j < colStart + size ; ++j){
			int localLine = i % size;
			int localCol = j % size;
			sudokuPart[localLine * size + localCol] = sudokuMap[i * squareSize + j];
		}
	}
	
	fclose(inFile);
	return sudokuPart;
	
}

int sudoku (int size , int line , int col , int * matrix , int * solutii){

	int i = 0; 
	int j = 0;
	
	if(line >= size){
		addResult(size , numarSolutii , matrix , solutii);
		numarSolutii++;
		return TRUE;
	}
	
	if(col == size)
		return sudoku(size , line + 1 , 0 , matrix , solutii);
	
	if(matrix[line * size + col] != 0){
		return sudoku(size , line , col + 1 , matrix , solutii);
	}
	
	for(i = 1 ; i < size * size + 1 ; ++i){
		if(isValid (size , line , col , i , matrix)){
			
			matrix[line * size + col] = i;
						
			int  result; 
			if(col == size)
				result = sudoku(size , line + 1 , 0 , matrix , solutii);
			result = sudoku (size , line , col + 1 , matrix , solutii);
			
			if(!result)
				return FALSE;
				
			matrix[line * size + col] = 0 ;
		}
	}
}

int isValid (int size , int line , int col , int value , int * matrix){
	
	int i;
	int j;
		
	for(i = 0 ; i < size ; ++i){
		if(matrix[line * size + i] == value){
			return FALSE;
		}
	}
	
	for(j = 0 ; j < size ; ++j){
		if(matrix[j * size + col] == value){
			return FALSE;
		}
	}
	
	for(i = 0 ; i < size ; ++i){
		for(j = 0 ; j < size ; j++){
			if(matrix[i * size + j] == value){
				return FALSE;
			}
		}
	}
		
	return TRUE;
}


int * computeLocalSolutii(int size ,int rank , int matrix[size][size]){
	
	int * solutions = (int *) calloc (size * size * size , sizeof(int));
	int current = 0 ;
	
	return NULL;
}

int getNumberOfNeighbors(int size , int rank , int parent , int topology[size][size]){
	
	int i = 0;
	int summ = 0; 

	if (rank == 0){
		for (i = 0; i < size; ++i){
			if(topology[rank][i] == 1){
				summ++;
			}
		}
		return summ;
	}
	
	for(i = 0 ; i < size ; ++i){
		if(topology[rank][i] == 1 && i != rank  && i != parent){
			summ ++;
		}
	}
	
	return summ;
}

void addResult(int size, int solutionCount , int * result , int * solutions){
	
	int start = solutionCount * size * size ;
	int resultSize = size * size; 
	
	//copy the matrix as array
	for(i = start ; i < start + resultSize ; ++i){
		solutions[i] = result[i % resultSize] ;
	}
	
}

void receiveSolution(int size , int rank , int parent , int * primite , int topology[size][size]){
	int i = 0;
	int j = 0;
	
	int * sol_copil = (int *) calloc (size * size , sizeof(int));
	
	if(getNumberOfNeighbors(size , rank , parent ,topology) > 0){
		int messagesCount = 0;
		
		for(i = 0 ; i < size ; ++i){
			if(topology[rank][i] == 1 && i != parent && i != rank){
				MPI_Recv(&messagesCount , 1 , MPI_INT , i , CONTROL_MESSAGE , MPI_COMM_WORLD , NULL);
			}
			
			for(j = 0 ; j < messagesCount ; ++j){
				MPI_Recv(sol_copil , size * size , MPI_INT , i , DATA_MESSAGE , MPI_COMM_WORLD , NULL);
				primesteSudoku (sol_copil , primite , size * size , numarPrimite);
				messagesCount --;
			}
		}
	}
	
	
}

int * extractSolution(int size , int rank , int * solutii){
	
	int i ;
	int j ;
	
	int * ultimaSolutie = (int *) calloc (size * size , sizeof(int));
	
	int start = (numarSolutii - 1) * size * size;
	for(i = 0 ; i < size * size ; ++ i){
		ultimaSolutie[i] = solutii[start + i] ;
	}	
	
	numarSolutii --;
	
	return ultimaSolutie;
}

//pun fiecare solutie in casuta respectiva rank-ului din matricea mare
void transportSolutieToAux(int size , int rank , int parent , int * solutii , int * aux){
	
	int i = 0;
	int j = 0;
	int k = 0;
	
	int globalSize = size * size ;
	
	int linStart = (rank / size) * size ;
	int colStart =  (rank % size ) * size ;

	for(i = 0 ; i < numarSolutii ; ++i){
		for(j = 0 ; j < size ; ++j){
			for(k = 0 ; k < size ; ++k){
				aux[i * globalSize * globalSize + (linStart + j) * globalSize + (colStart + k)] = solutii[i * size * size + j * size + k];
			}
		}
	}
}

void receiveMessagesFromChildren(int topologySize, int matrixSize , int rank , int parent , int neighborCount , int * primite , int topology[size][size]){
	
	int i= 0;
	int j= 0;
	int k= 0;
	int l= 0;
	int numarMatrici;
	int * matrix = (int *) calloc (topologySize * topologySize , sizeof(int));
	
	for(i = 0 ; i < matrixSize ; ++i){
		if(i != parent && topology[rank][i] == 1){
			MPI_Recv(&numarMatrici , 1 , MPI_INT , i , CONTROL_MESSAGE , MPI_COMM_WORLD , NULL);
			
			for(j = 0 ; j < numarMatrici ; ++j){
				MPI_Recv(matrix , topologySize * topologySize , MPI_INT , i , DATA_MESSAGE , MPI_COMM_WORLD , NULL);
				for(k = 0 ; k < topologySize * topologySize ; ++k){
					primite[numarPrimite * topologySize * topologySize + k] = matrix[k];
					numarPrimite ++;
				}
				
			}
		}
	}
}

void sendMessagesWithSolutions(int topoSize , int matrixSize , int numarSolutii , int rank , int parent , int * aux , int topology[size][size]){
	int i = 0;
	int j = 0;
	
	MPI_Send(&numarSolutii , 1 , MPI_INT , parent , CONTROL_MESSAGE , MPI_COMM_WORLD);
	
	for(i = 0 ; i < numarSolutii ; ++i){
		MPI_Send(&aux[i * topoSize * topoSize] , topoSize * topoSize , MPI_INT , parent , DATA_MESSAGE , MPI_COMM_WORLD);
	}	
}

void generateValidSolutions (int topoSize , int sqrtSize , int rank , int * primite , int * aux , int * deTrimis){
	
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	int value = 0;
	
	for(i = 0 ; i < numarSolutii ; ++i){
		for(j = 0 ; j < numarPrimite ; ++j){
			int * result = combineMatrixToMatrix(topoSize , &aux[i * topoSize * topoSize] , &primite[j * topoSize * topoSize]);
			
			if(validateSolution(topoSize , sqrtSize , result)){
				//copyMatrix
				for(k = 0 ; k < topoSize ; ++k){
					for(l = 0 ; l < topoSize ; ++l){
						deTrimis[numarDeTrimis * topoSize * topoSize + k * topoSize + l] = result[k *topoSize + l];
					}
				}
			numarDeTrimis ++;
			}				
		} 
	}
} 

int validateSolution(int size , int sqrtSize , int * matrix){
	
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	int value ;
		
	for (i = 0 ; i < topoSize ; ++i){
		for(j = 0 ; j < topoSize ; ++j){
			value = matrix[i * topoSize + j];		
			for(k = 0 ; k < topoSize ; ++k){
				for(l = 0 ; l < topoSize ; ++l){
					if(k == i && l == j)
						continue;
					if(matrix[k * topoSize + l] == value && value != 0)
						return FALSE;
				}
			}
		}
	}
	
	return TRUE;
}

int * combineMatrixToMatrix (int size , int * from , int * to){
	
	int i = 0;
	int * result = (int *) calloc (size * size , sizeof(int));
	
	for(i = 0; i < size ; ++i){
		for(j = 0 ; j < size ;++j){
			result[i * size + j] |= from[i * size + j] | to[i * size + j];
		}
	}
	
	return result;
}

void printVectorMatrix(int size , int * matrix){
	int i = 0;
	int j = 0;
	
	for(i = 0 ; i < size ; ++i){
		for(j = 0 ; j < size ; ++j){
			printf("%d " , matrix[i * size + j]);
		}
		
		printf("\n");
	}
	printf("\n");
}

void parseImages(char *filename, int size, int topology[size][size]){

	FILE * images_file = fopen(filename, "r+");
	char *operation;
	char *in_file;
	char *out_file;
	int num_images = 0;
	
	//read number of images
	read = getline(&line, &len, images_file);
	printf("Image file lien = %s\n",line);
	num_images = atoi(line);
	int k;

	for (k = 0; k < num_images; ++k){
		printf("num_images = %d \n", num_images);

		read = getline(&line, &len, images_file);
		operation = strtok(line, " \n\t");
		in_file = strtok(NULL, " \n\t");
		out_file = strtok(NULL, " \n\t");
		printf("Out file = %s \n", out_file);

		//just read the input file 
		FILE *image = fopen(in_file, "r");
		//erase the contents of the file and create the file if it does not exist
		FILE *out_image = fopen(out_file, "w+");
		//copy headers
		char *header_line;
		size_t header_len = 0;
		int header_read = 0;
		int header_i = 0;

		//process_headers
		int max_color = 0;
		int x = 0;
		int y = 0;
		for (header_i = 0; header_i < 4; ++header_i){
			header_read = getline(&header_line, &header_len, image);
			fputs(header_line, out_image);
			//size header
			if (header_i == 2){
				x = atoi(strtok(header_line," "));
				y = atoi(strtok(NULL, " "));
			}else if (header_i == 3){
				max_color = atoi(header_line);
			}
		}

		unsigned char *pixels = (unsigned char *)calloc(x * y, sizeof(unsigned char));
		int pixel_index;
		for (pixel_index = 0; pixel_index < x * y; pixel_index++){
			fscanf(image, "%hhu\n", &pixels[pixel_index]);
		}

		for (pixel_index = 0; pixel_index < x * y; pixel_index++){
			fprintf(out_image, "%hhu\n", pixels[pixel_index]);
		}

		int num_neighbors = getNumberOfNeighbors(size, 0, -1, topology);
		int chunk_size = (x * y)/num_neighbors;
		int current_neighbor = 0;

		for (i = 0; i < size; ++i){
			if (topology[0][i] == 1){

				//send effect matrix
				if (strcmp(operation, "blur") == 0){
					printf("Found blur \n");
					MPI_Send(&blur_factor, 1, MPI_INT, i, EFFECT_MESSAGE, MPI_COMM_WORLD);
				}else if (strcmp(operation, "sharpend") == 0){
					printf("Found sharpend \n");
					MPI_Send(&sharpen_factor, 1, MPI_INT, i, EFFECT_MESSAGE, MPI_COMM_WORLD);
				}else if (strcmp(operation, "smooth") == 0){
					printf("Found smooth \n");
					MPI_Send(&smooth_factor, 1, MPI_INT, i, EFFECT_MESSAGE, MPI_COMM_WORLD);
				}else if (strcmp(operation, "mean_removal") == 0){
					printf("Found mean_removal \n");
					MPI_Send(&mean_removal_factor, 1, MPI_INT, i, EFFECT_MESSAGE, MPI_COMM_WORLD);
				}
				// send x,y
				MPI_Send(&x, 1, MPI_INT, i, SIZE_MESSAGE, MPI_COMM_WORLD);
				MPI_Send(&y, 1, MPI_INT, i, SIZE_MESSAGE, MPI_COMM_WORLD);
				// //send top
				// //send bottom
				// //create bottom and top
				unsigned char *top;
				unsigned char *bottom;

				if (rank == 0 && current_neighbor == 0){
					//top = full of 0's
					top = (unsigned char *)calloc(x, sizeof(unsigned char));
					bottom = &pixels[(current_neighbor + 1)*chunk_size];
					//TODO Send it
					MPI_Send(top, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
					MPI_Send(bottom, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
					free(top);
				}  

				else if (rank == 0 && current_neighbor == num_neighbors -1){
					//bottom = full of 0's
					top = &pixels[current_neighbor * chunk_size - x];
					bottom = (unsigned char *)calloc(x, sizeof(unsigned char));
					//TODO Send it
					MPI_Send(top, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
					MPI_Send(bottom, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
					free(bottom);
				}

				else{
					top = &pixels[current_neighbor * chunk_size - x];
					bottom = &pixels[(current_neighbor + 1)*chunk_size];
					MPI_Send(top, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
					MPI_Send(bottom, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
				}

				//send chunk
				unsigned char *chunk = (unsigned char *)calloc(chunk_size, sizeof(unsigned char));
				int sent_chunks = 0;
				chunk = &pixels[current_neighbor * chunk_size];
				while(sent_chunks < (chunk_size/x)){
					MPI_Send(&chunk[sent_chunks * x], x, MPI_CHAR, i, CHUNK_MESSAGE, MPI_COMM_WORLD);
					sent_chunks++;
				}
				current_neighbor++;
			}
		}

		fflush(out_image);
		fclose(out_image);
	}
	fclose(images_file);
}

void startProcessing(int parent, int rank, int size, int topology[size][size]){
	MPI_Status status;
	int effect_type;
	MPI_Recv(&effect_type, 1, MPI_INT, parent, EFFECT_MESSAGE, MPI_COMM_WORLD, &status);

	int  i = 0;
	for (i = 0; i < size; ++i){
		if (topology[rank][i] == 1){
			//send the effect to neighbors
			MPI_Send(&effect_type, 1, MPI_INT, i, EFFECT_MESSAGE, MPI_COMM_WORLD);
		}
	}
	//recv x and y
	int x = 0;
	int y = 0;
	MPI_Recv(&x, 1, MPI_INT, parent, SIZE_MESSAGE, MPI_COMM_WORLD, &status);
	MPI_Recv(&y, 1, MPI_INT, parent, SIZE_MESSAGE, MPI_COMM_WORLD, &status);

	for (i = 0; i < size; ++i){
		if (topology[rank][i] == 1){
			//send x and y
			MPI_Send(&x, 1, MPI_INT, i, SIZE_MESSAGE, MPI_COMM_WORLD);
			MPI_Send(&y, 1, MPI_INT, i, SIZE_MESSAGE, MPI_COMM_WORLD);
		}
	}

	printf("Rank %d received x=%d and y=%d\n", rank, x, y);
	// //recv top & bottom
	unsigned char *top = (unsigned char *)calloc(x, sizeof(unsigned char));
	unsigned char *bottom = (unsigned char *)calloc(x, sizeof(unsigned char));

		MPI_Recv(top, x, MPI_CHAR, parent, TOP_MESSAGE, MPI_COMM_WORLD, &status);
		MPI_Recv(bottom, x, MPI_CHAR, parent, TOP_MESSAGE, MPI_COMM_WORLD, &status);

		for (int i = 0; i < size; ++i){
			if (topology[rank][i] ==1){
				MPI_Send(top, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
				MPI_Send(bottom, x, MPI_CHAR, i, TOP_MESSAGE, MPI_COMM_WORLD);
			}
		}

		// int received_chunks = 0;
		// unsigned char *chunk = (unsigned char*)calloc(x * y, sizeof(unsigned char));
		// while(received_chunks < y){
		// 	MPI_Recv(&chunk[received_chunks * x], x, MPI_CHAR, i, CHUNK_MESSAGE, MPI_COMM_WORLD, &status);
		// 	received_chunks++;
		// }
}

void send_chunks(int size,int topology[size][size],int x,int y,unsigned char *pixels,int rank,int parent){
	int current_neighbor = 0;
	/*
	 * Calculate chunks
	 */
	int num_neighbors = getNumberOfNeighbors(size, rank, parent, topology);
	int chunk_size = (x * y)/num_neighbors;
	int chunk_left_over = (x*y)%num_neighbors;
	/*
	 * Border the chunks
	 */
	while (num_neighbors > 0){
		if (topology[rank][current_neighbor] == 1){
			unsigned char *bordered_matrix = NULL;
			int start_index = current_neighbor * chunk_size;
			// put_pixels(chunk_size, pixels[start_index], x, y, &bordered_matrix);
		}
	}
	

}

void put_pixels(int chunk_size,unsigned char *pixels,int x,int y, unsigned char **bordered_matrix){
	/*
	 * Initialization
	 */
	if (*bordered_matrix == NULL){
		(*bordered_matrix) =(unsigned char *)calloc(chunk_size + x + x + y + y,sizeof(unsigned char));
	}
}