#include <stdio.h>
#include <string.h>
#include <stdlib.h>

inline char complement(char dnaBasePair){
    if(dnaBasePair == 'A')
	return 'T';
    else
	if(dnaBasePair == 'C')
	    return 'G';
	else
	    if(dnaBasePair == 'G')
		return 'C';
	    else
		if(dnaBasePair == 'T')
		    return 'A';
		else
		    if(dnaBasePair == 'a')
			return 't';
		    else
			if(dnaBasePair == 'c')
			    return 'g';
			else
			    if(dnaBasePair == 'g')
				return 'c';
			    else
				if(dnaBasePair == 't')
				    return 'a';				
    return dnaBasePair;
}

int counterCharsFirstLine( FILE *fin){
    char line[100];
    fgets(line, sizeof(line), fin);
    return (strlen(line)-1);
}

int main(int argc, char *argv[]){
    int start=-1;
    int end=-1;
    int startIndex;
    int endIndex;

    char fileToOpen[100];
    FILE *fin;
    int numberOfCharsPerLine;

    char defline[100];
    char deflineOfInputLength;

    char buffer[100];
    
    int i=0;
    char characterToPrint;
    int counterOfChars;

    int test;
    int lineNumber;
    int forward;

    int supressDefline=0;

    char * usage="usage:\nseq-substring -s [START] -e [END] -f [FASTA_FILE] -d (to supress the defline)\nwhere [START] and [END] are 1-based\n";


    /* CAPTURE ARGUMENTS */
    if(argc == 1){
	printf("%s",usage);
	return 0;
    }

    for(i=1;i<argc;i++){
	if(strcmp(argv[i],"-h") == 0 ){
	    printf("%s",usage);
	    return 0;
	}
	if(strcmp(argv[i],"-s") == 0 ){
	    start=atoi(argv[i+1]);
	}

	if(strcmp(argv[i],"-e") == 0 ){
	    end=atoi(argv[i+1]);
	}

	if(strcmp(argv[i],"-d") == 0 ){
	    supressDefline=1;
	}

	if(strcmp(argv[i],"-f") == 0 ){
	    fileToOpen[0]='\0';
	    strcpy(fileToOpen,argv[i+1]);
	}
    }
    i=0;   
    /* END CAPTURE ARGUMENTS */
    
    if(start<=0){
	printf("Error: wrong start coordinate\n%s",usage);
	return 0;
    }

    if(start<=end){
	forward=1;
    }else{
	forward=0;
    }

    fin = fopen(fileToOpen,"r" );
    fgets(buffer, sizeof(buffer), fin);

    if(buffer[0] != '>'){
	printf("Wrong fasta format\n\n");
	return 1;
    }
    deflineOfInputLength=(strlen(buffer));
    numberOfCharsPerLine=counterCharsFirstLine(fin);
    



    /* DEFINING DEFLINE */
    if(!supressDefline){
	defline[0]='\0';
	strncat(defline,buffer,strlen(buffer)-1);
	printf("%s SUB{ %d..%d }",defline,start,end);
    }
    /* END DEFINING DEFLINE */



    lineNumber=((start-1)/numberOfCharsPerLine);
    test=deflineOfInputLength+(lineNumber*(numberOfCharsPerLine+1))+((start-1)%numberOfCharsPerLine);

    startIndex = test;


    lineNumber=((end-1)/numberOfCharsPerLine);
    test=deflineOfInputLength+(lineNumber*(numberOfCharsPerLine+1))+((end-1)%numberOfCharsPerLine);

    if(forward){	  
	endIndex   = (test+1);	
    }else{
	endIndex   = (test-1);
    }

    i=startIndex;
    fseek (fin,i,SEEK_SET);

    counterOfChars=0;

    while(i != endIndex){
	characterToPrint=fgetc(fin);

	if(forward){	   
	    i+=1;
	}else{
	    characterToPrint=complement(characterToPrint);
	    i-=1;
	    fseek (fin,i,SEEK_SET);
	}

	if(characterToPrint != '\n'){
	    if( counterOfChars % 60 == 0   ){
		putc('\n',stdout);
	    }
	    putc(characterToPrint,stdout);
	    counterOfChars++;
	}	
    }
    
    putc('\n',stdout);
    fclose(fin);
    return 0;
}
