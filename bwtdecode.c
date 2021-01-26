#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
// #include <limits.h>

#define CHUNK_SIZE 4096
#define SYMBOL_COUNT 5
#define MAX_TABLE_SIZE 235000
#define BLOCK_SIZE 64
#define MAX_FILE_TABLE_SIZE 46880
#define FILE_BLOCK 320
#define SYS_CHUNK 300

/***************************************************************************
*   Function   : getCharIndex
*   Description: This function looksup the hardcoded index for a given
*                character in the alphabet for the corpus of documents.
*   Parameters : symbol - a character in the alphabet {'\n', 'A', 'C', 
*                'G', 'T'}.
*   Effects    : Looksup the index for provided symbol.
*   Returned   : Index for the provided symbol.
***************************************************************************/
int getCharIndex(char symbol){
    switch(symbol){
    case '\n':
        return 0;
    case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    case 'T':
        return 4;
    default:
        return -1;
    }
}

/***************************************************************************
*   Function   : getIndexChar
*   Description: This function looksup the character for a given
*                hardcoded index.
*   Parameters : index - a valid int that ranges from 0 to 4.
*   Effects    : Looksup the character for provided index.
*   Returned   : Character for the provided index.
***************************************************************************/
char getIndexChar(int index){
    switch(index){
    case 0:
        return '\n';
    case 1:
        return 'A';
    case 2:
        return 'C';
    case 3:
        return 'G';
    case 4:
        return 'T';
    default:
        return '\0';
    }
}

/***************************************************************************
*   Function   : setTables
*   Description: This function creates a long array that has 
*                SYMBOL_COUNT rows and MAX_TABLE_SIZE columns.
*   Parameters : NULL
*   Effects    : Mallocs SYMBOL_COUNT rows and MAX_TABLE_SIZE columns.
*   Returned   : Malloc'd long array.
***************************************************************************/
long **setTables(){
    long **tables  = calloc(SYMBOL_COUNT, sizeof(long *));
    assert(tables != NULL);
    for(int i = 0; i < SYMBOL_COUNT; i++){
        tables[i] = calloc(MAX_TABLE_SIZE, sizeof(long));
        assert(tables[i] != NULL);
    }
    return tables;
}

/***************************************************************************
*   Function   : freeTables
*   Description: This function frees a long array that has 
*                SYMBOL_COUNT rows and MAX_TABLE_SIZE columns.
*   Parameters : tables - Malloc'd long array.
*   Effects    : Frees malloc'd SYMBOL_COUNT rows and 
*                MAX_TABLE_SIZE columns.
*   Returned   : NULL
***************************************************************************/
void freeTables(long **tables){
    for(int j = 0; j < SYMBOL_COUNT; j++){
        free(tables[j]);
    }
    free(tables);
}

/***************************************************************************
*   Function   : chunkProcessing
*   Description: This function updates the count for the characters
*                in the alphabet and the occurrence table while 
*                processing the buffer and searching for the index 
*                of terminating character '\n'  
*   Parameters : file_ptr - The buffer currently read from the BWT file.
*                file_sz - The size of the BWT file.
*                index_insert - Current index in the checkpointed 
*                occurrence table.
*                curr_fl_size - The cumulative number of bit read from
*                the BWT file.
*                read_bits - The number of bits in the buffer.
*                term_index - The index of the terminating character '\n'
*                cntTab - The count table to keep track of the cummulative 
*                count of each character in the alphabet.
*                lTab - The occurrence table to keep track of the count
*                each character in the alphabet at a particular checkpoint.
*   Effects    : Updates the count table with the final count of each 
*                character in the alphabet, updates the final counts for all the 
*                checkpoints in the occurrence table and finds the index of
*                the terminating index.
*   Returned   : NULL
***************************************************************************/
void chunkProcessing(char *file_ptr, long file_sz, long *index_lTab, long *index_bwt, long curr_fl_size, size_t read_bits, long *term_index, long *cntTab, long **lTab, char *bwtTab){ 
    long gap_mod = 0;
    long file_mod = 0;
    for(long d = 0; d < read_bits; d++){
        if(file_ptr[d] == '\n'){
            *term_index = d + curr_fl_size;
        }
        if(curr_fl_size == 0 && d == 0){
            cntTab[getCharIndex(file_ptr[0])] += 1;
            lTab[getCharIndex(file_ptr[0])][0] += 1;
            bwtTab[0] = file_ptr[0];
            ++*(index_lTab);
            ++*(index_bwt);
        }
        else{
            int curr_char_idx = getCharIndex(file_ptr[d]);
            cntTab[curr_char_idx] += 1;
            if((gap_mod = (d + curr_fl_size) % BLOCK_SIZE) == 0 || (d + curr_fl_size) == file_sz - 1){
                for(int z = 0; z < SYMBOL_COUNT; z++){
                    lTab[z][*index_lTab] = cntTab[z];
                }
                ++*(index_lTab);
            }
            if((file_mod = (d + curr_fl_size) % FILE_BLOCK) == 0 || (d + curr_fl_size) == file_sz - 1){
                bwtTab[*index_bwt] = file_ptr[d];
                ++*(index_bwt);
            }
        }
    }
}

/***************************************************************************
*   Function   : characterLookup
*   Description: This function reads and returns the character at 
*                a particular index in the file.  
*   Parameters : file - The opened BWT file.
*                index - The index in the BWT file to read.
*   Effects    : Reads the character at the index specfied in the BWT file.
*   Returned   : The character read from the BWT file at the specifed index.
***************************************************************************/
char characterLookup(int file, char *bwtTab, long index, long file_size){
    long mod = index % FILE_BLOCK;
    long div = index / FILE_BLOCK;
    char result = '\0';
    if(index == file_size - 1){
        result = bwtTab[div + 1];
    }
    else if(mod == 0){
        result = bwtTab[div];
    }
    else{
        lseek(file, index * sizeof(char), SEEK_SET);
        // fread(&result, sizeof(char), 1, file);
        if(read(file, &result, 1)){};
    }

    return result;
}

/***************************************************************************
*   Function   : fileCharLookup
*   Description: This function looksup the count of the character in the F
*                column in the BW-Matrix in the occurrence table at a 
*                specified index.
*   Parameters : fle - The opened BWT file.
*                fle_size - The size of the BWT file.
*                lTab - The occurrence table that keeps track of the count
*                of each character in the alphabet at a particular checkpoint.
*                shift_position - The index to lookup in the occurrence
*                table.
*                lookup_char - The character in the F-column in the BW-Matrix.
*   Effects    : Calculates the count of the character in F column of the 
*                BW-Matrix at a particular index in the occurrence table.
*   Returned   : The count of a character at a particular index in the 
*                occurrence table.
***************************************************************************/
long fileCharLookup(int fle, long fle_size, long **lTab, long shift_position, char lookup_char){
    long count_char = 0;
    long final_count = 0;
    long checkpoint = 0;
    long last_index = 0;
    long mod = shift_position % BLOCK_SIZE;
    long div = shift_position / BLOCK_SIZE;

    if (shift_position == fle_size - 1){
        final_count = lTab[getCharIndex(lookup_char)][div + 1];
    }
    if(mod != 0){
        if((last_index = (div + 1) * BLOCK_SIZE) > fle_size){
            last_index = fle_size - 1;
        }

        //long low_diff = labs((div * BLOCK_SIZE) - shift_position);
        //long high_diff = labs(last_index - shift_position);
        long check_diff = 0;
        char result[BLOCK_SIZE];
        //if(low_diff <= high_diff){
            checkpoint = div;
            check_diff = shift_position - (checkpoint * BLOCK_SIZE);
            //result = calloc(BLOCK_SIZE, sizeof(char));
            lseek(fle, ((checkpoint * BLOCK_SIZE) + 1) * sizeof(char), SEEK_SET);
            if(read(fle, &result, BLOCK_SIZE)){};

            for(long v = 0; v < check_diff; v++){
                if(result[v] == lookup_char){
                    count_char += 1;
                }
            }

            final_count = lTab[getCharIndex(lookup_char)][div] + count_char;
        //}
        //else{
            //checkpoint = div + 1;
            //check_diff = last_index - shift_position;
            //result = calloc((check_diff + 1), sizeof(char));
            //lseek(fle, shift_position * sizeof(char), SEEK_SET);
            //read(fle, result, check_diff);

            //for(long v = 0; v < check_diff; v++){
                //if(result[v] == lookup_char){
                      //count_char += 1;
                  //}
              //}

              //if(characterLookup(fle, last_index) != lookup_char){
                  //final_count = lTab[getCharIndex(lookup_char)][div + 1] - count_char + 1;
              //}
              //else{
                  //final_count = lTab[getCharIndex(lookup_char)][div + 1] - count_char;
              //}
         //}
        // free(result);
    }
    else{
        final_count = lTab[getCharIndex(lookup_char)][div];
    }

    return final_count;
}

/***************************************************************************
*   Function   : bwtDecode
*   Description: This function decode the BWT string using the observation 
*                outlined by Ferragina and Manzini (2000) using the C table
*                and the occurrence table
*   Parameters : fl - The opened BWT file.
*                fl_out - The process id of the opened output file.
*                cntTab - The final count of all the characters in the
*                alphabet. 
*                cTab - The C table constructed from the final count table.
*                lTab - The occurrence table that keeps track of the count
*                of each character in the alphabet at a particular checkpoint.
*                fl_size - The size of the BWT file.
*                term_index - The index of terminating character '\n'.
*   Effects    : Backtracts from the terminating character using the
*                observation outlined by Ferragina and Manzini (2000), where
*                for a specific index in the BWT string you can obtain the
*                preceding character in the F column in the BW-Matrix by using
*                formula:
*                    F[i] = BWT[CT[BWT[i]] + (Occ(BWT[i], i) - 1)
*                In addition, each buffer decoded from the BWT file is written 
*                to the output file.
*   Returned   : NULL
***************************************************************************/
void bwtDecode(int fl_in, int fl_out, long *cntTab, long *cTab, long **lTab, char *bwtTab, long fl_size, long term_index){
    for(int r = 2; r < SYMBOL_COUNT; r++){
        cTab[r] = cTab[r - 1] + cntTab[r - 1];
    }

    long prepend_index = 0;
    char prepend_char = '\0';
    long buffer_len = 0;
    char *origin_string_buffer = calloc(CHUNK_SIZE, sizeof(char));
    char append_char[2];
    append_char[1] = '\0';


    append_char[0] = '\n';
    strncat(origin_string_buffer, append_char, 2);
    buffer_len += 1;

    prepend_index = cTab[getCharIndex('\n')] + (fileCharLookup(fl_in, fl_size, lTab, term_index, '\n') - 1);
    prepend_char = characterLookup(fl_in, bwtTab, prepend_index, fl_size);

    append_char[0] = prepend_char;
    strncat(origin_string_buffer, append_char, 2);
    buffer_len += 1;

    while(1){
        prepend_index = cTab[getCharIndex(prepend_char)] + (fileCharLookup(fl_in, fl_size, lTab, prepend_index, prepend_char) - 1);
        prepend_char = characterLookup(fl_in, bwtTab, prepend_index, fl_size);
        append_char[0] = prepend_char;
        if(prepend_char != '\n' && buffer_len != CHUNK_SIZE - 1){
            strncat(origin_string_buffer, append_char, 2);
            buffer_len += 1;
        }

        if(prepend_index == term_index){
            if(write(fl_out, origin_string_buffer, buffer_len)){};
            break;
        }
        else if(buffer_len == CHUNK_SIZE - 1){
            if(write(fl_out, origin_string_buffer, buffer_len)){};
            memset(origin_string_buffer, 0, CHUNK_SIZE);
            buffer_len = 0;
        }
    }

    free(origin_string_buffer);
}

/***************************************************************************
*   Function   : main
*   Description: This is the main function for this program, it validates
*                the command line input and, if valid, it will call
*                functions to decode a file using the BWT algorithm. 
*   Parameters : argc - number of parameters
*                argv - parameter list
*   Effects    : Decodes BWT file and saves to a specified output file. 
*   Returned   : 0 for success
***************************************************************************/
int main (int argc, char **argv){
    //char *buffer = malloc(sizeof(char) * CHUNK_SIZE);
    char buffer[CHUNK_SIZE];
    // FILE *file;
    
    size_t n_read;
    long file_size = 0;
    int file_in, file_out;
    long *count_table = calloc(SYMBOL_COUNT, sizeof(long));
    long *c_table = calloc(SYMBOL_COUNT, sizeof(long));
    long **occ_l_table = setTables();
    long current_size = 0;
    long term_char_index;
    struct stat file_stat = {};

    c_table[1] = 1;
    long *lf_table;
    long lTab_index = 0;
    long bwt_index = 0;

    if(argc == 3){
        // file = fopen(argv[1], "r");
        file_in = open(argv[1], O_RDONLY);
        file_out = open(argv[2], O_WRONLY | O_CREAT, 0777);
        if(file_in && file_out){
            // fseek(file, 0, SEEK_END);
            // file_size = ftell(file);
            // fseek(file, 0, SEEK_SET);
            fstat(file_in, &file_stat);
            file_size = file_stat.st_size;

            lf_table = malloc(MAX_TABLE_SIZE * sizeof(long));

            char *bwt_array = calloc(MAX_FILE_TABLE_SIZE, sizeof(char)); 

            while((n_read = read(file_in, buffer, CHUNK_SIZE)) > 0){
                chunkProcessing(buffer, file_size, &lTab_index, &bwt_index, current_size, n_read, &term_char_index, count_table, occ_l_table, bwt_array);
                current_size += n_read;
            }
            bwtDecode(file_in, file_out, count_table, c_table, occ_l_table, bwt_array, file_size, term_char_index);
            close(file_in);
            close(file_out);

            char sys_buffer[SYS_CHUNK];
            snprintf(sys_buffer, SYS_CHUNK, "rev %s | sponge %s", argv[2], argv[2]);
            if(system(sys_buffer)){};
            memset(sys_buffer, 0, SYS_CHUNK);

            snprintf(sys_buffer, SYS_CHUNK, "tac %s | sponge %s", argv[2], argv[2]);
            if(system(sys_buffer)){};

            free(lf_table);
            free(bwt_array);
        }
    }

    //free(buffer);
    free(count_table);
    free(c_table);
    freeTables(occ_l_table);

    return 0;
}
