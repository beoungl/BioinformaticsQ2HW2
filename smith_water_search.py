import numpy as np

#penalty box
match    = 3
mismatch = -3
gap      = -2

#sequences. seq2 should always be larger than seq1


#Original DNA
seq1 = 'TGTTACGG'
seq2 = 'GGTTGACTA'

#Human DNA
#seq1 = 'CACCCACGC'
#seq2 = 'CTGACAAAGCCT'

#Mouse DNA
#seq1 = 'CACCCACGC'
#seq2 = 'CTGACAAAGCCC'

#Tilapia DNA
seq1 = 'CACCCACTC'
seq2 = 'CACACAAAGCCT'

#Pufferfish DNA
#seq1 = 'CACCCACTC'
#seq2 = 'GAGACAAAGCCT'

#scoring matrix size
rows = len(seq1) + 1

cols = len(seq2) + 1

def create_score_matrix(rows, cols):
    '''
    Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            score_matrix[i][j] = score
    assert max_pos is not None, 'the x, y position with the highest score was not found'
    return score_matrix, max_pos

def calc_score(matrix, x, y):
    #Calculate score for a given x, y position in the scoring matrix.
    #The score is based on the up, left, and upper-left diagnol neighbors
    similarity = match if seq1[x - 1] == seq2[y - 1] else mismatch
    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap
    return max(0, diag_score, up_score, left_score)

def print_matrix(matrix):
    '''
    Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    print(np.matrix(matrix).T)

#add your function(s) to find a solution here.

def searching(matrix,start_x,start_y):
    #Search return position of the max points.
    starting = (start_x,start_y,'start')
    trace_list = []
    trace_list.append(starting)
    while matrix[start_x][start_y] != 0:
        start_x,start_y,pos = traceback(matrix,start_x,start_y)
        trace_list.append((start_x,start_y,pos))
    return trace_list
            
    

def traceback(matrix, x,y):
    #Gets the next position
    diagonal = matrix[x-1][y-1]
    up = matrix[x][y-1]
    left = matrix[x-1][y]
    next_match = max((diagonal,up,left))
    if next_match == diagonal:
        return x-1,y-1,'match'
    elif next_match == up:
        return x,y-1,'gap'
    elif next_match == left:
        return x-1,y,'gap'


def print_pos(list_trace):
    for i in list_trace:
        if (list_trace.index(i)+1) == len(list_trace):
            print(' i = ',i[0],' j = ',i[1],'\n')
            break
        print('i = ',i[0],' j = ',i[1],'-> ', end = '')


def align_string(trace_list, seq1, seq2):
    #Prints out resulting sequence
    sequence_1 = ''
    sequence_2 = ''
    for i in trace_list[::-1]:
        if i[2] == 'match':
            sequence_1 += seq1[i[0]]
            sequence_2 += seq2[i[1]]
        elif i[2] == 'gap':
            sequence_1 += '-'
            sequence_2 += seq2[i[1]]
    print('The resutling sequences')
    print('sequence 1 ',sequence_1)
    print('sequence 2 ',sequence_2)

#end of your function(s)

if __name__ == '__main__':
    #my main
    score_matrix, start_pos = create_score_matrix(rows, cols)
    start_x, start_y = start_pos

    print_matrix(score_matrix)
    list_trace = searching(score_matrix,start_x,start_y)
    #print(list_trace)
    print_pos(list_trace)
    align_string(list_trace, seq1,seq2)

