import numpy as np

#penalty box
match    = 3
transversion = -3
gap      = -2
transition = -1
mismatch = -4

#sequences. seq2 should always be larger than seq1


#Original DNA
#seq1 = 'TGTTACGG'
#seq2 = 'GGTTGACTA'


#Human DNA
#seq1 = 'CACCCACGC'
#seq2 = 'CTGACAAAGCCT'

#Mouse DNA
#seq1 = 'CACCCACGC'
#seq2 = 'CTGACAAAGCCC'

#Tilapia DNA
#seq1 = 'CACCCACTC'
#seq2 = 'CACACAAAGCCT'

#Pufferfish DNA
#seq1 = 'CACCCACTC'
#seq2 = 'GAGACAAAGCCT'

seq1 = 'ACTGGAACTGAGAGGCTTAC'
seq2 = 'ACTCCCGGAACGGAGAAGCTTAAAAC'

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
    #similarity = match if seq1[x - 1] == seq2[y - 1] else mismatch
    similarity = 0
    if seq1[x-1] == seq2[y-1]:
        similarity = match
    else:
        if seq1[x-1] == 'A':
            if seq2[y-1] == 'G':
                similarity = transition
            elif seq2[y-1] == 'C' or seq2[y-1] == 'T':
                similarity = transversion
        elif seq1[x-1] == 'T':
            if seq2[y-1] == 'C':
                similarity = transition
            elif seq2[y-1] == 'G' or seq2[y-1] == 'A':
                similarity = transversion
        elif seq1[x-1] == 'G':
            if seq2[y-1] == 'A':
                similarity = transition
            elif seq2[y-1] == 'T' or seq2[y-1] == 'C':
                similarity = transversion
        elif seq1[x-1] == 'C':
            if seq2[y-1] == 'T':
                similarity = transition
            elif seq2[y-1] == 'G' or seq2[y-1] == 'A':
                similarity = transversion
        else:
            similarity = mismatch
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
    print(trace_list)
    while matrix[start_x][start_y] != 0:
        start_x,start_y,sequence_1,sequence_2 = traceback(matrix,start_x,start_y)
        trace_list.append((start_x,start_y,sequence_1,sequence_2))
    return trace_list
            
    

def traceback(matrix, x,y):
    #Gets the next position
    diagonal = matrix[x-1][y-1]
    up = matrix[x][y-1]
    left = matrix[x-1][y]
    next_match = max((diagonal,up,left))
    if next_match == diagonal:
        return x-1,y-1,seq1[x-1],seq2[y-1]
    elif next_match == up:
        return x,y-1,seq1[x-1],'-'
    elif next_match == left:
        return x-1,y,'-',seq2[y-1]


def print_pos(list_trace):
    for i in list_trace:
        if (list_trace.index(i)+1) == len(list_trace):
            print(' i = ',i[0],' j = ',i[1],'\n')
            break
        print('i = ',i[0]
              ,' j = ',i[1],'-> ', end = '')


def align_string(trace_list, seq1, seq2):
    #Prints out resulting sequence
    sequence_1 = ''
    sequence_2 = ''
    seq1_list = []
    seq2_list = []
    match_list = trace_list[1:]
    print(trace_list)
    for i in match_list[::-1]:
        seq1_list.append(i[2])
        seq2_list.append(i[3])
    matching_list = ''.join(seq2_list)
    empty_list = [' ' for i in range(len(seq2))]
    for i in range(len(seq2)):
        if seq2[i] == matching_list[i]:
            break
        else:
            empty_list[i] = '-'
    last = -1
    try:
        last = empty_list[::-1].index('-')
    except ValueError:
        last = -1
    if last != -1:
        for i in range(len(seq2))[last:]:
            empty_list[i] = matching_list[matching_list-last]
    else:
        for i in range(len(seq2)):
            if i < len(matching_list):
                empty_list[i] = matching_list[i]
    for i in range(len(empty_list)):
        if empty_list[i] == ' ':
            empty_list[i] = '-'
    alignment = ''.join(empty_list)
    print(alignment)
    print(seq2)
    align_list = []
    for i,z in zip(alignment,seq2):
        if i == z:
            align_list.append('*')
        elif i == '-':
            align_list.append('-')
        else:
            if i == 'A':
                if z == 'G':
                    align_list.append('.')
                else:
                    align_list.append(':')
            elif i == 'C':
                if z == 'T':
                    align_list.append('.')
                else:
                    align_list.append(':')
            elif i == 'G':
                if z == 'A':
                    align_list.append('.')
                else:
                    align_list.append(':')
            elif i == 'T':
                if z == 'C':
                    align_list.append('.')
                else:
                    align_list.append(':')
            else:
                align_list.append(' ')
    print(''.join(align_list))
    #print('The resutling sequences')
    #print('sequence 1 ',sequence_1)
    #print('sequence 2 ',sequence_2)

#end of your function(s)

if __name__ == '__main__':
    #my main
    score_matrix, start_pos = create_score_matrix(rows, cols)
    start_x, start_y = start_pos

    print_matrix(score_matrix)
    list_trace = searching(score_matrix,start_x,start_y)
    #print(list_trace)
    #print(seq1,seq2)
    #print_pos(list_trace)
    align_string(list_trace, seq1,seq2)

