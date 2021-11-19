'''Efficient Implementation of the pairwise sequence alignment algorithm in pure Python (3.8+)'''

from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Dict, List, Tuple
)


@dataclass
class SubstitutionMatrix:
    
    matrix: Dict[str, Dict[str, int]] = field(default_factory = dict)
    matrix_type: str = field(default = 'BLOSUM62', compare = True)
    
    @property
    def allowed_matrices(self) -> List[str]:
        """List of allowed substitution matrices."""
        return ['BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'BLOSUM100', 'BLOSUM45',
                'pam30', 'pam70', 'pam250']
    
    @property
    def K(self):
        """Return the available substitution matrix options."""
        return self.matrix.keys()
    
    def __getitem__(self, key: List[str, str]) -> int:
        """Parse letters directly to the main object interface."""
        return self.matrix[key[0]][key[1]]
    
    def _content_parser(self, file_path: Path, splitter: str) -> None:
        """Parse the substitution matrix file."""
        with file_path.open('r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip().split(splitter)
                    self.matrix[line[0]][line[1]] = int(line[2])

    @classmethod
    def from_csv(cls, file_path: Path, matrix_type = 'BLOSUM62') -> 'SubstitutionMatrix':
        """Load a substitution matrix from a CSV file."""
        matrix = cls()
        if matrix_type not in matrix._allowed_matrices:
            raise ValueError(f'{matrix_type} is not a valid substitution matrix.')
        matrix._content_parser(file_path, ',')
        return matrix
    
    @classmethod
    def from_txt(cls, file_path: Path, matrix_type = 'BLOSUM62'):
        """Load a substitution matrix from a text file."""
        matrix = cls()
        if matrix_type not in matrix._allowed_matrices:
            raise ValueError(f'{matrix_type} is not a valid substitution matrix.')
        matrix._content_parser(file_path, '\t')
        return matrix


class PairwiseAlignment:
    
    def __init__(self, 
                 strategy: str = "global", 
                 gap_penalty: int = -2, 
                 substitution_matrix: SubstitutionMatrix | None = None):
        self.strategy = strategy if self._valid_strategy(strategy) else "global"
        self.gap_penalty = gap_penalty
        if substitution_matrix is None:
            self._look_for_matrices()
        else:
            self.sub_matrix = substitution_matrix
    
    @property
    def M(self):
        if not hasattr(self, 'sequence_1'):
            raise AttributeError(
                'No alignment has been performed so no Matrix has been constructed.')
        return self._score_matrix()
    
    def _valid_strategy(self, 
                        strategy: str) -> bool:
        """Check if the strategy is valid."""
        return strategy in ['global', 'semiglobal', 'local']
    
    def _look_for_matrices(self) -> None:
        """Look for substitution matrices in the current directory."""
        for possible_matrix in SubstitutionMatrix.allowed_matrices:
            try:
                
                self.sub_matrix = SubstitutionMatrix.from_txt(
                    Path(__file__).parent.joinpath(f'{possible_matrix}.txt')
                )
                
                return
            
            except FileNotFoundError: pass
    
    def _score_matrix(self) -> List[List[int]]:
        """Construct a score_matrix template, where
        seq1 = X = i and seq2 = Y = j"""
        M = [ [0]*self.X for _ in range(self.Y) ]
        
        if self.strategy == 'global':
            M[0] = [self.gap_penalty * i for i in range(self.X)]
            [M[j].insert(0, self.gap_penalty * j) for j in range(self.Y)]         
        
        return M
    
    def _matrix_pathing(self, 
                        M: List[List[int]], 
                        i: int, 
                        j: int) -> Tuple[int, int]:
        """Return the maximum score and the path taken to reach it."""
        directions = [
            (M[i - 1][j] + self.gap_penalty,   3),  # High road
            (M[i - 1][j - 1] + self.sub_matrix[
                self.sequence_1[i - 1]][self.sequence_2[j - 1]], 2),  # Middle road
            (M[i][j - 1] + self.gap_penalty,   1),  # Low road
        ] 
        
        if self.strategy == 'local': directions += [(0, 0)]
        
        # Current priority is the high-road (value `3` in directions)
        return max(directions)
    
    def _get_max_score_index(self, 
                       M: List[List[int]]) -> Tuple[int, int]:
        """Index the position of a max value."""
        MAX = max([max(row) for row in M])
        for idx, row in enumerate(M):
            # Reverse the list such that we get the first occurance of the value
            # as seen from the back of the list.
            # We keep iterating top-down due to the high-road priority.
            if MAX in row: return idx, -( ( row.reverse().index(MAX) ) + 1 )
    
    def _traceback(self, 
                   M: List[List[int]],
                   coordinates: Tuple[int, int],
                   alignment_1: Tuple[str, str],
                   alignment_2: Tuple[str, str],) -> Tuple[str, str]:
        """Using matrix M, produce the alignment strings when starting from the coordinates."""
        seq1, aligned_seq1, seq2, aligned_seq2 = alignment_1, alignment_2
        i, j = coordinates  # reassign to assure we always have the correct `i` and `j`.
        while i > 0 and j > 0:
            
            # We can simply reuse the function as it looks at all
            # previous cells to determine the path. The direction
            # has been given in the second element of the tuple.
            path = self._matrix_pathing(M, i, j)[1]
            
            if path == 0:
                # This `0` can only be a hard-zero.
                break
            
            elif path == 1:  # Lower road
                j -= 1
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = seq2[j] + aligned_seq2
            
            elif path == 2:  # Middle road
                i -= 1
                j -= 1
                aligned_seq1 = seq1[i] + aligned_seq1
                aligned_seq2 = seq2[j] + aligned_seq2
            
            elif path == 3:  # High road
                i -= 1
                aligned_seq1 = seq1[i] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2            

        if self.strategy != 'local':
            # local only requires the direct alignment.
            # We fill the rest according to the left-over value in either i or j.
            # The other value (= 0) will not add any gap symbols or sequence elements.
            
            # For however many seq1 elements are left, we need to add the same amount of gaps
            # to the other sequence and visa versa:
            aligned_seq1 = ('-' * j) + (seq1[:i] + aligned_seq1)
            aligned_seq2 = ('-' * i) + (seq2[:j] + aligned_seq2)
        
        return (aligned_seq1, aligned_seq2)
    
    def _align(self, 
               sequence_one: str, 
               sequence_two: str) -> Tuple[str, str, int]:
        sequence_1, sequence_2 = sequence_one, sequence_two
        self.X, self.Y = len(sequence_1) + 1, len(sequence_2) + 1
        
        # Create a matrix template:
        M = self._score_matrix()
        
        # Fill in the Matrix M:
        [
            M[i].insert(j, self._matrix_pathing(M, i, j)[0]) 
            for i in range(1, self.X) 
            for j in range(1, self.Y)
        ]
        
        # Set Starting point:
        START, alignment_1, alignment_2 = (), '', ''
        
        if self.strategy == 'global': START = (self.X - 1, self.Y - 1)
        
        elif self.strategy == 'semiglobal':
            START = (self.X - 1, self.Y - 1)
            MAX = I = -1
            
            for idx, row in enumerate(M):
                if row[-1] > MAX:
                    I, MAX = idx, row[-1]
                    
            J = self.Y - 1
            
            if max(M[-1]) < MAX:
                
                for idx in range(M, I + 1, -1):
                    alignment_1 = sequence_1[idx - 2] + alignment_1
                    alignment_2 += '-'
            else:  
                I = -1
                MAX = max(M[-1])
                
                while M[-1][J] != MAX:
                    J -= 1
                    alignment_1 += '-'
                    alignment_2 = sequence_2[J] + alignment_2
            
            # Depending on the if-else, I and J will be assigned accordingly.
            START = (I, J)
        
        else: START = self._get_max_score_index(M)  # Local alignment
        
        # Return -> (seq1, seq2, score)
        return self._traceback(M, START, alignment_1, alignment_2), M[START[0]][START[1]]