"""
    Pattern Count Problem: Count all occurrences of a pattern in a string.
        Input: A DNA string Pattern.
        Output: The number of occurrences of a pattern in the DNA string.
"""


def PatternCount(Pattern, Text):
    count = 0

    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1

    return count


"""
Frequent Words Problem:  Find the most frequent k-mers in a string.
     Input: A string Text and an integer k.
     Output: All most frequent k-mers in Text.
"""


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)

    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0

    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] += 1

    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())

    for key in freq:
        if freq[key] == m:
            words.append(key)

    return words


"""
    Reverse Complement Problem: Find the reverse complement of a DNA string.
        Input: A DNA string Pattern.
        Output: The reverse complement of Pattern.
"""


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)

    return Pattern


def Reverse(Pattern):
    rev = ""

    for char in Pattern:
        rev = char + rev

    return rev


def Complement(Pattern):
    com = ""

    for char in Pattern:
        if char == 'A':
            com += 'T'
        elif char == 'T':
            com += 'A'
        elif char == 'C':
            com += 'G'
        elif char == 'G':
            com += 'C'

    return com


"""
    Pattern Matching Problem: Find all occurrences of a pattern in a string.
        Input: Strings Pattern and Genome.
        Output: All starting positions in Genome where Pattern appears as a substring.
"""


def PatternMatching(Pattern, Genome):
    positions = []

    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)

    return positions


"""
    Symbol Array Problem: Find the number of occurrences of a symbol encountered in each window of the genome
        Input: Symbol and Genome.
        Output: Symbol array of Genome corresponding to symbol.
"""


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])

    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        array[i] = array[i-1]

        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1

    return array


"""
    Skew Array Problem: Find the skew array of genome.
        Input: String Genome.
        Output: The skew array of Genome as a list.
"""


def SkewArray(Genome):
    skew = [0]

    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            skew.append(skew[i])

        elif Genome[i] == 'G':
            skew.append(skew[i] + 1)

        elif Genome[i] == 'C':
            skew.append(skew[i] - 1)

    return skew


"""
    Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
        Input: A DNA string Genome.
        Output: All integer(s) i minimizing Skew[i] among all values of i (from 0 to len(Genome)).
"""


def MinimumSkew(Genome):
    positions = []
    count = 0

    skew = SkewArray(Genome)
    minimum = min(skew)

    for i in skew:
        if i == minimum:
            positions.append(count)
        count += 1

    return positions


"""
    Hamming Distance Problem:  Compute the Hamming distance between two strings.
        Input: Two strings of equal length.
         Output: The Hamming distance between these strings.
"""


def HammingDistance(p, q):
    count = 0

    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1

    return count


"""
    Approximate Pattern Matching Problem:  Find all approximate occurrences of a pattern in a string.
        Input: Strings Pattern and Text along with an integer d.
        Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
"""


def ApproximatePatternMatching(Pattern, Text, d):
    positions = []

    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)

    return positions


def ApproximatePatternCount(Text, Pattern, d):
    count = 0

    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1

    return count
