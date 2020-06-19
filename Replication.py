"""
    Pattern Count Problem: Count all occurrences of a pattern in a string.
        Input: A DNA string Pattern.
        Output: The number of occurrences of a pattern in the DNA string.
"""


def PatternCount(Text, Pattern):
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
