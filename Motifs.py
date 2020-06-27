"""
    Input:  A set of kmers Motifs
    Output: The count matrix of  Motifs, as a dictionary
"""


def Count(Motifs):
    count = {}

    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


"""
    Input:  A list of kmers Motifs
    Output: the profile matrix of Motifs, as a dictionary of lists.
"""


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    profile = count

    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j] = count[symbol][j] / float(t)

    return profile


"""
    Input:  A set of kmers Motifs
    Output: A consensus string of Motifs.
"""


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


"""
    Input:  A set of k-mers Motifs
    Output: The score of these k-mers.
"""


def Score(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    consensus = Consensus(Motifs)
    count = Count(Motifs)
    score = 0

    for j in range(k):
        symbol = consensus[j]
        value = count[symbol][j]
        score += t - value

    return score
