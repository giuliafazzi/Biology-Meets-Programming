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


"""
    Input:  String Text and profile matrix Profile
    Output: Pr(Text, Profile)
"""


def Pr(Text, Profile):
    prob = 1

    for i, symbol in enumerate(Text):
        prob *= Profile[symbol][i]

    return prob


"""
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
        Input: A string Text, an integer k, and a 4 x k matrix Profile.
        Output: A Profile-most probable k-mer in Text.
"""


def ProfileMostProbableKmer(Text, k, Profile):
    pr = {}
    most_prob = []
    n = len(Text)

    for i in range(n-k+1):
        k_mer = Text[i:i+k]
        probability = Pr(k_mer, Profile)
        pr[k_mer] = probability

    m = max(pr.values())

    for i, value in pr.items():
        if value == m:
            most_prob.append(i)

    return most_prob[0]


"""
    Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
    Output: GreedyMotifSearch(Dna, k, t)
"""


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []

    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])

    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


"""
    Input:  A set of kmers Motifs
    Output: CountWithPseudocounts(Motifs)
"""


def CountWithPseudocounts(Motifs):
    count = {}
    t = len(Motifs)
    k = len(Motifs[0])

    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    for symbol in "ACGT":
        for j in range(k):
            count[symbol][j] += 1

    return count


"""
    Input:  A set of kmers Motifs
    Output: ProfileWithPseudocounts(Motifs)
"""


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    total = 0
    profile = {}
    count = CountWithPseudocounts(Motifs)
    profile = count

    for symbol in "ACGT":
        total += count[symbol][0]

    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j] /= float(total)

    return profile


"""
    Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
    Output: GreedyMotifSearch(Dna, k, t)
"""


# Include these lines to have WithPseudocounts versions be used:
Count = CountWithPseudocounts
Profile = ProfileWithPseudocounts


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    return GreedyMotifSearch(Dna, k, t)
