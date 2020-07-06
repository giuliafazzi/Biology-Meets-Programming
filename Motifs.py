import random

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


"""
    Input:  A profile matrix Profile and a list of strings Dna
    Output: Motifs(Profile, Dna)
"""


def Motifs(Profile, Dna):
    t = len(Dna)
    most_prob = []

    for i in range(t):
        most_prob.append(ProfileMostProbableKmer(
            Dna[i], len(Profile['A']), Profile))

    return most_prob


"""
    Input:  A list of strings Dna, and integers k and t
    Output: RandomMotifs(Dna, k, t)
"""


def RandomMotifs(Dna, k, t):
    s = len(Dna[0])
    random_mer = []

    for i in range(t):
        r = random.randint(1, s - k)
        random_mer.append(Dna[i][r:r+k])

    return random_mer


"""
    Input:  Positive integers k and t, followed by a list of strings Dna
    Output: RandomizedMotifSearch(Dna, k, t)
"""


def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


"""
    Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
    Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
"""


def Normalize(Probabilities):
    total = 0
    normalized = {}

    for symbol in Probabilities:
        total += Probabilities[symbol]

    for symbol in Probabilities:
        normalized[symbol] = Probabilities[symbol] / total

    return normalized


"""
    Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
    Output: A randomly chosen k-mer with respect to the values in Probabilities
"""


def WeightedDie(Probabilities):
    kmer = ''
    r = random.uniform(0, 1)
    count = 0

    for key in Probabilities:
        count += Probabilities[key]

        if r < count:
            kmer = key
            return kmer

    return kmer


"""
    Input:  A string Text, a profile matrix Profile, and an integer k
    Output: ProfileGeneratedString(Text, profile, k)
"""


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}

    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


"""
    Input:  Integers k, t, and N, followed by a collection of strings Dna
    Output: GibbsSampler(Dna, k, t, N)
"""


def GibbsSampler(Dna, k, t, N):
    BestMotifs = []
    motifs = RandomMotifs(Dna, k, t)
    BestMotifs = motifs

    for j in range(N):
        i = random.randint(0, t - 1)
        del motifs[i]
        prof = ProfileWithPseudocounts(motifs)
        prof_string = ProfileGeneratedString(Dna[i], prof, k)
        motifs.insert(i, prof_string)

        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs

    return BestMotifs
