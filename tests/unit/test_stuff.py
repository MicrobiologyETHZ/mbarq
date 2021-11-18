# import tnseq.sequence
# import tnseq.commandline
# import collections
#
#
#
# def mutate(sequenceToMutate):
#
#     mutatedSequences = []
#
#     length = len(sequenceToMutate)
#     for pos in range(0, length):
#         prefix = sequenceToMutate[:pos]
#         suffix = sequenceToMutate[pos + 1:]
#         mutatedSequences.append(prefix + 'A' + suffix)
#         mutatedSequences.append(prefix + 'C' + suffix)
#         mutatedSequences.append(prefix + 'T' + suffix)
#         mutatedSequences.append(prefix + 'G' + suffix)
#     return mutatedSequences
#
# k = 4
# kmer_2_pos = collections.defaultdict(list)
# stream = tnseq.sequence.stream_fa('../test/Salmonella_genome_FQ312003.1_SL1344.fasta.gz')
#
# for line in stream:
#     sequence = line.sequence
#     for startpos in range(0, len(sequence) - k):
#         kmer = sequence[startpos: startpos + k]
#         kmer_2_pos[kmer].append(startpos)
#
# print(len(kmer_2_pos))
#
# stream = src.sequence.stream_fa('/Users/hans/Desktop/tnseq/derep.fasta')
#
# for line in stream:
#     sequence = line.sequence
#     kmers = []
#     for startpos in range(0, len(sequence) - k):
#         kmer = sequence[startpos: startpos + k]
#         kmers.append(kmer)
#     print(kmers)
#     break
#
#
