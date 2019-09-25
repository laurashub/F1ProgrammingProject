import DNABindingMap as dbm

DNA = "aaaaaaaaaagtcgcagcgtgggaccgtagctgaGTaattaCGgcagcgcac"
seq = dbm.DNABindingMap()
seq.set_sequence(DNA)
seq.findBindingProteins()
print(seq.get_sequence())
seq.showResults()


