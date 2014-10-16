__author__ = 'blais'

count = 0
fasta = open('/home/blais/Downloads/BAX.fasta', 'w')
with open('/home/blais/Downloads/MD6675SE.SEQ') as file:
    for line in file:
        if line[:11] == '; SEQ ID NO':
            fasta.write(('>' + line[2:].rstrip().replace(" ", "_")))
            count = 1
        elif line[:13] == ';   ORGANISM:':
            fasta.write('|' + line[14:].rstrip().replace(" ", "_"))
        elif line[:22] == ';   OTHER INFORMATION:':
            fasta.write('|' + line[23:].rstrip().replace(" ", "_"))
        elif line[:3] == 'CA2':
            fasta.write('|' + line)
        elif line[0] != r';'and count == 1:
            fasta.write(line.replace("1", ""))


