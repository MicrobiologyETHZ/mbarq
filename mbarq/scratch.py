from  Bio.Blast.Applications import NcbimakeblastdbCommandline

import Bio.Application as ba

import subprocess
cline = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file = '/Users/ansintsova/git_repos/mbarq/tests/mbarq_test_data/mapping/test_genome.fna.gz')
                                   #input_file="/Users/ansintsova/git_repos/mbarq/tests/mbarq_test_data/mapping/test_genome.fna")


if __name__ == '__main__':
    try:
        cline()
    except ba.ApplicationError as e:
        print(e.stderr)
        if 'does not match input format type, default input type is FASTA' in e.stderr:
            print('aha')
    #subprocess.Popen(cline, shell=True, stdout=subprocess.PIPE)