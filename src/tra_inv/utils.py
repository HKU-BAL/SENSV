import shlex
from os.path import isfile
from sys import stderr
from subprocess import PIPE, Popen, DEVNULL
from statistics import mean

major_contigs = (
    ["chr" + str(a) for a in (list(range(1, 23)) + ["X", "Y"])] +
    [str(a) for a in list(range(1, 23))] + ["X", "Y"]
)
major_contigs_set = set(major_contigs)


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def is_file_exists(path):
    return isfile(path)


def get_contig_list(bam, samtools):
    print("... Reading contig list ...")

    contig_list = []

    with subprocess_popen(shlex.split("%s view -H %s" % (samtools, bam))) as p:
        for row in p.stdout:
            columns = row.strip().split(maxsplit=2)
            if columns[0] != "@SQ":
                continue

            contig = columns[1].split(":", maxsplit=2)[1]
            if contig in major_contigs_set:
                contig_list.append(contig)

    print("... %d contigs extracted ..." % len(contig_list))

    return contig_list


def get_reference_from_fasta(ref, chrom, left_pos, right_pos):
    cmd = ['samtools', 'faidx', ref, '%s:%s-%s' % (chrom, left_pos, right_pos)]

    reads = []
    with subprocess_popen(cmd, stderr=DEVNULL) as process:
        for row in process.stdout:
            reads.append(row.strip())

    return "".join(reads[1:])


def depth_stat_from(bam, samtools, breakpoint, left_pos, right_pos):
    left_depths, right_depths = [], []
    chromosome, position = breakpoint.chr, breakpoint.position

    idx = left_pos
    with subprocess_popen(shlex.split("%s depth -a -aa -r %s:%d-%d %s" % (samtools, chromosome, left_pos, right_pos-1, bam))) as p:
        for row in p.stdout:
            columns = row.strip().split(maxsplit=2)

            if idx < position:
                left_depths.append(int(columns[2]))
            else:
                right_depths.append(int(columns[2]))
            idx += 1

    return mean(left_depths), mean(right_depths)
