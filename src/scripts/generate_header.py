from sys import stdin


def date_str_for_header():
    from datetime import datetime
    current_date = datetime.today()
    date_str = (
        current_date.strftime('%Y-%m-%d') +
        f'|{current_date.hour}:{current_date.minute}:{current_date.second}PM|HKT|+0800'
    )

    return f'##fileDate={date_str}'


def main():
    from textwrap import dedent

    print('##fileformat=VCFv4.2')
    print(date_str_for_header())
    print('##source=SENSV-v1.0.1')

    for row in stdin.readlines():
        columns = row.strip().split("\t")
        contig_name, contig_size = columns[0], columns[1]
        # make sure no chr prefix here
        if contig_name[:3] == "chr":
            contig_name = contig_name[3:]
        print("##contig=<ID=%s,length=%s>" % (contig_name, contig_size))

    print(dedent("""\
        ##ALT=<ID=DEL,Description="Deletion">
        ##ALT=<ID=INV,Description="Inversion">
        ##ALT=<ID=DUP,Description="Duplication">
        ##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
        ##ALT=<ID=DUP:INT,Description="Interspersed Duplication">
        ##ALT=<ID=INS,Description="Insertion">
        ##ALT=<ID=BND,Description="Breakend">
        ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
        ##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description="Genomic origin of interspersed duplication seems to be deleted">
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
        ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
        ##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">
        ##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description="Standard deviation in span of merged SV signatures">
        ##INFO=<ID=STD_POS,Number=1,Type=Float,Description="Standard deviation in position of merged SV signatures">
        ##INFO=<ID=STD_POS1,Number=1,Type=Float,Description="Standard deviation of breakend 1 position">
        ##INFO=<ID=STD_POS2,Number=1,Type=Float,Description="Standard deviation of breakend 2 position">
        ##FILTER=<ID=q5,Description="Score below 5">
        ##FILTER=<ID=hom_ref,Description="Genotype is homozygous reference">
        ##FILTER=<ID=not_fully_covered,Description="Tandem duplication is not fully covered by a single read">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample"""
    ))


if __name__ == "__main__":
    # usage: less "<fai_absolute_file_path>" | python generate_header.py
    # no chr prefix would be outputted
    main()
