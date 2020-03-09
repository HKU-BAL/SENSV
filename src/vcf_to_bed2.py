import re
from argparse import ArgumentParser

vcfHeaders = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
outHeaders = ['CHROM', 'POS', 'CHR2', 'SVEND', 'SVTYPE', 'REMARK']

if __name__ == "__main__":
    parser = ArgumentParser(description='detect SV')
    parser.add_argument('-v', '--vcf', help='input vcf file', required=True)
    parser.add_argument('-t', '--defaultSvType', help='default SV type', required=False)
    parser.add_argument('-a', '--ans', help='answer', required=False)

    options = parser.parse_args()

    prevSvEnd = 0
    with open(options.vcf, 'r') as f_vcf:
        for line in f_vcf:
            if line.startswith('#'):
                continue

            result = {}
            fields = line.rstrip().split('\t')

            try:
                infos = fields[7].split(';')
            except:
                print("not enough fields\n", fields)
                continue

            for i, info in enumerate(infos, 1):
                try:
                    key, value = info.split('=')
                except ValueError:
                    key = info
                    value = 1

                result[key] = value

            if 'SVTYPE' in result:
                if result['SVTYPE'] == 'BND':
                    result['SVTYPE'] = 'TRA'
                    arr = re.split(r'\[|\]|:', fields[4])
                    result['CHR2'] = arr[1]
                    result['END'] = arr[2]
            else:
                if fields[4] in ['<DEL>', '<INV>', '<DUP>', '<INS>', '<TRA>', '<TRA_INV>']:
                    result['SVTYPE'] = fields[4][1:-1]
                elif '[' in fields[4] or ']' in fields[4]:
                    arr = re.split(r'\[|\]|:', fields[4])
                    result['SVTYPE'] = 'TRA'
                    result['CHR2'] = arr[1]
                    result['END'] = arr[2]
                else:
                    print('check 1', line)
                    continue

            for i, col in enumerate(vcfHeaders[:7]):
                result[col] = fields[i]

            if options.ans and not result['ANS'] == options.ans:
                continue

            if 'END' in result and not 'SVEND' in result:
                result['SVEND'] = result['END']

            if not 'SVLEN' in result:
                result['REMARK'] = str(int(result['SVEND']) - int(result['POS']))
            else:
                result['REMARK'] = str(abs(int(result['SVLEN'])))

            """
            if result['SVTYPE'] in ['TRA','TRA_INV']:
                try:
                    result['REMARK'] = result['CHR2']
                except:
                    print('check 1', line)
            """

            # if int(result['SVLEN']) < 1000:
            #    continue

            """
            if int(result['POS']) < prevSvEnd + 1000:
                continue

            prevSvEnd = int(result['SVEND'])
            """

            print('\t'.join([result[header] for header in outHeaders]))
