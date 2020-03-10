import shlex
import sys
from multiprocessing import Process, Queue, Manager

from utils import subprocess_popen, is_file_exists
from read import Read


def gather_result(q, pipe, has_tag=None):
    for line in pipe:
        read = Read(line)

        if type(has_tag) == list:
            flag = True
            for tag in has_tag:
                if not tag in read.tags:
                    flag = False
                    break
            if flag:
                q.put(read)

        elif type(has_tag) == str:
            if has_tag in read.tags:
                q.put(read)

        else:
            q.put(read)


class SamtoolsViewViewer:

    def __init__(self, has_tag=None, bam="", samtools="samtools", ctg_list=[], ctg_range=None, mapq_filter=10):
        if not is_file_exists(bam):
            print("... Error: Bam file does not exist, please check the input file path ...")
            sys.exit(1)

        if not is_file_exists(bam + ".bai"):
            print("... Bam file is not indexed yet, the index file will now be generated ...")

            bam_index = subprocess_popen(shlex.split("%s index %s" % (samtools, bam)))
            ret_code = bam_index.wait()

            if ret_code != 0:
                print("... Error: An error occured during indexing bam file, please check if the bam file is valid ...")
                sys.exit(1)

        queue = Manager().Queue()
        process_list = []
        subprocess_list = []

        for contig in ctg_list:
            subprocess = subprocess_popen(
                shlex.split(
                    f'{samtools} view -q {mapq_filter} {bam} {contig}{f":{ctg_range}" if ctg_range else ""}'
                )
            )
            subprocess_list.append(subprocess)
            process = Process(target=gather_result, args=(queue, subprocess.stdout, has_tag))
            process_list.append(process)

        for process in process_list:
            process.start()
        for process in process_list:
            process.join()
        for subprocess in subprocess_list:
            subprocess.stdout.close()

        # self.reads = list(queue.queue)
        self.reads = []
        while not queue.empty():
            self.reads.append(queue.get())
            if len(self.reads) % 100000 == 0:
                print(f'... {len(self.reads)} reads processed ...')

        print(f'... Finished reading bam, extracted {len(self.reads)} reads ...')
