import csv
from multiprocessing import cpu_count
from time import sleep
from tempfile import NamedTemporaryFile
import arrow
import psutil
from sh import oligoarray_cl as oligoarray



class OligoarrayDesigner(object):
    def run(self,
            input_file,
            output_file="temp/output.txt",
            rejected_file="temp/rejected.txt",
            log_file="temp/log.txt",
            max_oligos=10000,
            min_length=30,
            max_length=30,
            max_dist=10000,
            min_tm=85,
            max_tm=90,
            secondary_struct_temp=65.0,
            cross_hyb_temp=65,
            min_gc=40,
            max_gc=60,
            prohibited_seqs=[6 * 'C', 6 * 'G', 6 * 'T', 6 * 'A'],
            num_processors=-1,
            min_dist=2,
            timeout=10):
        min_spacing = max_length + min_dist
        prohibited_seqs = ";".join(prohibited_seqs)
        if num_processors == -1:
            from math import ceil
            num_processors = int(ceil(cpu_count() * 2.5))
        passed = False
        start_time = arrow.utcnow()
        # Time after which to error out
        end_time = start_time.replace(minutes=timeout, )
        results = OligoArrayResults()
        call = oligoarray(
            i=input_file,
            d=self.blast_db,
            o=results.file.name,
            r=rejected_file,
            R=log_file,
            n=max_oligos,
            l=min_length,
            L=max_length,
            D=max_dist,
            t=min_tm,
            T=max_tm,
            s=secondary_struct_temp,
            x=cross_hyb_temp,
            p=min_gc,
            P=max_gc,
            N=num_processors,
            g=min_spacing,
            _bg=True,
            m=prohibited_seqs,
            _tty_in=True  # This was required after adding progressbar.  Don't know why.  Fishy
        )
        # This loop doesn't let oligoarray run too long
        # Oligoarray has an issue with rogue blastall processes running forever
        no_blast = False
        sleep_time = 0
        while arrow.utcnow() < end_time:
            # if sleep_time % 60 == 0 and sleep_time != 0:
            #     print("Blast running for: {:.01f} minutes".format(sleep_time /
            #                                                       60))
            if call._process_completed:
                break
            for proc in psutil.process_iter():
                try:
                    if proc.name() == 'blastall':
                        create_time = arrow.get(proc._create_time)
                        # Kill any blast thats run more than two minutes
                        # This number will probably need to get adjusted for large DBs and slow cpus
                        if arrow.utcnow() > create_time.replace(minutes=3,
                                                                seconds=00):
                            proc.kill()
                            print("Killing {}".format(proc))
                except:
                    pass  # Sometimes pid gets killed before this can happen
            # Check that no blasts are running
            # If so, oligoarray probably done
            try:
                # Sometimes this fails
                blast_proc = [proc for proc in psutil.process_iter()
                              if proc.name() == 'blastall']
            except:
                pass
            if not blast_proc:
                # Make sure no blast two times in a row
                if no_blast:
                    break
                no_blast = True
            else:
                no_blast = False
            sleep(2)
            sleep_time += 2
        else:
            call.kill()
            for proc in psutil.process_iter():
                if proc.name() == 'blastall':
                    proc.kill()
            print("Errored out of loop")
        tot_time = (arrow.utcnow() - start_time).total_seconds()
        # print("Design took : {:.02f}".format(tot_time))

        res = results.parse()
        results.close()
        return res

    def __init__(self,
                 blast_db="/home/eric/blastdb/old_format/mouse_refseq_rnaDB"):
        self.blast_db = blast_db


class OligoArrayResults(object):
    def close(self):
        self.file.close()

    def parse(self):
        results = []
        failed = []
        n = -1
        with open(self.file.name, "r") as f:
            tsvin = csv.reader(f, delimiter='\t')
            for n, row in enumerate(tsvin):
                off_targets = row[7].count('; ')
                if off_targets == 0:
                    name = row[0]
                    # Incase probes are chunked into blocks
                    if "=Chunk:" in name:
                        name = name.split('=Chunk:')[0]
                    # Check to make sure designated target is in target list
                    if name not in row[7]:
                        failed.append(row)
                        continue
                    results.append({
                        "target": name,
                        "seq": row[-1],
                    })
                else:
                    failed.append(row)
        return results

    def __init__(self, ):
        self.file = NamedTemporaryFile('w')

