import subprocess
import csv
import os
import time
import click


@click.group()
def cli():
    pass


@cli.command()
def process_result():
    def parse_file_name(fname):
        if not fname.endswith('.log'):
            return None

        parts = fname[:-4].split('-')
        if len(parts) != 6:
            return None
        
        if parts[2] != 'random' or parts[4] != 'repeat':
            return None
        
        return {
            'constraints': parts[0],    # dfj or mtz
            'method': parts[1],         # lazy or greedy
            'size': parts[3],
            'repeat': parts[5], 
        }
    
    def parse_log_file(fname):
        info = {}
        with open('./dev/' + fname) as f:
            stopwatch = None
            for row in f:
                if 'peptides above threshold, breakdown' in row:
                    info['peptides'] = int(row.split()[0])
                elif '==== Stopwatch' in row:
                    stopwatch = []
                elif stopwatch is not None:
                    duration = float(row.split(':')[1].split()[0])
                    stopwatch.append(duration)
        
        if not stopwatch or not info:
            return None

        info.update({
            'total': stopwatch[0],
            'preproc': stopwatch[1],
            'creation': stopwatch[2],
            'solving': stopwatch[3],
        })
        return info

    writer = None
    with open('./dev/benchmark.csv', 'w') as f:
        for fname in os.listdir('./dev/'):
            run_info = parse_file_name(fname)
            if not run_info:
                continue

            stopwatch = parse_log_file(fname)
            if not stopwatch:
                continue

            run_info.update(stopwatch)
            if writer is None:
                writer = csv.DictWriter(f, run_info.keys())
                writer.writeheader()
            writer.writerow(run_info)


@cli.command()
def run_benchmark():
    base_args = [
        'python', 'mosaic_test.py', 'resources/hivgen.fasta',
        '-v', '-a', '35', '-T'
    ]

    for i in range(50):
        for r in [5, 10, 15, 20, 25, 30]:
            for lazy in [None, 'mtz', 'dfj']:
                args = base_args + ['-r', str(r)]
                if lazy is not None:
                    args.extend(['-l', lazy])
                    method = lazy + '-lazy'
                else:
                    method = 'mtz-greedy'
                
                fname = './dev/%s-random-%d-repeat-%d.log' % (method, r, i)
                print(' '.join(args))
                print(fname)
                with open(fname, 'w') as fout:
                    start_time = time.time()
                    subprocess.call(args, stdout=fout, stderr=fout)
                    end_time = time.time()
                
                print('%d;%d;%s;%.2f' % (r, i, lazy, end_time - start_time))


if __name__ == '__main__':
    cli()