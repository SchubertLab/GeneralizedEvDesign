from __future__ import print_function
import numpy as np
import random
import subprocess
import os
import time
import click
import pandas as pd
import pymc3 as pm

@click.group()
def cli():
    pass


@cli.command()
def analyze_results():
    def do_plots(summary_df):
        sns.set()

        for col, marker, color in zip(summary_df.columns, 'o^s', '012'):
            peptides, _, time = summary_df[col].dropna().reset_index().values.T
            plt.plot(peptides, time, 'C%s%s' % (color, marker), label=' '.join(col))

        plt.legend()
        plt.xlabel('Candidate Peptides')
        plt.ylabel('Time (min.)')
        plt.title('Total Time to Design a 35 Aminoacids Vaccine')
        plt.tight_layout()
        plt.savefig('./dev/benchmark.png')

    df = pd.read_csv('./dev/benchmark.csv')
    df['total'] = df['total'] / 60.0
    summary_df = df.pivot_table('total', ['peptides', 'repeat'], ['method', 'constraints'])

    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception as exc:
        print('Cannot import matplotlib or seaborn, will not produce plots! Error is:', exc)
    else:
        do_plots(summary_df)


@cli.command()
def process_results():
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

    results = []
    for fname in os.listdir('./dev/'):
        run_info = parse_file_name(fname)
        if not run_info:
            continue

        stopwatch = parse_log_file(fname)
        if not stopwatch:
            continue

        run_info.update(stopwatch)
        results.append(run_info)

    df = pd.DataFrame(results)
    df.to_csv('./dev/benchmark.csv', index=False)


@cli.command()
def run_benchmark():
    def get_new_file_name(method, size):
        trials = [
            int(fname[:-4].split('-')[5])
            for fname in os.listdir('./dev')
            if fname.startswith('%s-random-%d-repeat-' % (method, size))
        ]

        index = (max(trials) + 1) if trials else 0
        return '%s-random-%d-repeat-%d.log' % (method, size, index)

    base_args = [
        'python', 'mosaic_test.py', 'resources/hivgen.fasta',
        '-v', '-a', '35', '-T'
    ]

    while True:
        size = random.randint(5, 31)
        lazy = random.choice([None, 'mtz', 'dfj'])

        args = base_args + ['-r', str(size)]
        if lazy is not None:
            args.extend(['-l', lazy])
            method = lazy + '-lazy'
        else:
            method = 'mtz-greedy'
        
        fname = './dev/' + get_new_file_name(method, size)
        print(' '.join(args), '>', fname)

        try:
            with open(fname, 'w') as fout:
                start_time = time.time()
                subprocess.call(args, stdout=fout, stderr=fout)
                end_time = time.time()
        except KeyboardInterrupt:
            print('Interrupted - removing incomplete log file!')
            os.remove(fname)
            return


if __name__ == '__main__':
    cli()