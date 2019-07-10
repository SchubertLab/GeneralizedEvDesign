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
    def do_plots(summary_df, res):
        sns.set()

        for col, marker, color in zip(summary_df.columns, 'o^s', '012'):
            peptides, time = summary_df[col].dropna().reset_index().values.T

            test_xs = np.linspace(300, 2700, 100)
            test_ys = res['time_sd'] * (res[col]['cc'] + res[col]['aa'] * np.exp(
                res[col]['bb'] * (test_xs - res['peptides_mean']) / res['peptides_sd']
             ))

            plt.plot(peptides, time, 'C%s%s' % (color, marker))
            plt.plot(test_xs, test_ys, label=col)

        plt.legend()
        plt.xlabel('Candidate Peptides')
        plt.ylabel('Time (min.)')
        plt.title('Total Time to Design a 35 Aminoacids Vaccine')
        plt.tight_layout()
        plt.savefig('./dev/benchmark.png')

    def exponential_fit(summary_df):
        peptides_mean = summary_df.index.values.mean()
        peptides_sd = summary_df.index.values.std()
        all_times = summary_df.values.ravel()
        time_sd = all_times[~np.isnan(all_times)].std()

        with pm.Model() as model:
            for column in summary_df.columns:
                name = '_'.join(column)

                peptides, time = summary_df[column].dropna().reset_index().values.T

                data_xs = (peptides - peptides_mean) / peptides_sd
                data_ys = time / time_sd

                alpha = pm.Gamma('%s_aa' % name, 3, 2)
                beta = pm.Gamma('%s_bb' % name, 3, 2)
                gamma = pm.Normal('%s_cc' % name, 0, 1)

                zs = pm.Deterministic('%s_zs' % name, gamma + alpha * pm.math.exp(beta * data_xs))
                ys = pm.Normal('%s_ys' % name, zs, pm.InverseGamma('%s_ee_sd' % name, 5, 5),
                               observed=data_ys)

            trace = pm.sample(1000, nuts_kwargs={'target_accept': 0.99}, tune=1000)

        fit_results = {
            'peptides_mean': peptides_mean,
            'peptides_sd': peptides_sd,
            'time_sd': time_sd,
        }

        for column in summary_df.columns:
            name = '_'.join(column)
            fit_results[column] = {
                'aa': trace['%s_aa' % name].mean(),
                'aa_sd': trace['%s_aa' % name].std(),
                'bb': trace['%s_bb' % name].mean(),
                'bb_sd': trace['%s_bb' % name].std(),
                'cc': trace['%s_cc' % name].mean(),
                'cc_sd': trace['%s_cc' % name].std(),
            }

        return fit_results
    
    df = pd.read_csv('./dev/benchmark.csv')
    df['total_minutes'] = df['total'] / 60.0
    # here I group by peptides, so there is no need to deal with repeats since it's unlikely
    # two runs have the same number of peptides. and if that happens we just take the mean, so no big deal
    summary_df = df.pivot_table('total_minutes', 'peptides', ['method', 'constraints'])
    fit_results = exponential_fit(summary_df)
    for col in summary_df.columns:
        print(' '.join(col), ':', 'alpha = %.2f(%.2f) - beta = %.2f(%.2f)' % (
            fit_results[col]['aa'], fit_results[col]['aa_sd'],
            fit_results[col]['bb'], fit_results[col]['bb_sd'],
        ))

    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception as exc:
        print('Cannot import matplotlib or seaborn, will not produce plots! Error is:', exc)
    else:
        do_plots(summary_df, fit_results)


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