import subprocess
import time


def main():
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
                
                fname = 'dev/%s-random-%d-repeat-%d.log' % (method, r, i)
                print(' '.join(args))
                print(fname)
                with open(fname, 'w') as fout:
                    start_time = time.time()
                    subprocess.call(args, stdout=fout, stderr=fout)
                    end_time = time.time()
                
                print('%d;%d;%s;%.2f' % (r, i, lazy, end_time - start_time))


if __name__ == '__main__':
    main()