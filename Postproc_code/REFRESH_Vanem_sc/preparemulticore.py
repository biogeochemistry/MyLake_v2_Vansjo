'''this script should be scenario-dependent. command line argument is
the number of CPUs'''
import sys
import math
import os

startdate = '1990-04-01' ## should be 1990-04-01
enddate = '2012-09-25'

ii2 = range(1, 46)
mylakein2 = ['../../../MyLake_inputs/MyLake_Vanemfjorden/c%s_m%s_input' %
             ([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4][si],
              [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3][si]
              )
             for si in range(12)]
mylakeinit = '../../../MyLake_inputs/MyLake_Vanemfjorden/' + \
    'VNM_init_v12_2012-02-15seereadme.txt'
mylakebasepar = '../../../MyLake_inputs/MyLake_Vanemfjorden/' + \
    'VANEM_para_v12_1b_2012-11-19.txt'
mylakeposteriorparameters = '../../../MyLake_inputs/MyLake_Vanemfjorden/' + \
    'posterior_parameters_using_median_INCA-P/uniqueresampledparameters.csv'

ncores = int(sys.argv[1])
nruns = len(ii2) * len(mylakein2)
nchainspercore = int(math.ceil(nruns / float(ncores)))
nshortcores = nchainspercore * ncores - nruns
n = [nchainspercore] * (ncores - nshortcores) + \
    [nchainspercore - 1] * nshortcores
print(n)

## create all combinations
ii = ii2 * len(mylakein2)
mylakein = []
for k in range(len(mylakein2)):
    mylakein += [mylakein2[k]] * len(ii2)
##

shscriptnames = ['core%d.sh' % i for i in range(ncores)]

f = open('startmulticore.sh', 'w')
f.write(' & '.join(['sh %s' % script for script in shscriptnames]) + ' &')
f.write('\n')
f.close()

dirnames = ['runs/%03d-%s' %
            (ii[j], mylakein[j].split('/')[-1].split('_input')[0])
            for j in range(len(ii))]

for dirname in dirnames:
    print(dirname)
    os.mkdir(dirname)
for corei in range(ncores):
    shscriptname = shscriptnames[corei]
    assignedjj = range(sum(n[:corei]) + 1 - 1, sum(n[:corei + 1]) + 1 - 1)
    commands = ['cd %s ; octave ../../../MyLake_scripts/mylake10pars.m %s %s %s %s %s %s %s ; Rscript ../../../MyLake_scripts/interpolatetemperature7depths.R ; cd -' % \
                (dirnames[j],
                 mylakeinit,
                 mylakebasepar,
                 mylakeposteriorparameters,
                 ii[j],
                 mylakein[j],
                 startdate,
                 enddate)
                for j in assignedjj]
    out = ' ;\n\n'.join(commands)
    print('\n\n')
    print(shscriptname)
    print('\n')
    print(out)
    f = open(shscriptname, 'w')
    f.write(out)
    f.write('\n')
    f.close()
