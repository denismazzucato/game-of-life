import subprocess
import time
import json
import os.path

################################################################################

mpi_cores = [6]
prun_cores = [8, 16, 32]
cores = prun_cores # TODO: change here
def build_cmd(file_name, param, np):
  return build_mpi(file_name, param, np) # TODO: change here

################################################################################

seq = 'gol-seq'
par = 'gol-par'
opt = 'gol-par-opt'
treshold = 5000 * 5000 * 1000 # se sotto il treshold  fallo 3 times
tests = {
  seq: {
    1: [(1000, 1000, 500)],
  },
  par: {
    3: [(1000, 1000, 500)],
    6: [(1000, 1000, 500)],
  },
  opt: {
    3: [(1000, 1000, 500)],
    6: [(1000, 1000, 500)],
  }
}


################################################################################

def build_prun(file_name, param, np):
  return "prun -1 -np {} -script $PRUN_ETC/prun-openmpi ./gol/{} {} {} {} 0".format(np, file_name, param[0], param[1], param[2])

def build_mpi(file_name, param, np):
  return "mpirun -np {} ./gol/{} {} {} {} 0".format(np, file_name, param[0], param[1], param[2])

def perform_cmd(cmd):
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  retval = p.wait()
  stream = p.stdout.read().decode("utf-8")
  return stream

def build_out_name(file_name, param, np):
  return "./out/test_{}_np{}_{}{}{}_.out".format(file_name, np, param[0], param[1], param[2])

def append(file_name, param, np, ctx):
  out_path = build_out_name(file_name, param, np)
  with open(out_path) as f:
    data = json.load(f)

  new_data = {
    'exec_time': (data['exec_time'] + ctx['exec_time'])/2,
  }
  (data['exec_time'] + ctx['exec_time'])/2

  with open('test.json', 'w') as f:
      json.dump(data, f)

def read(file_name, param, np, ctx):
  in_path = build_out_name(file_name, param, np)
  with open(in_path, 'r') as in_file:
    data = json.load(in_file)
  return data

def build_json(stream):
  lines = list(filter(lambda x: x != '', stream.split('\n')))
  lc = -1
  ex = -1
  for line in lines:
    if "Number" == line.split()[0]:
      lc = line.split()[-1]
    if "Game" == line.split()[0]:
      ex = line.split()[-2]

  return { 'live_cells': lc, 'exec_time': ex }

def perform_test(file_name, param, np):
  output = perform_cmd(build_cmd(file_name, param, np))
  # append output in file
  out_ctx = build_json(output)
  append(file_name, param, np, out_ctx)

def exec_tests(tests, use_saved_data):
  for file_name in tests:
    for np in tests[file_name]:
      for param in tests[file_name][np]:
        perform_test(file_name, param, np)

millis = int(round(time.time() * 1000))
out_ctx = "Output run {}\n".format(millis)

out_obj = {}

for fil in files:
  out_obj[fil] = {}
  for core in cores:
    out_obj[fil][core] = {}
    for param in params:
      print("\nStarted {} {} {} ({} out of {})".format(fil, core, param, i, tests))
      out_obj[fil][core][param] = {}
      for k in range(step):
        out_obj[fil][core][param][k] = {}
        p = subprocess.Popen(build_cmd(fil, param, core), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        retval = p.wait()
        (lc, ex) = extrapolate_info(p.stdout.read())
        out_obj[fil][core][param][k]['live_cells'] = int(lc)
        out_obj[fil][core][param][k]['exec_time'] = float(ex)
        out_ctx += "[{}, {}, {}, {}]\n > (live cells: {}; exec time: {})\n".format(fil, core, param, k, lc, ex)

      i += step
      print("Ok")

fil = 'gol-seq'
out_obj[fil] = {}
core = 1
out_obj[fil][core] = {}
for param in params:
  out_obj[fil][core][param] = {}
  for k in range(step):
    out_obj[fil][core][param][k] = {}
    p = subprocess.Popen(build_cmd(fil, param, core), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    retval = p.wait()
    (lc, ex) = extrapolate_info(p.stdout.read())
    out_obj[fil][core][param][k]['live_cells'] = int(lc)
    out_obj[fil][core][param][k]['exec_time'] = float(ex)
    out_ctx += "[{}, {}, {}, {}]\n > (live cells: {}; exec time: {})\n".format(fil, core, param, k, lc, ex)

  i += step
  print("Ok {} {} {} ({} out of {})".format(fil, core, param, i, tests))

# print(out_obj)

def get_avg_ex(out_obj, fil, core, param):
  sum = 0
  for k in range(step):
    sum += out_obj[fil][core][param][k]['exec_time']
  return sum / step

def get_speedup(out_obj, fil, core, param):
  seq_time = get_avg_ex(out_obj, 'gol-seq', 1, param)
  this_time = get_avg_ex(out_obj, fil, core, param)
  return seq_time / this_time

def get_eff(out_obj, fil, core, param):
  return get_speedup(out_obj, fil, core, param) / core

def get_comm_weight(out_obj, fil, core, param):
  if fil == 'gol-par':
    return param[1]
  else:
    return int(param[1] / 31) + 1

def get_tot_msg(out_obj, fil, core, param):
  return get_comm_weight(out_obj, fil, core, param) * 2 * param[2] * core

def get_problem_size(out_obj, fil, core, param):
  return (int(param[0] / core) + 1) * param[1]

# speed up for cores
stat = {}
for fil in files:
  stat[fil] = {}
  for core in cores:
    stat[fil][core] = {}
    for param in params:
      stat[fil][core][str(param)] = {
        'seq_avg_ex' : get_avg_ex(out_obj, 'gol-seq', 1, param),
        'avg_ex': get_avg_ex(out_obj, fil, core, param),
        'speedup' : get_speedup(out_obj, fil, core, param),
        'efficiency': get_eff(out_obj, fil, core, param),
        'comm_weight': get_comm_weight(out_obj, fil, core, param),
        'tot_msg': get_tot_msg(out_obj, fil, core, param),
        'problem_size': {'total':param[0]*param[1], 'partitioned':get_problem_size(out_obj, fil, core, param)},
      }

# output update
out_path = "./out/test_script_{}.out".format(millis)
out_file = open(out_path, 'w')
out_file.write(out_ctx)
out_file.close()

# stat update
stat_path = "./out/test_script_stat_{}.json".format(millis)
stat_file = open(stat_path, 'w')
json.dump(stat, stat_file)
stat_file.close()