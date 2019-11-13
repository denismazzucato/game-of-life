import subprocess
import time
import os.path
import json

################################################################################

def build_cmd(file_name, param, np):
  return build_prun(file_name, param, np) # TODO: change here

################################################################################

seq = 'gol-seq'
par = 'gol-par'
opt = 'gol-par-opt'
tests = {
  seq: {
    1: [
      (10000, 10000, 5000), # ~33 min
      (5000, 5000, 1000),
      (5000, 5000, 1000),
    ],
  },
  par: {
    1: [(5000, 5000, 1000),(5000, 5000, 1000),]
    8: [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),],
    16: [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),],
    32 : [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),]
  },
  opt: {
    8: [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),],
    16: [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),],
    32: [(10000, 10000, 5000),(5000, 5000, 1000),(5000, 5000, 1000),],
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
  return "./out/test_{}_np{}_{}-{}-{}.data".format(file_name, np, param[0], param[1], param[2])

def append(file_name, param, np, ctx):
  out_path = build_out_name(file_name, param, np)

  if os.path.isfile(out_path):
    with open(out_path, 'r') as f:
      old_data = json.load(f)
  else:
    old_data = { "exec_time": ctx, "times": 0 }
  ex = float(old_data["exec_time"])
  times = int(old_data["times"])

  new_data = {
    "exec_time": ((ex * times) + ctx) / (times + 1) ,
    "times": times + 1,
  }

  print("> [writing ] {} to {} file".format(new_data, out_path))

  with open(out_path, 'w') as f:
    json.dump(new_data, f)

def extrapolate_execution_time(stream):
  lines = list(filter(lambda x: x != '', stream.split('\n')))
  for line in lines:
    if "Game" == line.split()[0]:
      print("> [exe_time] {}".format(line.split()[-2]))
      return float(line.split()[-2])

def perform_test(file_name, param, np):
  output = perform_cmd(build_cmd(file_name, param, np))
  ex = extrapolate_execution_time(output)
  append(file_name, param, np, ex)

def exec_tests(tests):
  i = 0
  for file_name in tests:
    for np in tests[file_name]:
      for param in tests[file_name][np]:
        i += 1
        print("performing test number #{}: {} {} {}".format(i, file_name, np, str(param)))
        perform_test(file_name, param, np)
        print("> [finished] wrote output to '{}' file\n".format(build_out_name(file_name, param, np)))

exec_tests(tests)