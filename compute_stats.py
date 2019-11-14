import subprocess
import time
import json
import os.path
import os

ctx = {}
def read_and_save(filepath):
  path_info = filepath.split("/")[1].split(".")[0].split("_")
  file_name = path_info[1]
  np = int(path_info[2][2:])
  params = tuple(map(lambda x: int(x), path_info[3].split("-")))
  par_st = str(params)

  if par_st not in ctx:
    ctx[par_st] = {}
  if np not in ctx[par_st]:
    ctx[par_st][np] = {}
  if file_name not in ctx[par_st][np]:
    ctx[par_st][np][file_name] = {}
  else:
    print("> [error   ] data already inserted: {}".format(filepath))
    return -1

  print("> [compute ] computing test ({}, {}, {})".format(file_name, np, par_st))

  with open(filepath, 'r') as f:
    data = json.load(f)

  ctx[par_st][np][file_name] = {
    'exec_time': data["exec_time"],
  }

  if "seq" in file_name:
    ctx[par_st][np][file_name] = {
      'exec_time': data["exec_time"],
      'problem_size_each_node': params[0]*params[1],
    }
  elif "opt" in file_name:
    ctx[par_st][np][file_name] = {
      'exec_time': data["exec_time"],
      'speedup' : get_speedup(ctx, file_name, np, params),
      'efficiency': get_eff(ctx, file_name, np, params),
      'comm_weight': get_comm_weight(ctx, file_name, np, params),
      'tot_msg': get_tot_msg(ctx, file_name, np, params),
      'increment': get_increment(ctx, file_name, np, params),
      'problem_size_each_node': get_problem_size(ctx, file_name, np, params),
    }
  else:
    ctx[par_st][np][file_name] = {
      'exec_time': data["exec_time"],
      'speedup' : get_speedup(ctx, file_name, np, params),
      'efficiency': get_eff(ctx, file_name, np, params),
      'comm_weight': get_comm_weight(ctx, file_name, np, params),
      'tot_msg': get_tot_msg(ctx, file_name, np, params),
      'problem_size_each_node': get_problem_size(ctx, file_name, np, params),
    }

def iterate_over_data(rootdir):
  for subdir, dirs, files in os.walk(rootdir):
    files.sort()
    files.reverse()
    for file in files:
      #print os.path.join(subdir, file)
      filepath = subdir + os.sep + file

      if filepath.endswith(".data"):
        read_and_save(filepath)

def get_avg_ex(out_obj, fil, core, params):
  if fil not in out_obj[str(params)][core]:
    return -1
  else:
    return out_obj[str(params)][core][fil]["exec_time"]

def get_speedup(out_obj, fil, core, param):
  if 1 not in out_obj[str(param)] or 'gol-seq' not in out_obj[str(param)][1]:
    return -1
  else:
    seq_time = get_avg_ex(out_obj, 'gol-seq', 1, param)
    this_time = get_avg_ex(out_obj, fil, core, param)
    if this_time > seq_time:
      return -1
    else:
      return seq_time / this_time

def get_eff(out_obj, fil, core, param):
  if 1 not in out_obj[str(param)] or 'gol-seq' not in out_obj[str(param)][1]:
    return -1
  else:
    speedup = get_speedup(out_obj, fil, core, param)
    if speedup == -1:
      return -1
    else:
      return get_speedup(out_obj, fil, core, param) / core

def get_comm_weight(out_obj, fil, core, param):
  if fil == 'gol-par':
    return param[1]
  else:
    return int(param[1] / 31) + 1 # optimized

def get_tot_msg(out_obj, fil, core, param):
  return get_comm_weight(out_obj, fil, core, param) * 2 * param[2] * core

def get_problem_size(out_obj, fil, core, param):
  return (int(param[0] / core) + 1) * param[1]

def get_increment(out_obj, fil, core, param):
  if 'gol-par' not in out_obj[str(param)][core]:
    return -1
  else:
    par = get_avg_ex(out_obj, "gol-par", core, param)
    opt = get_avg_ex(out_obj, fil, core, param)

    return (par-opt)*100/par

iterate_over_data("out")

with open("stats.json", 'w') as f:
  json.dump(ctx, f)