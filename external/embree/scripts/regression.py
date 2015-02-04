#!/usr/bin/python

## ======================================================================== ##
## Copyright 2009-2012 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

# Embree Regression Test Script
# ===============================

# Windows
# -------

# Prerequisites:
#   Install Python 3.2+
#   Install Visual Studio 2013
#   Install Intel C++ Compiler
#   Check out Embree into <embree_dir>

# Instructions:
#   Open the "Visual Studio x64 Cross Tools Command Prompt (2013)"
#   cd <embree_dir>
#   <python_dir>\python.exe <embree_dir>\scripts\benchmark.py run     windows test_dir
#   <python_dir>\python.exe <embree_dir>\scripts\benchmark.py compile windows test_dir

# Linux and OS X
# --------------

# Prerequisites:
#   Install Python 2.6+
#   Install Intel C++ Compiler
#   Check out Embree into <embree_dir>

# Instructions:
#   Open a shell
#   cd <embree_dir>
#   mkdir TEST
#   ./scripts/benchmark.py run     linux <model_dir> test_dir
#   ./scripts/benchmark.py compile linux <model_dir> test_dir

import sys
import os
import re

dash = '/'

########################## configuration ##########################

#compilers_win = ['V120']
#compilers_win = ['ICC']
compilers_win  = ['V120', 'ICC']
#compilers_win  = ['V100', 'V110', 'V120', 'ICC']
#compilers_unix = ['ICC']
compilers_unix = ['GCC', 'CLANG']
#compilers_unix = ['GCC', 'CLANG', 'ICC']
compilers      = []

#platforms_win  = ['Win32']
#platforms_win  = ['x64']
platforms_win  = ['Win32', 'x64']
platforms_unix = ['x64']
platforms      = []

#builds_win = ['Debug']
builds_win = ['Release']
#builds_win = ['Release', 'Debug']
#builds_unix = ['Debug']
builds_unix = ['Release']
#builds_unix = ['Release', 'Debug']
builds = []

#ISAs_win  = ['SSE2']
ISAs_win  = ['SSE2', 'SSE4.2', 'AVX', 'AVX2']
#ISAs_unix = ['AVX2']
ISAs_unix = ['SSE2', 'SSE4.2', 'AVX', 'AVX2']
ISAs = []

supported_configurations = [
  'V120_Win32_Debug_SSE2',   'V120_Win32_Debug_SSE4.2',   'V120_Win32_Debug_AVX',   'V120_Win32_Debug_AVX2', 
  'V120_Win32_Release_SSE2', 'V120_Win32_Release_SSE4.2', 'V120_Win32_Release_AVX', 'V120_Win32_Release_AVX2', 
  'V120_x64_Debug_SSE2',     'V120_x64_Debug_SSE4.2',     'V120_x64_Debug_AVX',     'V120_x64_Debug_AVX2', 
  'V120_x64_Release_SSE2',   'V120_x64_Release_SSE4.2',   'V120_x64_Release_AVX',   'V120_x64_Release_AVX2', 
  'ICC_Win32_Debug_SSE2',    'ICC_Win32_Debug_SSE4.2',    'ICC_Win32_Debug_AVX',    'ICC_Win32_Debug_AVX2', 
  'ICC_Win32_Release_SSE2',  'ICC_Win32_Release_SSE4.2',  'ICC_Win32_Release_AVX',  'ICC_Win32_Release_AVX2', 
  'ICC_x64_Debug_SSE2',      'ICC_x64_Debug_SSE4.2',      'ICC_x64_Debug_AVX',      'ICC_x64_Debug_AVX2', 
  'ICC_x64_Release_SSE2',    'ICC_x64_Release_SSE4.2',    'ICC_x64_Release_AVX',    'ICC_x64_Release_AVX2', 
  'GCC_x64_Debug_SSE2',      'GCC_x64_Debug_SSE4.2',      'GCC_x64_Debug_AVX',      'GCC_x64_Debug_AVX2', 
  'GCC_x64_Release_SSE2',    'GCC_x64_Release_SSE4.2',    'GCC_x64_Release_AVX',    'GCC_x64_Release_AVX2', 
  'CLANG_x64_Debug_SSE2',    'CLANG_x64_Debug_SSE4.2',    'CLANG_x64_Debug_AVX',    'CLANG_x64_Debug_AVX2',  
  'CLANG_x64_Release_SSE2',  'CLANG_x64_Release_SSE4.2',  'CLANG_x64_Release_AVX',  'CLANG_x64_Release_AVX2',  
  ]


models = {}
models['Win32'] = [ 'conference', 'sponza', 'headlight', 'crown', 'bentley' ]
models['x64'  ] = [ 'conference', 'sponza', 'headlight', 'crown', 'bentley', 'xyz_dragon', 'powerplant' ]

modelDir  = ''
testDir = ''

def configName(OS, compiler, platform, build, isa, tutorial, scene, flags):
  cfg = OS + '_' + compiler + '_' + platform + '_' + build + '_' + isa
  if tutorial != '':
    cfg += '_' + tutorial
  if scene != '':
    cfg += '_' + scene
  if flags != '':
    cfg += '_' + flags
  return cfg

########################## compiling ##########################

def compile(OS,compiler,platform,build,isa):

  base = configName(OS, compiler, platform, build, isa, 'build', '', '')
  logFile = testDir + dash + base + '.log'

  if OS == 'windows':

    full_compiler = compiler
    if (compiler == 'ICC'): full_compiler = '"Intel C++ Compiler XE 14.0" '

    # generate build directory
    if os.path.exists('build'):
      ret = os.system('rm -rf build && mkdir build')
      if ret != 0:
        sys.stdout.write("Cannot delete build folder!")
        return ret
    else:	
      os.system('mkdir build')

    # generate solution files using cmake
    command = 'cmake -L '
    command += ' -G "Visual Studio 12 2013"'
    command += ' -T ' + full_compiler
    command += ' -A ' + platform
    command += ' -D COMPILER=' + compiler
    command += ' -D XEON_ISA=' + isa
    command += ' -D RTCORE_RAY_MASK=OFF'
    command += ' -D RTCORE_BACKFACE_CULLING=OFF'
    command += ' -D RTCORE_INTERSECTION_FILTER=ON'
    command += ' -D RTCORE_BUFFER_STRIDE=ON'
    command += ' -D RTCORE_STAT_COUNTERS=OFF'
    command += ' ..'
    os.system('echo ' + command + ' > ' + logFile)
    ret = os.system('cd build && ' + command + ' >> ../' + logFile)
    if ret != 0: return ret

    # compile Embree
    command =  'msbuild build\embree.sln' + ' /m /nologo /p:Platform=' + platform + ' /p:Configuration=' + build + ' /t:rebuild /verbosity:n' 
    os.system('echo ' + command + ' >> ' + logFile)
    return os.system(command + ' >> ' + logFile)
  
  else:

    if (platform != 'x64'):
      sys.stderr.write('unknown platform: ' + platform + '\n')
      sys.exit(1)

    # compile Embree
    command = 'mkdir -p build && cd build && cmake > /dev/null'
    command += ' -D COMPILER=' + compiler
    command += ' -D CMAKE_BUILD_TYPE='+build
    command += ' -D XEON_ISA=' + isa
    command += ' -D RTCORE_RAY_MASK=OFF'
    command += ' -D RTCORE_BACKFACE_CULLING=OFF'
    command += ' -D RTCORE_INTERSECTION_FILTER=ON'
    command += ' -D RTCORE_BUFFER_STRIDE=ON'
    command += ' -D RTCORE_STAT_COUNTERS=OFF'
    command += ' .. && make clean && make -j 8'
    command += ' &> ../' + logFile
    return os.system(command)

def compileLoop(OS):
    for compiler in compilers:
      for platform in platforms:
        for build in builds:
          for isa in ISAs:
            if (compiler + '_' + platform + '_' + build + '_' + isa) in supported_configurations:
              sys.stdout.write(OS + ' ' + compiler + ' ' + platform + ' ' + build + ' ' + isa)
              sys.stdout.flush()
              ret = compile(OS,compiler,platform,build,isa)
              if ret != 0: sys.stdout.write(" [failed]\n")
              else:        sys.stdout.write(" [passed]\n")

########################## rendering ##########################

def render(OS, compiler, platform, build, isa, tutorial, scene, flags):
  sys.stdout.write("  "+tutorial)
  if scene != '': sys.stdout.write(' '+scene)
  if flags != '': sys.stdout.write(' '+flags)
  sys.stdout.flush()
  base = configName(OS, compiler, platform, build, isa, tutorial, scene, flags)
  logFile = testDir + dash + base + '.log'
  imageFile = testDir + dash + base + '.tga'
  if os.path.exists(logFile):
    sys.stdout.write(" [skipped]\n")
  else:
    if OS == 'windows': command = 'build' + '\\' + build + '\\' + tutorial + ' '
    else:               command = 'build' + '/' + tutorial + ' '
    if tutorial[0:10] == 'tutorial10':
      command += '-i tutorials/tutorial10/' + scene + '.xml '
    elif scene != '':
      command += '-c ' + modelDir + dash + scene + dash + scene + '_regression.ecs '
    if tutorial == 'regression':
      command += '-regressions 2000 '
    if tutorial[0:8] == 'tutorial':
      command += '-rtcore verbose=2 -size 1024 1024 -o ' + imageFile
    command += ' > ' + logFile
    ret = os.system(command)
    if ret == 0: sys.stdout.write(" [passed]\n")
    else       : sys.stdout.write(" [failed]\n")

def processConfiguration(OS, compiler, platform, build, isa, models):
  sys.stdout.write('compiling configuration ' + compiler + ' ' + platform + ' ' + build + ' ' + isa)
  sys.stdout.flush()
  ret = compile(OS,compiler,platform,build,isa)
  if ret != 0: sys.stdout.write(" [failed]\n")
  else:        
    sys.stdout.write(" [passed]\n")
                    
    render(OS, compiler, platform, build, isa, 'verify', '', '')
    render(OS, compiler, platform, build, isa, 'benchmark', '', '')

    render(OS, compiler, platform, build, isa, 'tutorial00', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial01', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial02', '', '')
    for model in models:
      render(OS, compiler, platform, build, isa, 'tutorial03', model, 'static')
      render(OS, compiler, platform, build, isa, 'tutorial03', model, 'dynamic')
      render(OS, compiler, platform, build, isa, 'tutorial03', model, 'high_quality')
      render(OS, compiler, platform, build, isa, 'tutorial03', model, 'robust')
      render(OS, compiler, platform, build, isa, 'tutorial03', model, 'compact')

    render(OS, compiler, platform, build, isa, 'tutorial04', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial05', '', '')
    for model in models:
      render(OS, compiler, platform, build, isa, 'tutorial06', model, '')
    render(OS, compiler, platform, build, isa, 'tutorial07', '', '')

    render(OS, compiler, platform, build, isa, 'tutorial08', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial09', '', '')

    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv0', 'static')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv1', 'dynamic')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv2', 'robust')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv3', 'high_quality')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv4', 'static')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv5', 'dynamic')
    render(OS, compiler, platform, build, isa, 'tutorial10', 'subdiv6', 'robust')
			    
    render(OS, compiler, platform, build, isa, 'tutorial00_ispc', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial01_ispc', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial02_ispc', '', '')
    for model in models:
      render(OS, compiler, platform, build, isa, 'tutorial03_ispc', model, '')
    render(OS, compiler, platform, build, isa, 'tutorial04_ispc', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial05_ispc', '', '')
    for model in models:
      render(OS, compiler, platform, build, isa, 'tutorial06_ispc', model, '')
    render(OS, compiler, platform, build, isa, 'tutorial07_ispc', '', 'static')
    render(OS, compiler, platform, build, isa, 'tutorial07_ispc', '', 'dynamic')
    render(OS, compiler, platform, build, isa, 'tutorial07_ispc', '', 'high_quality')

    render(OS, compiler, platform, build, isa, 'tutorial08_ispc', '', '')
    render(OS, compiler, platform, build, isa, 'tutorial09_ispc', '', '')

    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv0', 'static')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv1', 'dynamic')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv2', 'robust')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv3', 'high_quality')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv4', 'static')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv5', 'dynamic')
    render(OS, compiler, platform, build, isa, 'tutorial10_ispc', 'subdiv6', 'robust')

def renderLoop(OS):
    for compiler in compilers:
      for platform in platforms:
        for build in builds:
          for isa in ISAs:
            if (compiler + '_' + platform + '_' + build + '_' + isa) in supported_configurations:
              processConfiguration(OS, compiler, platform, build, isa, models[platform])

########################## command line parsing ##########################

def printUsage():
  sys.stderr.write('Usage: ' + sys.argv[0] + ' compile <os> <testDir>\n')
  sys.stderr.write('       ' + sys.argv[0] + ' run     <os> <testDir> <modelDir>\n')
  sys.exit(1)

if len(sys.argv) < 3: printUsage()
mode = sys.argv[1]
OS = sys.argv[2]

if OS == 'windows':
  dash = '\\'
  compilers = compilers_win
  platforms = platforms_win
  builds = builds_win
  ISAs = ISAs_win
  modelDir = '%HOMEPATH%\\models\\embree'

else:
  dash = '/'
  compilers = compilers_unix
  platforms = platforms_unix
  builds = builds_unix
  ISAs = ISAs_unix
  modelDir = '~/models/embree'

if mode == 'run':
  if len(sys.argv) < 4: printUsage()
  testDir = sys.argv[3]
  if not os.path.exists(testDir):
    os.system('mkdir '+testDir)
  if len(sys.argv) > 4: 
    modelDir = sys.argv[4]
  renderLoop(OS)
  sys.exit(1)

if mode == 'compile':
  if len(sys.argv) < 4: printUsage()
  testDir = sys.argv[3]
  if not os.path.exists(testDir):
    os.system('mkdir '+testDir)
  compileLoop(OS)
  sys.exit(1)

