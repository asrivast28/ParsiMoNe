##
# @file SConstruct
# @brief Provides general functionality for building.
# @author Ankit Srivastava <asrivast@gatech.edu>
#
# Copyright 2020 Georgia Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
import platform

topDir = os.path.abspath(os.getcwd())
cpp = '' if ARGUMENTS.get('CCACHE', 0) in [0, '0'] else 'ccache '
buildDir = None
targetSuffix = ''


cppPaths = [
            topDir,
            ]

cppFlags = [
            ]

cppDefs = [
           ]

libPaths = [
            ]

allLibs = [
           'boost_filesystem',
            ]

linkFlags = [
             ]

# Flag for building in debug mode. Defaults to release build.
releaseBuild = ARGUMENTS.get('DEBUG', 0) in [0, '0']
# Comma separated list of non-default library locations
localLibPaths = ARGUMENTS.get('LOCALLIBS')
if localLibPaths is not None:
  libPaths.extend(localLibPaths.split())
# Comma separated list of non-default include directory locations
localIncludePaths = ARGUMENTS.get('LOCALINCLUDES')
if localIncludePaths is not None:
  cppPaths.extend(localIncludePaths.split())
# Location of the mxx directory
mxxDir = ARGUMENTS.get('MXX', os.path.join(topDir, 'mxx'))
if os.path.exists(mxxDir):
  cppPaths.append(os.path.join(mxxDir, 'include'))
  cppPaths.append(os.path.join(mxxDir, 'gtest'))
  cppPaths.append(os.path.join(mxxDir, 'ext'))


if platform.system() in ['Darwin', 'Linux']:
  cpp += 'mpic++'
  cppFlags.extend([
              '-Wall',
              '-Wextra',
              '-std=c++14',
              ])
  if releaseBuild:
    cppFlags.append('-O3')
    cppDefs.append('NDEBUG')
  else:
    cppFlags.append('-g')

# For OS X
if platform.system() == 'Darwin':
  libPaths.append('/opt/local/lib')
  cppPaths.append('/opt/local/include')

if not releaseBuild:
  targetSuffix += '_debug'

enableLogging = ARGUMENTS.get('LOGGING', None)
if enableLogging is not None:
  # use the user provided option for enabling/disabling logging
  enableLogging = False if enableLogging in [0, '0'] else True
else:
  # otherwise, enable logging only in debug build
  enableLogging = not releaseBuild
if enableLogging:
  # logging libraries
  cppDefs.append('LOGGING')
  allLibs.extend(['boost_log', 'boost_log_setup', 'boost_system', 'boost_thread', 'pthread'])

enableTimer = ARGUMENTS.get('TIMER', None)
if enableTimer is not None:
  # use the user provided option for enabling/disabling timer
  enableTimer = False if enableTimer in [0, '0'] else True
if enableTimer:
  # logging libraries
  cppDefs.append('TIMER')

# Flag for enabling profiling
profiler = ARGUMENTS.get('PROFILER')
if profiler is not None:
  if platform.system() == 'Linux':
    if profiler == 'gprof':
      cppFlags.append('-pg')
      linkFlags.append('-pg')
    elif profiler == 'hpctoolkit':
      # generate debug symbols
      cppFlags.insert(0, '-g')
    else:
      print('ERROR: Profiler "%s" is not supported' % profiler)
      Exit(1)
    targetSuffix += '_%s' % profiler
  else:
    print('WARNING: Profiling is not supported on', platform.system())

if ARGUMENTS.get('SANITIZE', 0) not in [0, '0']:
  if platform.system() == 'Linux':
    cppFlags.append('-fsanitize=address')
    #cppFlags.append('-fsanitize=leak')
    #cppFlags.append('-fsanitize=undefined')
    cppFlags.append('-fno-omit-frame-pointer')
    linkFlags.append('-fsanitize=address')

suffix = ARGUMENTS.get('SUFFIX')
if suffix is not None:
  targetSuffix += suffix

if ARGUMENTS.get('LOCALENVIRON', 1) not in [0, '0']:
  # Use all the paths specified by the environment
  for path in os.environ.get('CPATH', '').split(os.pathsep):
    cppPaths.append(path)
  for path in os.environ.get('LIBRARY_PATH', '').split(os.pathsep):
    libPaths.append(path)

  # Use all the flags specified by the environment
  for flag in os.environ.get('CPPFLAGS', '').split():
    if flag.startswith('-I'):
      cppPaths.append(flag[2:])
    else:
      cppFlags.append(flag)
  for flag in os.environ.get('LDFLAGS', '').split():
    if flag.startswith('-L'):
      libPaths.append(flag[2:])
    else:
      linkFlags.append(flag)


env = Environment(ENV=os.environ, CXX=cpp, CXXFLAGS=cppFlags, CPPPATH=cppPaths, CPPDEFINES=cppDefs, LIBPATH=libPaths, LINKFLAGS=linkFlags)

def CheckCXXBuilder(context):
  '''
  Custom test for determining if the C++ build environment works.
  '''
  # First check if the compiler works
  result = context.sconf.CheckCXX()
  if result:
    # Now, check if the linker works
    context.Message('Checking whether the C++ linker works...')
    # Use the file generated during the compiler testing
    source = context.sconf.lastTarget.sources[0]
    result = context.TryLink(source.get_contents().decode('utf-8'), source.get_suffix())
    context.Result(result)
  return result
# Add the custom test for testing C++ build environment in the configuration
conf = Configure(env, custom_tests = {'CheckCXXBuilder' : CheckCXXBuilder})

# Check if the initial build environment works
if not conf.CheckCXXBuilder():
  Exit(1)

conf.env.Append(CXXFLAGS='-march=native')
# Check for functions and build options required by the external libraries
extBuiltins = ['ctzll', 'popcountll']
for builtin in extBuiltins:
  if not conf.CheckDeclaration('__builtin_%s' % builtin):
    print('ERROR: __builtin_%s is required by bit_util.hpp' % builtin)
    Exit(1)

# Check for boost header location
if not conf.CheckCXXHeader('boost/version.hpp'):
  Exit(1)
# Check for all the libraries
for lib in allLibs:
  if not conf.CheckLib(lib, language='C++'):
    Exit(1)
# Use the modified environment
env = conf.Finish()

env.targetSuffix = targetSuffix
env.topDir = topDir

if targetSuffix:
  buildDir = targetSuffix
  if buildDir.startswith('_'):
    buildDir = buildDir[1:]
else:
  buildDir = 'release'
buildDir = os.path.join('builds', buildDir)

if not SConscript('SConscript', exports='env', src_dir='.', variant_dir=buildDir, duplicate=0):
  print('Executable was not built')
  Exit(1)
