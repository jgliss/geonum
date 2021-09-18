# -*- coding: utf-8 -*-
#
# Pyplis is a Python library for the analysis of UV SO2 camera data
# Copyright (C) 2017 Jonas Gliß (jonasgliss@gmail.com)
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License a
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import glob
from os.path import basename
from traceback import format_exc
from SETTINGS import OPTPARSE

paths = sorted(glob.glob('ex*.py'))

# init arrays, that store messages that are printed after execution of all 
# scripts
test_err_messages = []
skipped_messages = []
passed_messages = []

for path in paths:
    print(path)
    try:
        exec(open(path).read())
        passed_messages.append("All tests passed in script: %s" %basename(path))
    except AssertionError:
        msg = ("\n"
               "--------------------------------------------------------\n"
               "Tests in script %s failed.\n"
               "Reason:\n %s\n"
               "--------------------------------------------------------"
               "\n"
               %(basename(path), format_exc()))
        test_err_messages.append(msg)
    except ImportError as e:
        msg = ("\n"
               "--------------------------------------------------------\n"
               "Feature missing for executing {}.\n"
               "Reason:\n{}\n"
               "--------------------------------------------------------"
               "\n".format(basename(path), repr(e)))
        skipped_messages.append(msg)
    except Exception as e:
        msg = ("\n"
               "--------------------------------------------------------\n"
               "Error when executing {}.\n"
               "Reason:\n{}\n"
               "--------------------------------------------------------"
               "\n".format(basename(path), repr(e)))
        test_err_messages.append(msg)
(options, args) = OPTPARSE.parse_args()

# If applicable, do some tests. This is done only if TESTMODE is active: 
# testmode can be activated globally (see SETTINGS.py) or can also be 
# activated from the command line when executing the script using the 
# option --test 1
if int(options.test):
    print("\n----------------------------\n"
          "T E S T  F A I L U R E S"
          "\n----------------------------\n")  
    if test_err_messages:   
        for msg in test_err_messages:
            print(msg)
    else:
        print("None")
    print("\n----------------------------\n"
          "T E S T S  S K I P P E D"
          "\n----------------------------\n")  
    if skipped_messages:   
        for msg in skipped_messages:
            print(msg)
    else:
        print("None")
    print("\n----------------------------\n" 
          "T E S T  S U C C E S S" 
          "\n----------------------------\n")
    if passed_messages:
        for msg in passed_messages:
            print(msg)
    else:
        print("None")