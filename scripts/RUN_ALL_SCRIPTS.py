import glob
from os.path import basename
from traceback import format_exc

paths = sorted(glob.glob('ex*.py'))

for path in paths:
    try:
        exec(open(path).read())
        print('SUCCESS: ', path)
    except Exception as e:
        print('FAILED: ', path)
        print(format_exc())