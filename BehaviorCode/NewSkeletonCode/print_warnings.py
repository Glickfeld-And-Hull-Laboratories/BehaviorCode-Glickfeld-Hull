import sys

mw_python_path = '/Library/Application Support/MWorks/Scripting/Python'
if mw_python_path not in sys.path:
    sys.path.insert(0, mw_python_path)

from mworks.data import MWKFile


def print_warnings(filename):
    with MWKFile(filename) as fp:
        for evt in fp.get_events_iter(codes=['#announceMessage']):
            data = evt.data
            if isinstance(evt.data, dict):
                msg = data.get('message', '')
                if msg.startswith('WARNING: '):
                    print(evt.time, msg)


if __name__ == '__main__':
    print_warnings(sys.argv[1])
