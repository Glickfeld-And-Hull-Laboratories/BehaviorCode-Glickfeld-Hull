import os
import sys

mw_python_path = '/Library/Application Support/MWorks/Scripting/Python'
if mw_python_path not in sys.path:
    sys.path.insert(0, mw_python_path)

from mworks import file_hash
from mworks.data import MWKFile


def get_image_hashes(filename):
    hashes = {}

    with MWKFile(filename) as fp:
        for evt in fp.get_events_iter(codes=['#stimDisplayUpdate']):
            for info in evt.data:
                if info.get('type') == 'frame_list':
                    info = info['current_stimulus']
                    if isinstance(info, dict):
                        assert info.get('type') == 'image'
                        image_hash = info['file_hash']
                        image_path = info['filename']
                        # If we've seen this hash before, verify that the path
                        # is the same, too
                        assert hashes.get(image_hash, image_path) == image_path
                        hashes[image_hash] = image_path

    return hashes


def check_image_hashes(image_dir, hashes):
    images = []

    for dirpath, dirnames, filenames in os.walk(image_dir):
        for name in filenames:
            image_path = os.path.join(dirpath, name)
            image_hash = file_hash(image_path)

            image_path = os.path.relpath(image_path, image_dir)
            images.append(image_path)

            if image_hash not in hashes:
                print('No matching hash for image', image_path)
            elif not hashes[image_hash].endswith(image_path):
                print('Name mismatch for image', image_path)

    if len(images) != len(hashes):
        print('Wrong number of image files')


if __name__ == '__main__':
    hashes = get_image_hashes(sys.argv[1])
    check_image_hashes(sys.argv[2], hashes)
