from __future__ import absolute_import
from __future__ import print_function
import util
if util.is_python3():
    import pickle
else:
    import cPickle as pickle

class MSObject(object):
    '''MSObject is the base class of all HCI classes.  This is enable us to dump and load objects.
    Example:
    yap.dump('yap', 'cache') # dump yap HCI object to cache/yap.picke
    new_yap=hci.HCI.load('cache/yap.pickle') # create new_yap HCI object from previous dump
    hci.HCIObject.dump([obj1, obj2], 'mylist', 'cache') # dump a python list'''
    CACHE="."

    def __init__(self):
        pass

    @staticmethod
    def dump_object(obj, s_name='untitled', s_cache_dir="."):
        """use to dump any object
        s_name: str, filename, defaults to 'untitled', filename will be appended with suffix '.pickle'
        s_cache_dir: str, filepath, defaults to current folder '.'
        no return"""
        s_cache_dir=s_cache_dir if s_cache_dir is not None else ""
        import os
        s_file=os.path.join(s_cache_dir, s_name)
        util.dump_object(obj, s_file)

    def dump(self, s_name='untitled', s_cache_dir="."):
        """dump instance inherited from HCIObject
        s_name: str, filename, defaults to 'untitled', filename will be appended with suffix '.pickle'
        s_cache_dir: str, filepath, defaults to current folder '.'
        no return"""
        import sys
        try:
            MSObject.dump_object(self, s_name=s_name, s_cache_dir=s_cache_dir)
        except:
            print("Unexpected error:", sys.exc_info()[0])

    @staticmethod
    def load(s_file):
        return util.load_object(s_file)
