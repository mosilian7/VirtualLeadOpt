class A:
    def __init__(self, size, name=None):
        self.name = name
        self.size = size


class B:
    @staticmethod
    def my_size(sth: A):
        return sth.size

    @staticmethod
    def my_size_2(sth):
        return 2 * sth.size


a = A(5, "alice")
l = ["11", "ss"]

import os
PATH = os.environ.get('SCHRODINGER')
print(PATH)
