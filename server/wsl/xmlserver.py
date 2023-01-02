# prerequisite: path of gnina is added to PATH


from xmlrpc.server import SimpleXMLRPCServer
from xmlrpc.server import SimpleXMLRPCRequestHandler
from socketserver import ThreadingMixIn
import subprocess
import time
import os
import sys


# Restrict to a particular path.
class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths = ('/MRPC', '/RPC2')


class ThreadXMLRPCServer(ThreadingMixIn, SimpleXMLRPCServer):
    pass


class MyFunctions:
    """
        my functions
    """

    @staticmethod
    def run_args(args):
        """
        Run args with command line.
        :param args: list of arguments
        :return: stdout of the process
        """
        cp = subprocess.run(args, shell=False, capture_output=True, encoding="utf-8", errors="ignore")

        curr_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        if len(cp.stdout) > 0:
            print(f"\033[32m[{curr_time}] {args[0]} STDOUT:\033[0m" + cp.stdout)
        if len(cp.stderr) > 0:
            print(f"\033[31m[{curr_time}] {args[0]} STDERR:\033[0m" + cp.stderr)

        return cp


# Create server
with ThreadXMLRPCServer(('0.0.0.0', 8111),
                        requestHandler=RequestHandler, allow_none=True) as server:
    server.register_introspection_functions()
    server.register_instance(MyFunctions())
    server.serve_forever()
