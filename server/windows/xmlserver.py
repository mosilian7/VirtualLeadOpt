# python environment: see https://github.com/wengong-jin/hgraph2graph


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
        :param logging: if logging, stdout and stderr will be writen to log
        :param log: opened log file
        :return: stdout of the process
        """
        cp = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")

        curr_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        if len(cp.stdout) > 0:
            print(f"\033[32m[{curr_time}] {args[0]} STDOUT:\033[0m" + cp.stdout)
        if len(cp.stderr) > 0:
            print(f"\033[31m[{curr_time}] {args[0]} STDERR:\033[0m" + cp.stderr)

        return cp


# Create server
with ThreadXMLRPCServer(('localhost', 8112),
                        requestHandler=RequestHandler, allow_none=True) as server:
    server.register_introspection_functions()
    server.register_instance(MyFunctions())
    server.serve_forever()
