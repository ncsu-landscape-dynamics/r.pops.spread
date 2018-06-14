#!/usr/bin/env python
import os
import socket
import sys
import ssl
from thread import start_new_thread
from threading import Event
import select
import time

import grass.script as gscript


HOST = ''   # Symbolic name, meaning all available interfaces
PORT = 8888  # Arbitrary non-privileged port
#PORT_C = 8000

TMP_DIR = '/tmp/'

PROCESS = None

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
print 'Socket created'

# Bind socket to local host and port
try:
    s.bind((HOST, PORT))
except socket.error as msg:
    print 'Bind failed. Error Code : ' + str(msg[0]) + ' Message ' + msg[1]
    sys.exit()


print 'Socket bind complete'

# Start listening on socket
s.listen(10)

print 'Socket now listening'

PROCESS = None

COMMANDS = [
    'step',
    'wait 2',
    'step',
    'wait 2',
    'play',
    'wait 2',
    'end'
]


def run_model():
    model = 'r.spread.sod.steering'
    params = {}
    params['output_series'] = 'output'
    params['output'] = 'output'
    params['random_seed'] = 42
    params['nprocs'] = 1
    params['ip_address'] = 'localhost'
    params['port'] = 8888
    params['moisture_file'] = '/home/anna/Documents/Projects/SOD2/chapter/weather_files/moisture_file_future.txt'
    params['temperature_file'] = '/home/anna/Documents/Projects/SOD2/chapter/weather_files/temperature_file_future.txt'
    params['species'] = 'lide_den_int'
    params['lvtree'] = 'all_den_int'
    params['infected'] = 'inf_2016'
    params['start_time'] = '2001'
    params['end_time'] = '2005'
    params['wind'] = 'NE'
    global PROCESS
    PROCESS = gscript.start_command(model, overwrite=True, flags='l', **params)


def run(conn):
    for command in COMMANDS:
        if command.startswith('start'):
            if PROCESS is None:
                run_model()
        elif command.startswith('end'):
            conn.sendall('cmd:stop')
            global PROCESS
            PROCESS.wait()
            PROCESS = None
            conn.close()
            del conn
        elif command.startswith('step'):
            conn.sendall('cmd:step')
        elif command.startswith('play'):
            conn.sendall('cmd:play')
        elif command.startswith('wait'):
            time.sleep(float(command.split()[-1]))


run_model()

while True:
    read, write, error = select.select([s], [s], [])
    for r in read:
        conn, addr = r.accept()
        if r == s:
            start_new_thread(run, (conn, ))

