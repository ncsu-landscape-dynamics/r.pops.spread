#!/usr/bin/env python
import os
import socket
import sys
import ssl
from thread import start_new_thread
from threading import Event
import select

import grass.script as gscript

HOST = ''   # Symbolic name, meaning all available interfaces
PORT = 8889  # Arbitrary non-privileged port
PORT_C = 8000

TMP_DIR = '/tmp/'

PROCESS = None

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s_c = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
print 'Sockets created'

# Bind socket to local host and port
try:
    s.bind((HOST, PORT))
except socket.error as msg:
    print 'Bind failed. Error Code : ' + str(msg[0]) + ' Message ' + msg[1]
    sys.exit()

# Bind socket to local host and port
try:
    s_c.bind((HOST, PORT_C))
except socket.error as msg:
    print 'Bind failed. Error Code : ' + str(msg[0]) + ' Message ' + msg[1]
    sys.exit()

print 'Sockets bind complete'

# Start listening on socket
s.listen(10)
s_c.listen(10)
print 'Socket now listening'
connections = {}


def run_model(settings):
    model = 'r.spread.sod.steering'
    params = {}
    params['species'] = 'UMCA_den_100m@PERMANENT'
    params['lvtree'] = 'TPH_den_100m@PERMANENT'
    params['infected'] = 'init_2000_cnt@PERMANENT'
    params['output'] = 'output'
    params['output_series'] = 'output'
    params['wind'] = 'NE'
    params['start_time'] = 2000
    params['end_time'] = 2010
    params['random_seed'] = 42
    params['ip_address'] = 'localhost'
    params['port'] = 8000
    params['moisture'] = '/home/anna/Documents/Projects/SOD2/sonoma_weather_Mcoef.txt'
    params['temperature'] = '/home/anna/Documents/Projects/SOD2/sonoma_weather_Ccoef.txt'
    params.update(settings)
    print params
    global PROCESS
    PROCESS = gscript.start_command(model, overwrite=True, **params)


def clientGUI(conn, connections, event):
    # Sending message to connected client
    conn.sendall('Welcome to the server.\n')

    # infinite loop so that function do not terminate and thread do not end.
    while True:
        # receiving from client
        data = conn.recv(1024)
        message = data.split(':')
        if message[0] == 'clientfile':
            # receive file
            fsize, path = int(message[1]), message[2]
            conn.sendall(data)
            f = open('/tmp/test_file.py', 'wb')
            data = conn.recv(1024)
            total_received = len(data)
            f.write(data)
            while(total_received < fsize):
                data = conn.recv(1024)
                total_received += len(data)
                f.write(data)
            f.close()
            conn.sendall('{} received: {} bytes'.format(path, os.path.getsize('/tmp/test_file.py')))
            if 'computation' in connections:
                connections['computation'].sendall('load:{}'.format(path))
        if message[0] == 'serverfile':
            fsize, path = int(message[1]), message[2]
            with open(path, 'rb') as f:
                data = f.read()
                try:
                    conn.sendall(data)
                except socket.error:
                    print 'erroro sending file'
                event.set()
        if message[0] == 'cmd':
            if message[1] == 'start':
                params = {}
                if len(message) == 3:  # additional parameters
                    for each in message[2].split('|'):
                        key, val = each.split('=')
                        params[key] = val
                if 'computation' not in connections:
                    run_model(params)
            elif message[1] == 'stop':
                print "server: get stop from GUI"
                if 'computation' in connections:
                    print "server: send stop from GUI to OSD"
                    connections['computation'].sendall('cmd:stop')
                    global PROCESS
                    PROCESS.wait()
                    PROCESS = None
                    connections['computation'].close()
                    del connections['computation']
            elif message[1] == 'play':
                if 'computation' in connections:
                    connections['computation'].sendall('cmd:play')
            elif message[1] == 'pause':
                if 'computation' in connections:
                    connections['computation'].sendall('cmd:pause')
            elif message[1] == 'stepf':
                if 'computation' in connections:
                    connections['computation'].sendall('cmd:stepf')
            elif message[1] == 'stepb':
                if 'computation' in connections:
                    connections['computation'].sendall('cmd:stepb')

        # client closed
        if not data:
            break
    # came out of loop
    conn.shutdown(socket.SHUT_WR)
    conn.close()
    del connections['GUI']


def clientComputation(conn, connections, event):
    # Sending message to connected client
    conn.sendall('Welcome to the server.\n')  # send only takes string
    # this event blocks sending messages to GUI
    # when GUI expects files
    event.set()
    while True:
        event.wait(2000)
        data = conn.recv(200)
        message = data.split('|')
        for m in message:
            lm = m.split(':')
            event.wait(2000)
            if lm[0] == 'output':
                # r.pack
                pack_path = TMP_DIR + lm[1] + '.pack'
                gscript.run_command('r.pack', input=lm[1], output=pack_path, overwrite=True)
                if 'GUI' in connections:
                    event.clear()
                    connections['GUI'].sendall('serverfile:{}:{}'.format(os.path.getsize(pack_path), pack_path))
            elif lm[0] == 'info':
                if lm[1] == 'last':
                    connections['GUI'].sendall(m)

        if not data:
            break

    # came out of loop
    conn.close()

event = Event()
while True:
    read, write, error = select.select([s, s_c], [s, s_c], [])
    for r in read:
        conn, addr = r.accept()
        if r == s:
#            conn = ssl.wrap_socket(conn, server_side=True, cert_reqs=ssl.CERT_REQUIRED,
#                                 ca_certs="/etc/ssl/certs/SOD.crt",
#                                 certfile="/etc/ssl/certs/server.crt",
#                                 keyfile="/etc/ssl/private/server.key")
            connections['GUI'] = conn
            start_new_thread(clientGUI, (conn, connections, event))
        else:
            connections['computation'] = conn
            start_new_thread(clientComputation, (conn, connections, event))
