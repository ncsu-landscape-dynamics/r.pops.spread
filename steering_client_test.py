# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 21:58:31 2018

@author: anna
"""

import wx
import os
import socket
import threading
import Queue

import wx.lib.newevent

import grass.script as gscript

TMP_DIR = '/tmp/fromserver'

updateDisplay, EVT_UPDATE_DISPLAY = wx.lib.newevent.NewEvent()


class SteeringFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent=parent, title=title,  size=(350, 200))
        
        
        self.socket = None
        self.urlSteering = 'localhost:8889'
        self.resultsToDisplay = Queue.Queue()
        

        panel = wx.Panel(self)
        box = wx.BoxSizer(wx.VERTICAL)

        btnStart = wx.Button(panel, label="Start")
        btnStop = wx.Button(panel, label="Stop")
        btnPlay = wx.Button(panel, label="Play")
        btnPause = wx.Button(panel, label="Pause")
        btnForward = wx.Button(panel, label="Step forward")
        btnBack = wx.Button(panel, label="Step back")

        btnStart.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:start'))
        btnStop.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:stop'))
        btnPlay.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:play'))
        btnPause.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:pause'))
        btnForward.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:stepf'))
        btnBack.Bind(wx.EVT_BUTTON, lambda evt: self.socket.sendall('cmd:stepb'))
        
        self.Bind(EVT_UPDATE_DISPLAY, self._update)

        self.infoText = wx.StaticText(panel, label=''*20)
        # btnChangeInput = wx.Button(panel, label="Change input")

        box1 = wx.BoxSizer(wx.HORIZONTAL)
        box1.Add(btnStart)
        box1.Add(btnStop)
        box.Add(box1)

        box2 = wx.BoxSizer(wx.HORIZONTAL)
        box2.Add(btnBack)
        box2.Add(btnPlay)
        box2.Add(btnPause)
        box2.Add(btnForward)
        box.Add(box2)

        box3 = wx.BoxSizer(wx.HORIZONTAL)
        box3.Add(self.infoText)
        box.Add(box3)

        panel.SetSizer(box)
        panel.Layout()
        self._connectSteering()

    def _update(self, event):
        self.infoText.SetLabel(event.value)

    def _connectSteering(self):
        if self.socket:
            return
        urlS = self.urlSteering.split(':')
        self.socket = socket.socket()
#        self.s = ssl.wrap_socket(self.s, cert_reqs=ssl.CERT_REQUIRED,
#                                 certfile="/etc/ssl/certs/SOD.crt",
#                                 keyfile="/etc/ssl/private/ssl-cert-snakeoil.key",
#                                 ca_certs="/etc/ssl/certs/test_certificate.crt")
        try:
            self.socket.connect((urlS[0], int(urlS[1])))
        except socket.error, exc:
            print "Error connecting to steering server: {}".format(exc)
            self.socket = None
            return

        self.isRunningClientThread = True
        self.clientthread = threading.Thread(target=self._client, args=(self.resultsToDisplay, ))
        self.clientthread.start()

    def _client(self, resultsToDisplay):
        while self.isRunningClientThread:
            data = self.socket.recv(1024)
            if not data:
                # GUI received close from server
                # finish while loop
                self.socket.close()
                continue
            message = data.split(':')
            if message[0] == 'clientfile':
                _, fsize, path = message
                with open(message[2], 'rb') as f:
                    data = f.read()
                    try:
                        self.socket.sendall(data)
                    except socket.error:
                        print 'erroro sending file'
            elif message[0] == 'serverfile':
                # receive file
                fsize, path = int(message[1]), message[2]
                self.socket.sendall(data)
                data = self.socket.recv(1024)
                total_received = len(data)
                if not os.path.exists(TMP_DIR):
                    os.mkdir(TMP_DIR)
                new_path = os.path.join(TMP_DIR, os.path.basename(path))
                f = open(new_path, 'wb')
                f.write(data)
                while(total_received < fsize):
                    data = self.socket.recv(1024)
                    total_received += len(data)
                    f.write(data)
                f.close()
                ##########
#                gscript.run_command('r.unpack', input=new_path, overwrite=True, quiet=True)
#                name = os.path.basename(path).strip('.pack')
#                resultsToDisplay.put(name)
                ##########
                    
                gscript.run_command('r.unpack', input=new_path, overwrite=True, quiet=True)
                name = os.path.basename(path).strip('.pack')
                    # avoid showing aggregate result
                    #if len(name.split('_')) == 6 or len(name.split('_')) == 7:
                resultsToDisplay.put(name)
                evt = updateDisplay(value=name)
                wx.PostEvent(self, evt)

                ##########
            elif message[0] == 'info':
                if message[1] == 'last':
                    name = message[2]
                    evt = updateDisplay(value='last: ' + name)
                    wx.PostEvent(self, evt)

    def OnClose(self, event):
        # first set variable to skip out of thread once possible
        self.isRunningClientThread = False
        try:
            # send message to server that we finish sending
            # then we receive empty response, see above
            if self.socket:
                self.socket.shutdown(socket.SHUT_WR)
        except socket.error, e:
            print e
            pass
        # wait for ending the thread
        if self.clientthread and self.clientthread.isAlive():
            self.clientthread.join()
        # allow clean up in main dialog
        event.Skip()


if __name__ == '__main__':
    app = wx.App()
    top = SteeringFrame(parent=None, title="Steering Client")
    top.Show()
    app.MainLoop()
