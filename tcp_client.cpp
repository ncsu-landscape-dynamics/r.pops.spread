/*
 * SOD model - TCP client
 *
 * Copyright (C) 2017 by the authors.
 *
 * Authors: Anna Petrasova (kratochanna gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

#include "tcp_client.h"

using std::cout;
using std::endl;
using std::string;

tcp_client::tcp_client()
{
    sock = -1;
    port = 0;
    address = "";
}

/**
    Connect to a host on a certain port number
*/
bool tcp_client::conn(string address, int port)
{
    // create socket if it is not already created
    if(sock == -1)
    {
        // create socket
        sock = socket(AF_INET , SOCK_STREAM , 0);
        if (sock == -1)
        {
            perror("Could not create socket");
        }

        cout<<"Socket created\n";
    }
    // else nothing

    // setup address structure
    if(inet_addr(address.c_str()) == -1)
    {
        struct hostent *he;
        struct in_addr **addr_list;

        // resolve the hostname, its not an IP address
        if ((he = gethostbyname(address.c_str())) == NULL)
        {
            // gethostbyname failed
            herror("gethostbyname");
            cout << "Failed to resolve hostname\n";

            return false;
        }

        // cast the h_addr_list to in_addr, since h_addr_list also has
        // the IP address in long format only
        addr_list = (struct in_addr **) he->h_addr_list;

        for(int i = 0; addr_list[i] != NULL; i++)
        {
            server.sin_addr = *addr_list[i];

            cout << address << " resolved to " << inet_ntoa(*addr_list[i]) << endl;

            break;
        }
    }
    else
    {
        // plain IP address
        server.sin_addr.s_addr = inet_addr(address.c_str());
    }

    server.sin_family = AF_INET;
    server.sin_port = htons( port );

    // connect to remote server
    if (connect(sock, (struct sockaddr *)&server, sizeof(server)) < 0)
    {
        perror("connect failed. Error");
        return 1;
    }

    cout << "Connected\n";
    return true;
}

/**
    Send data to the connected host
*/
bool tcp_client::send_data(string data)
{
    // send some data
    if(send(sock, data.c_str(), strlen(data.c_str()), 0) < 0)
    {
        perror("Send failed : ");
        return false;
    }
    cout << "Data send\n";

    return true;
}

/**
    Receive data from the connected host
*/
string tcp_client::receive(int size, int &error)
{
    char buffer[size];
    memset(&buffer[0], 0, sizeof(buffer));
    string reply;

    // receive a reply from the server
    error = recv(sock, buffer, sizeof(buffer), 0);
    reply = buffer;
    return reply;
}

void tcp_client::close_socket() {
    shutdown(sock, 0);
}
