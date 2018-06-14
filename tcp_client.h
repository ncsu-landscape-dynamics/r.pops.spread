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

#ifndef TCP_CLIENT_H
#define TCP_CLIENT_H

#include <string>
#include <sys/socket.h>
#include <arpa/inet.h>

/**
    TCP Client class
*/
class tcp_client
{
private:
    int sock;
    std::string address;
    int port;
    struct sockaddr_in server;

public:
    tcp_client();
    bool conn(std::string address, int);
    bool send_data(std::string data);
    std::string receive(int size, int &error);
    void close_socket();
};

#endif /* TCP_CLIENT_H */
