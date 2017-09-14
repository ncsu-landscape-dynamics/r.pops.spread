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


#include<string>  //string
#include<sys/socket.h>    //socket
#include<arpa/inet.h> //inet_addr

using std::string;

/**
    TCP Client class
*/
class tcp_client
{
private:
    int sock;
    string address;
    int port;
    struct sockaddr_in server;
    
public:
    tcp_client();
    bool conn(string, int);
    bool send_data(string data);
    string receive(int, int &error);
    void close_socket();
};

